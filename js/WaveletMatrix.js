import { BitVector } from './BitVector.js';
import { RLEBitVector } from './RLEBitVector.js';
import { reverseBits, trailing0, popcount, binarySearchAfter, binarySearchBefore } from './util.js';
import { ScratchSpace } from './ScratchSpace.js';
import { ArrayWalker } from './ArrayWalker.js';

// todo
// - perform index size checks when allocating scratch space, and use float 64 buffers when last >= 2**32
// - enforce len(multiplicity) == len(data)
// - enforce positive multiplicity
// - check scratch space does not get overlarge due to multiplicity
// - consider removing multiplicityIndex and cumulativeMultiplicity and doing that all on the outside,
//   since there may be uses for the run-length encoding that do not involve queries precisely on run boundaries
// - check that all symbols are < alphabetSize
// - check that multiplicity is a typed array (has to be, for arraywalker reset)
// - idea: do not use bitwise operations if we can avoid them inside eg. counts; that way can scale to values up to 2^53!
//   - might not be possible (eg. we do | and &); could try a bigint64...
// - implement SparseBitVector?
// -  explore the idea of storing the complement whenever 1 density exceeds 50%; then rank0 is rank1 and same for select.
// error if ignoreBits > numLevels
// later
// - implement range_next_value, range_intersect, and fingered range quantile from
//   the paper "New algorithms on wavelet trees and applications to information retrieval".
// - use flatqueue; keep a single queue around and use .shrink() between calls for top frequent
// - think about what an API might look like that allows us to specify first and last in batches
// api design
// - mistakenly used count instead of countSymbol; really wanted a rank; maybe rename countSymbol to rank and countSymbolBatch to rankBatch??
// - document that the top of the highest subcode need not be marked; this is why, as a special case, the default value 0 works.
// - make first/last also kwargs?
// - we sometimes accept last == first, but other times not. I think it would be good to support empty ranges if we can; but we need to verify the results are correct.
// Implements a binary wavelet matrix that splits on power-of-two alphabet
// boundaries, rather than splitting based on the true alphabet midpoint.
export class WaveletMatrix {
  // This implements Algorithm 1 (seq.pc) from the paper "Practical Wavelet Tree Construction".
  // The return value is the array of wavelet matrix levels. Adapting the algorithm to construct
  // a wavelet tree instead requires changing the borders computation (see section 5.3).
  // todo: check that all symbols are < alphabetSize
  // todo: pass in maxSymbol, with alphabetSize = maxSymbol + 1?
  constructor(data /* call this symbols? */, maxSymbol, opts = {}) {
    if (maxSymbol > 2 ** 32 - 1) {
      throw new Error('cannot represent integers greater than 2^32 - 1 at the moment');
    }
    // data is an array of integer values in [0, alphabetSize)
    // maxSymbol is the maximum *actual* symbol in this WM;
    // not the maximum *possible* symbol, which is 2**numLevels - 1
    this.maxSymbol = maxSymbol;
    this.alphabetSize = maxSymbol + 1;
    this.numLevels = Math.ceil(Math.log2(this.alphabetSize));
    this.maxLevel = this.numLevels - 1;
    // this.allOnes = 2 ** this.numLevels - 1; // was thinking this could be useful as a selector...
    const { largeAlphabet = 2 ** this.numLevels > data.length, counts } = opts;
    // The more efficient construction algorithm does not scale well to large alphabets,
    // and it cannot handle multiplicities because it constructs the bitvectors out-of-order.
    // It also requires O(2^numLevels) space. So, if conditions are unfavorable, use the
    // more straightforward construction algorithm that fills in levels from top-to-bottom
    // by iterating and incrementally sorting the data.
    // Each of the construction algorithms return
    if (data.length === 0) {
      this.levels = [];
      this.length = 0;
    } else if (counts) {
      this.levels = this.constructLargeAlphabetWithMultiplicity(data, counts);
    } else if (largeAlphabet || data.length >= 2 ** 32) {
      this.levels = this.constructLargeAlphabet(data);
    } else {
      this.levels = this.construct(data);
    }

    // Compute the number of zeros at each level
    const { numLevels, levels } = this;
    const numZeros = new Uint32Array(numLevels);
    for (let i = 0; i < numLevels; i++) {
      numZeros[i] = levels[i].rank0(levels[i].length);
    }
    this.numZeros = numZeros;
    // Initialize a scratch space to hold intermediate computations
    // which lets us reduce the number of allocations during use.
    this.scratch = new ScratchSpace();
    // The length of the wavelet matrix is the same as the
    // length of its bitvectors, all of which are the same
    // length since they represent bits of the same sequence.
    this.length = levels[0].length;
  }

  // Implements Algorithm 1 (seq.pc) from the paper "Practical Wavelet Tree Construction".
  construct(data) {
    console.assert(data.length < 2 ** 32);
    const { numLevels, maxLevel } = this;
    const hist = new Uint32Array(2 ** numLevels);
    const borders = new Uint32Array(2 ** numLevels);

    const levels = new Array(numLevels);
    // Initialize the level bit vectors
    for (let i = 0; i < numLevels; i++) {
      levels[i] = new BitVector(data.length);
    }

    // Compute the histogram of the data
    const level = levels[0];
    const levelBitMask = 1 << maxLevel;
    for (let i = 0; i < data.length; i++) {
      const d = data[i];
      hist[d] += 1;
      // Fill the first bitvector (MSBs in data order)
      if (d & levelBitMask) level.one(i);
    }
    level.finish();

    // Construct the other levels bottom-up
    for (let l = maxLevel; l > 0; l--) {
      const m = 2 ** l;

      // Compute the histogram based on the previous level's one
      for (let i = 0; i < m; i++) {
        // Update the histogram in-place
        hist[i] = hist[2 * i] + hist[2 * i + 1];
      }

      // Get starting positions of intervals from the new histogram
      borders[0] = 0;
      for (let i = 1; i < m; i++) {
        // Update the positions in-place
        const prevIndex = reverseBits(i - 1, l);
        borders[reverseBits(i, l)] = borders[prevIndex] + hist[prevIndex];
      }

      // Fill the bit vector of the current level
      const level = levels[l];
      const levelBit = maxLevel - l;
      const levelBitMask = 1 << levelBit;
      const bitPrefixMask = 0xfffffffe << levelBit;
      const bitPrefixShift = levelBit + 1;
      for (let i = 0; i < data.length; i++) {
        const d = data[i];
        // Get and update position for bit by computing its bit prefix,
        // which encodes the path from the root to the node at level l
        // containing this bit
        // Note: This shift operation is why we can't handle arrays with
        // more than 2^32 elements using this approach.
        const nodeIndex = (d & bitPrefixMask) >>> bitPrefixShift;
        const p = borders[nodeIndex];
        borders[nodeIndex] += 1;
        // Set the bit in the bitvector
        if (d & levelBitMask) level.one(p);
      }
      level.finish();
    }
    return levels;
  }

  // Alternative construction algorithm for the 'sparse' case when the alphabet size
  // is larger than the number of symbols that actually occur in the data.
  constructLargeAlphabet(data) {
    const { numLevels, maxLevel } = this;
    const len = data.length;

    data = new Uint32Array(data);
    let nextData = new Uint32Array(len);
    let tmp; // used for swapping current/next values

    // Initialize the level bit vectors
    const levels = new Array(numLevels);
    for (let i = 0; i < numLevels; i++) {
      levels[i] = new BitVector(len);
    }
    const walk = new ArrayWalker(0, len);
    // For each level, sort the data point by its bit value at that level.
    // Zero bits get sorted left, one bits get sorted right. This amounts
    // to a bucket sort with two buckets.
    // We sort into `nextData`, then swap `nextData` and `data`.
    for (let l = 0; l < maxLevel; l++) {
      const level = levels[l];
      const levelBit = maxLevel - l;
      const levelBitMask = 1 << levelBit;
      for (let i = 0; i < len; i++) {
        const d = data[i];
        if (d & levelBitMask) {
          const ni = walk.nextBackIndex();
          nextData[ni] = d;
          level.one(i);
        } else {
          const ni = walk.nextFrontIndex();
          nextData[ni] = d;
        }
      }
      walk.reset(true, nextData);
      // swap data and nextData
      (tmp = data), (data = nextData), (nextData = tmp);
      level.finish();
    }

    // For the last level we don't need to build anything but the bitvector
    const level = levels[maxLevel];
    const levelBitMask = 1 << 0;
    for (let i = 0; i < len; i++) {
      if (data[i] & levelBitMask) level.one(i);
    }
    level.finish();
    return levels;
  }

  // Extended version of the large-alphabet algorithm supporting element multiplicities.
  constructLargeAlphabetWithMultiplicity(data, counts) {
    const { numLevels, maxLevel } = this;
    const len = data.length;

    data = new Uint32Array(data);
    counts = new Uint32Array(counts);
    let nextData = new Uint32Array(len);
    let nextCounts = new Uint32Array(len);
    let tmp; // used for swapping current/next values

    // Initialize the level bit vectors
    const levels = new Array(numLevels);
    for (let i = 0; i < numLevels; i++) {
      levels[i] = new RLEBitVector();
    }
    const walk = new ArrayWalker(0, len);
    // For each level, sort the data point by its bit value at that level.
    // Zero bits get sorted left, one bits get sorted right. This amounts
    // to a bucket sort with two buckets.
    // We sort into `nextData`, then swap `nextData` and `data`.
    for (let l = 0; l < maxLevel; l++) {
      const level = levels[l];
      const levelBit = maxLevel - l;
      const levelBitMask = 1 << levelBit;
      for (let i = 0; i < len; i++) {
        const d = data[i];
        const m = counts[i];
        if (d & levelBitMask) {
          const ni = walk.nextBackIndex();
          nextData[ni] = d;
          nextCounts[ni] = m;
          level.run(0, m);
        } else {
          const ni = walk.nextFrontIndex();
          nextData[ni] = d;
          nextCounts[ni] = m;
          level.run(m, 0);
        }
      }
      walk.reset(true, nextData, nextCounts);
      // swap data and nextData
      (tmp = data), (data = nextData), (nextData = tmp);
      // swap counts and nextCounts
      (tmp = counts), (counts = nextCounts), (nextCounts = tmp);
      // finish levels as we're done with them, since this can sometimes
      // give us memory savings (eg. with the optimizations in RLEBitVector
      // for vectors with all length-1 one-runs)
      level.finish();
    }

    // For the last level we don't need to build anything but the bitvector
    const level = levels[maxLevel];
    const levelBitMask = 1 << 0;
    for (let i = 0; i < len; i++) {
      if (data[i] & levelBitMask) {
        level.run(0, counts[i]);
      } else {
        level.run(counts[i], 0);
      }
    }
    level.finish();
    return levels;
  }

  access(index) {
    if (typeof index !== 'number') throw new Error('index must be a number');
    if (index < 0 || index > this.length) throw new Error('access: out of bounds');
    let symbol = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      // todo: combine rank and access queries since they likely access the same block
      const index1 = level.rank1(index - 1);
      if (level.access(index) === 0) {
        // go left
        index = index - index1; // = index0
      } else {
        // go right
        const nz = this.numZeros[l];
        index = nz + index1;
        // update symbol
        const levelBitMask = 1 << (this.maxLevel - l);
        symbol |= levelBitMask;
      }
    }
    return symbol;
  }

  accessBatch(indices) {
    for (let i = 0; i < indices.length; i++) {
      const index = indices[i];
      if (typeof index !== 'number') throw new Error('index must be a number');
      if (index < 0 || index > this.length) throw new Error('access: out of bounds');
    }
    // indices get mutated as we go, so make a copy
    indices = new Uint32Array(indices);
    // allocate space for the output symbols
    const symbols = new Uint32Array(indices.length);
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      for (let i = 0; i < indices.length; i++) {
        const index = indices[i];
        const symbol = symbols[i];
        // todo: combine rank and access queries since they likely access the same block
        const index1 = level.rank1(index - 1);
        if (level.access(index) === 0) {
          // go left
          indices[i] = index - index1; // = index0
        } else {
          // go right
          const nz = this.numZeros[l];
          indices[i] = nz + index1;
          // update symbol
          const levelBitMask = 1 << (this.maxLevel - l);
          symbols[i] |= levelBitMask;
        }
      }
    }
    return symbols;
  }

  // Returns the number of occurrences of `symbol` in the range [first, last).
  rank(symbol, { first = 0, last = this.length, ignoreBits = 0 } = {}) {
    const indices = this.symbolRange(symbol, { first, last, ignoreBits });
    return indices.last - indices.first;
  }

  // Returns the index in this wavelet matrix of the nth occurrence of `symbol` in the range [first, last).
  select(symbol, n = 1, { first = 0, last = this.length, ignoreBits = 0 } = {}) {
    if (symbol < 0 || symbol >= this.alphabetSize) return -1;
    if (n < 1 || n > this.length) return -1;
    const indices = this.symbolRange(symbol, { first, last, ignoreBits });
    if (indices.last - indices.first < n) return -1; // in analogy with select
    let index = indices.first + n - 1;
    for (let l = this.numLevels; l > 0; ) {
      l -= 1;
      const level = this.levels[l];
      const nz = this.numZeros[l];
      if (index < nz) {
        // this position was mapped from a zero at the previous level
        const n = index + 1; // this was the nth zero on the that level
        index = level.select0(n); // locate the corresponding index on the preceding level
      } else {
        // this position was mapped from a one at the previous level
        const n = index - nz + 1; // this was the nth one on the that level
        index = level.select1(n); // locate the corresponding index on the preceding level
      }
    }
    return index;
  }

  // Internal function returning the (first, last] index range covered by this symbol on the virtual bottom level,
  // or on a higher level if ignoreBits > 0. `last - first` gives the symbol count within the provided range.
  // If the symbol does not appear in the range, an arbitrary empty range will be returned.
  symbolRange(symbol, { first = 0, last = this.length, ignoreBits = 0 } = {}) {
    const symbolGroupSize = 1 << ignoreBits;
    // actually, this is useful when the bottom bits represent something, ie. the code is a concatenation of subcods
    // if (symbol % symbolGroupSize !== 0) {
    //   // note: could be done with bit math (check that low bits are zero)
    //   throw new Error('symbol must evenly divide the block size implied by ignoreBits');
    // }
    if (symbol >= this.alphabetSize) throw new Error('symbol must be < alphabetSize');
    if (first >= last) return { first, last: first };
    if (first < 0) throw new Error('first must be >= 0');
    if (last > this.length) throw new Error('last must be <= length');
    if (this.numLevels > 0 && first - last > 0) return { first, last };

    const numLevels = this.numLevels - ignoreBits;
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const first1 = level.rank1(first - 1);
      const last1 = level.rank1(last - 1);
      const levelBitMask = 1 << (this.maxLevel - l);

      if ((symbol & levelBitMask) === 0) {
        // go left
        first = first - first1; // = first0
        last = last - last1; // = last0
      } else {
        // go right
        const nz = this.numZeros[l];
        first = nz + first1;
        last = nz + last1;
      }
    }
    return { first, last };
  }

  // Returns the number of values with symbol strictly less than the given symbol.
  // This is kind of like a ranged rank operation over a symbol range.
  // Adapted from https://github.com/noshi91/Library/blob/0db552066eaf8655e0f3a4ae523dbf8c9af5299a/data_structure/wavelet_matrix.cpp#L25
  // We could make a version that returns individual symbol counts for all symbols
  // less than `symbol`; could be useful as an indirect way of determing the number
  // of symbols preceding `symbol`.
  countLessThan(symbol, { first = 0, last = this.length } = {}) {
    if (first >= last) return 0;
    if (first < 0) throw new Error('first must be >= 0');
    if (last > this.length) throw new Error('last must be <= length');
    if (symbol <= 0) return 0;
    if (symbol >= this.alphabetSize) return last - first;
    let count = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const first1 = level.rank1(first - 1);
      const last1 = level.rank1(last - 1); // inclusive
      const levelBitMask = 1 << (this.maxLevel - l);
      if ((symbol & levelBitMask) === 0) {
        // go left
        first = first - first1; // = first0
        last = last - last1; // = last0
      } else {
        // update count before going right
        count += last - last1 - (first - first1); // = last0 - first0
        // go right
        const nz = this.numZeros[l];
        first = nz + first1;
        last = nz + last1;
      }
    }
    return count;
  }

  // Returns the number of occurrences of symbols [lower, upper)
  // in the index range [first, last). It is possible to implement
  // this function roughly twice as efficiently in terms of number
  // of rank calls but at the cost of increased implementation complexity.
  // See the paper "New algorithms on wavelet trees and applications to
  // information retrieval" for details.
  // 🌶 NOTE: Upper is currently exclusive, since that's how countLessThan works...
  count({ first = 0, last = this.length, lower = 0, upper = this.alphabetSize } = {}) {
    // todo: check and error if options or separator or sort are specified
    return this.countLessThan(upper, { first, last }) - this.countLessThan(lower, { first, last });
  }

  rankBatch(sortedSymbols, { first = 0, last = this.length, ignoreBits = 0 } = {}) {
    // todo: error if sortedsymbols.length >= 2^32 (binary search)
    // todo: check and error if separator or sort are specified
    for (let i = 1; i < sortedSymbols.length; i++) {
      if (!(sortedSymbols[i - 1] <= sortedSymbols[i])) throw new Error('sortedSymbols must be sorted');
    }
    const symbolGroupSize = 1 << ignoreBits;
    // for (const symbol of sortedSymbols) {
    //   if (symbol % symbolGroupSize !== 0)
    //     // note: could be done with bit math (check that low bits are zero)
    //     throw new Error('symbol must evenly divide the block size implied by ignoreBits');
    // }
    if (first >= last) return this.emptyResult();
    if (last > this.length) throw new Error('last must be <= wavelet matrix length');

    // account for duplicate sorted symbols and the fact that there cannot be more outputs than elements
    const scratchLength = Math.min(this.alphabetSize, sortedSymbols.length, last - first);
    this.scratch.reset();
    const F = this.scratch.allocU32(scratchLength); // firsts
    const L = this.scratch.allocU32(scratchLength); // lasts
    const S = this.scratch.allocU32(scratchLength); // symbols
    const C = this.scratch.allocU32(scratchLength); // counts
    const walk = new ArrayWalker(sortedSymbols.length === 0 ? 0 : 1, scratchLength);
    const nextIndex = walk.nextFrontIndex();
    F[nextIndex] = first;
    L[nextIndex] = last;
    C[nextIndex] = sortedSymbols.length;
    S[nextIndex] = 0;
    walk.reset(false, F, L, C, S);
    let nRankCalls = 0;
    const numLevels = this.numLevels - ignoreBits;
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBitMask = 1 << (this.maxLevel - l); // determines whether a symbol goes to the left or right
      const symbolsPerNode = (levelBitMask << 1) - 1;
      const symbolBitMask = 0xffffffff << (this.maxLevel - l); // clears the low bits from a symbol, giving the left edge of its node

      let k = sortedSymbols.length; // march k over the sorted symbols array
      for (let i = walk.length; i > 0; ) {
        i--;
        const first = F[i];
        const first1 = level.rank1(first - 1);
        const first0 = first - first1;

        const last = L[i];
        const last1 = level.rank1(last - 1);
        const last0 = last - last1;
        nRankCalls += 2;

        const symbolCount = C[i];
        const symbol = S[i];
        const a = symbol & symbolBitMask; // leftmost symbol in this node
        const b = a + symbolsPerNode; // rightmost symbol in this node
        const m = (a + b) >>> 1;

        // perform two binary searches over the sorted symbols for this node to determine
        // the number in the left child node. We do a single linear pass across sortedSymbols
        // over the course of the walk.
        // Note that this relies on the fact that the symbols in `S` are in sorted order,
        // since we're marching `k` across the sorted symbols array to reduce the search
        // space for the binary search calls to those symbols corresponding to this node.
        // This means that we have to use `walk.nextBackIndex()` for both left and right children,
        // since that's what gives rise to the sorted order of `S` its sorted order.
        // [lo, hi) is the range of sorted offsets covered by this node
        const lo = k - symbolCount;
        const hi = k;
        const splitIndex = binarySearchAfter(sortedSymbols, m, lo, hi);
        const numGoLeft = splitIndex - lo;
        const numGoRight = symbolCount - numGoLeft;
        k = lo;

        if (numGoRight > 0) {
          // go right
          const nextIndex = walk.nextBackIndex();
          F[nextIndex] = nz + first1;
          L[nextIndex] = nz + last1;
          C[nextIndex] = numGoRight;
          S[nextIndex] = symbol | levelBitMask;
        }

        if (numGoLeft > 0) {
          // go left
          const nextIndex = walk.nextBackIndex();
          F[nextIndex] = first0;
          L[nextIndex] = last0;
          C[nextIndex] = numGoLeft;
          S[nextIndex] = symbol;
        }
      }
      walk.reset(false, F, L, C, S);
    }

    for (let i = 0; i < walk.length; i++) L[i] -= F[i];
    const counts = L.subarray(0, walk.length).slice();
    const symbols = S.subarray(0, walk.length).slice();
    return { symbols, counts, nRankCalls };
  }

  // Returns all  distinct symbols [lower, upper] (inclusive) in the range [first, last)
  // together with their number of occurrences. Symbols are grouped and processed
  // in groups of size 2^ignoreBits (symbols are grouped together when they
  // differ only in their lowest `ignoreBits` bits)
  // Each distinct group is labeled by its lowest element, which represents
  // the group containing symbols in the range [symbol, symbol + 2^ignoreBits).
  // todo: call this ranks??
  counts({
    first = 0,
    last = this.length,
    lower = 0,
    upper = this.maxSymbol,
    ignoreBits = 0,
    subcodeSeparator = 0,
    sort = false,
  } = {}) {
    // todo: validate lower/upper bounds wrt alphabet size
    const symbolGroupSize = 1 << ignoreBits;
    // todo: handle lower === upper
    // these error messages could be improved, explaining that ignore bits tells us the power of two
    // that lower and upper need to be multiples of.
    // if (lower % symbolGroupSize !== 0)
    //   throw new Error('lower must evenly divide the symbol block size implied by ignoreBits');
    // if ((upper + 1) % symbolGroupSize !== 0)
    //   throw new Error('(upper + 1) must evenly divide the symbol block size implied by ignoreBits');
    if (first >= last) return this.emptyResult();
    if (first < 0) throw new Error('first must be >= 0');
    if (last > this.length) throw new Error('last must be <= length');

    // if (upper >= this.alphabetSize) throw new Error('upper must be < alphabetSize ([lower, upper] is inclusive)');
    // ^ we now allow this so that subcode stuff works without us having to be unrealistic about the true alphabet size
    // (eg. allow querying code consisting of all maximum subcodes)
    const numLevels = this.numLevels - ignoreBits;

    // todo: bound this more closely. slightly involved to upper-bound due to subcodes;
    // need to compute the product of the subcode ranges since that's the maximum possible
    // number of unique symbols.
    const scratchLength = Math.min(2 ** this.numLevels - ignoreBits, upper - lower + 1, this.length);
    this.scratch.reset();
    const F = this.scratch.allocU32(scratchLength); // firsts
    const L = this.scratch.allocU32(scratchLength); // lasts
    const S = this.scratch.allocU32(scratchLength); // symbols
    const walk = new ArrayWalker(1, scratchLength);
    const reverse = !sort; // for walk.reset(reverse, ...)
    const nextIndex = walk.nextFrontIndex();
    F[nextIndex] = first;
    L[nextIndex] = last;
    S[nextIndex] = 0;
    walk.reset(reverse, F, L, S);

    let nRankCalls = 0;
    // start with all 'extra' high bits set so that we properly handle an
    // `upper` value above ceil(log2(alphabetSize))
    let subcodeMask = 0xffffffff << numLevels;
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBit = this.maxLevel - l;
      const levelBitMask = 1 << levelBit;
      // Usually, the entire code is treated as a single integer, and the [lower, upper] range
      // limits the range of returned codes.
      // It can be useful to instead treat the code as representing a concatenation of subcodes,
      // and the [lower, upper] values as representing a concatenation of the ranges of those
      // subcodes. This behavior can be specified by the `subcodeSeparator` argument, which is a
      // bitmask in which a 1 bit indicates the onset of a new subcode and a 0 implies the continuation
      // of the current subcode. All range comparisons are done within a subcode, and the default
      // subcodeSeparator of 0 gives us the default behavior in which the full code is treated as
      // a single subcode.
      if ((subcodeSeparator & levelBitMask) === 0) subcodeMask |= levelBitMask;
      else subcodeMask = levelBitMask;
      const subcodeLower = lower & subcodeMask;
      const subcodeUpper = upper & subcodeMask;

      // if we want sorted outputs, iterate in reverse to ensure that we don't
      // overwrite unprocessed elements when writing from the back of the array.
      const start = sort ? walk.length - 1 : 0;
      const step = sort ? -1 : 1;
      const end = sort ? -1 : walk.length;
      for (let i = start; i != end; i += step) {
        const first = F[i];
        const first1 = level.rank1(first - 1);
        const first0 = first - first1;

        const last = L[i];
        const last1 = level.rank1(last - 1);
        const last0 = last - last1;

        const symbol = S[i];
        nRankCalls += 2;

        const rightCount = last1 - first1;
        if (rightCount > 0) {
          // go right if the right node range [a, b] overlaps [lower, upper]
          const a = (symbol | levelBitMask) & subcodeMask;
          const b = (a | (levelBitMask - 1)) & subcodeMask;
          if (intervalsOverlapInclusive(a, b, subcodeLower, subcodeUpper)) {
            const nextIndex = walk.nextBackIndex();
            F[nextIndex] = nz + first1;
            L[nextIndex] = nz + last1;
            S[nextIndex] = symbol | levelBitMask;
          }
        }

        const leftCount = last0 - first0;
        if (leftCount > 0) {
          // go left if the left node range [a, b] overlaps [lower, upper]
          // or if there is no range mask. We use an inclusive symbol range since we allow permit all
          // bit patterns as codes, including the maximum value.
          const a = symbol & subcodeMask;
          const b = (a | (levelBitMask - 1)) & subcodeMask;
          if (intervalsOverlapInclusive(a, b, subcodeLower, subcodeUpper)) {
            const nextIndex = sort ? walk.nextBackIndex() : walk.nextFrontIndex();
            F[nextIndex] = first0;
            L[nextIndex] = last0;
            S[nextIndex] = symbol;
          }
        }
      }
      walk.reset(reverse, F, L, S);
    }

    for (let i = 0; i < walk.length; i++) L[i] -= F[i];
    const counts = L.subarray(0, walk.length).slice();
    const symbols = S.subarray(0, walk.length).slice();
    return { symbols, counts, nRankCalls };
  }

  emptyResult() {
    return { symbols: new Uint32Array(), counts: new Uint32Array(), nRankCalls: 0 };
  }

  // Adapted from https://github.com/noshi91/Library/blob/0db552066eaf8655e0f3a4ae523dbf8c9af5299a/data_structure/wavelet_matrix.cpp#L76
  // Range quantile query returning the kth largest symbol in A[i, j).
  // offset is a "sorted index" within the provided index range. If we take all of
  // the elements in [first, last) and sort them, the minimum element lies at offset 0
  // and the maximum element lies at offset last - first - 1.
  quantile(offset, { first = 0, last = this.length } = {}) {
    if (first >= last) return this.emptyResult();
    if (last > this.length) throw new Error('last must be <= wavelet matrix length');
    if (offset < 0 || offset >= last - first)
      throw new Error('offset cannot be less than zero or exceed length of range [first, last)');
    let symbol = 0;
    let nRankCalls = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const first0 = level.rank0(first - 1);
      const last0 = level.rank0(last - 1);
      const count = last0 - first0;
      nRankCalls += 2;
      if (offset < count) {
        // go left
        first = first0;
        last = last0;
      } else {
        // go right
        const nz = this.numZeros[l];
        first = nz + (first - first0); // = nz + first1
        last = nz + (last - last0); // = nz + last1
        // update symbol and new target offset in the child node
        const levelBitMask = 1 << (this.maxLevel - l);
        symbol |= levelBitMask;
        offset -= count;
      }
    }
    return { symbol, count: last - first, nRankCalls };
  }

  // todo: make a minMax function that is like quantile/batchQuantile, but finds specifically
  // the first and last sorted index. right now quantile batch is slow, so this can be used to
  // speed up the enumeration of all symbols in a range together with their counts. though when
  // i put it that way, hold up... why not counts?

  // todo: investigate performance of this and other batch functions.
  // quickly comparing quantileBatch([index]) to quantile(index) seemed to show a possile >2x overhead.
  quantileBatch(sortedOffsets, { first = 0, last = this.length } = {}) {
    // todo: error if sortedoffsets.length >= 2^32 (binary search)
    // these error messages could be improved, explaining that ignore bits tells us the power of two
    // that lower and upper need to be multiples of.
    if (first >= last) return Object.assign(this.emptyResult(), { numSortedOffsets: new Uint32Array() });
    if (last > this.length) throw new Error('last must be <= wavelet matrix length');
    if (sortedOffsets[0] < 0 || sortedOffsets[sortedOffsets.length - 1] >= last - first)
      throw new Error('sorted offsets out of range for [first, last); should be in the range [0, first - last)');
    for (let i = 1; i < sortedOffsets.length; i++) {
      if (!(sortedOffsets[i - 1] <= sortedOffsets[i])) throw new Error('sortedOffsets must be sorted');
    }
    // todo: error if there are more sortedOffsets than last-first, since we copy them into a scratch space (or ensure the space can hold the size we need)
    if (sortedOffsets[0] < 0 || sortedOffsets[sortedOffsets.length - 1] >= last)
      throw new Error('sortedIndex cannot be less than zero or exceed length of range [first, last)');

    // account for duplicate sorted offsets and the fact that there cannot be more outputs than elements
    const scratchLength = Math.min(this.alphabetSize, sortedOffsets.length, last - first);
    this.scratch.reset();
    const F = this.scratch.allocU32(scratchLength); // firsts
    const L = this.scratch.allocU32(scratchLength); // lasts
    const S = this.scratch.allocU32(scratchLength); // symbols
    const C = this.scratch.allocU32(scratchLength); // counts
    const I = this.scratch.allocU32(sortedOffsets.length); // sorted offsets
    const walk = new ArrayWalker(sortedOffsets.length === 0 ? 0 : 1, scratchLength);
    const nextIndex = walk.nextFrontIndex();
    F[nextIndex] = first;
    L[nextIndex] = last;
    S[nextIndex] = 0;
    C[nextIndex] = sortedOffsets.length; // number of sortedOffsets represented by node 0
    walk.reset(true, F, L, S, C);
    // copy sorted offsets into a scratch space since they are mutated as we go
    I.set(sortedOffsets);

    let nRankCalls = 0;

    // note: grouping by LSB computes approximate quantiles where
    // the count of symbols assigned to each range is given by I,
    // and I think the ranges are [symbol, symbol+2^ignoreBits).
    const numLevels = this.numLevels;
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBitMask = 1 << (this.maxLevel - l);
      let k = sortedOffsets.length; // march k over the sorted offsets array
      for (let i = walk.length; i > 0; ) {
        i--;
        const first = F[i];
        const first1 = level.rank1(first - 1);
        const first0 = first - first1;

        const last = L[i];
        const last1 = level.rank1(last - 1);
        const last0 = last - last1;
        const leftChildCount = last0 - first0; // left child count

        const symbol = S[i];
        const sortedOffsetCount = C[i]; // number of sorted offsets inside this node
        nRankCalls += 2;

        // Determine the number of nodes that wants to be mapped to the right child of this node,
        // then subtract the count of left children from all of the nodes matched to the right child,
        // to account for the elements counted in the left counted.
        // [lo, hi) is the range of sorted symbols covered by this node
        const lo = k - sortedOffsetCount;
        const hi = k;

        const splitIndex = binarySearchBefore(I, leftChildCount, lo, hi);
        const numGoLeft = splitIndex - lo;
        const numGoRight = sortedOffsetCount - numGoLeft;
        k = lo;
        // adjust count for quantiles mapped to the right child based on the left count,
        // so that we look only for the remaining count of elements in the child node.
        for (let n = splitIndex; n < hi; n++) I[n] -= leftChildCount;

        if (numGoRight > 0) {
          // go right
          const nextIndex = walk.nextBackIndex();
          F[nextIndex] = nz + (first - first0); // = nz + first1
          L[nextIndex] = nz + (last - last0); // = nz + last1
          S[nextIndex] = symbol | levelBitMask;
          C[nextIndex] = numGoRight;
        }

        if (numGoLeft > 0) {
          // go left
          const nextIndex = walk.nextBackIndex();
          F[nextIndex] = first0;
          L[nextIndex] = last0;
          S[nextIndex] = symbol;
          C[nextIndex] = numGoLeft;
        }
      }

      walk.reset(false, F, L, S, C);
    }

    for (let i = 0; i < walk.length; i++) L[i] -= F[i];
    const counts = L.subarray(0, walk.length).slice();
    const symbols = S.subarray(0, walk.length).slice();
    // numSortedOffsets indicates how many entries in sortedOffsets are assigned
    // to each symbol. this is a more economical representation than a dense array of
    // length sortedOffsets when multiple sortedOffsets point to the same symbol.
    const numSortedOffsets = C.subarray(0, walk.length).slice();
    return { symbols, counts, numSortedOffsets, nRankCalls };
  }

  // The approach below is from by "New algorithms on wavelet trees and applications to information retrieval"
  quantiles({ first = 0, last = this.length, firstOffset = 0, lastOffset = last - first } = {}) {
    if (last > this.length) throw new Error('last must be <= wavelet matrix length');
    if (first >= last) return this.emptyResult();
    if (firstOffset > lastOffset) throw new Error('firstOffset must be <= lastOffset');
    if (firstOffset < 0 || lastOffset > last - first)
      throw new Error('sortedIndex cannot be less than zero or exceed length of range [first, last)');
    const scratchLength = Math.min(this.alphabetSize, lastOffset - firstOffset);
    this.scratch.reset();
    const F = this.scratch.allocU32(scratchLength); // firsts
    const L = this.scratch.allocU32(scratchLength); // lasts
    const S = this.scratch.allocU32(scratchLength); // symbols
    const C = this.scratch.allocU32(scratchLength); // counts up to the first offset
    const C2 = this.scratch.allocU32(scratchLength); // counts up to the last offset
    const walk = new ArrayWalker(firstOffset === lastOffset ? 0 : 1, scratchLength);
    const nextIndex = walk.nextFrontIndex();
    F[nextIndex] = first;
    L[nextIndex] = last;
    S[nextIndex] = 0;
    C[nextIndex] = firstOffset;
    C2[nextIndex] = lastOffset;
    walk.reset(true, F, L, S, C, C2);
    let nRankCalls = 0;

    const numLevels = this.numLevels;
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBitMask = 1 << (this.maxLevel - l);
      for (let i = walk.length; i > 0; ) {
        i--;
        const first = F[i];
        const first1 = level.rank1(first - 1);
        const first0 = first - first1;

        const last = L[i];
        const last1 = level.rank1(last - 1);
        const last0 = last - last1;

        const symbol = S[i];
        const c = C[i];
        const c2 = C2[i];

        nRankCalls += 2;

        const leftChildCount = last0 - first0;

        if (c2 > leftChildCount) {
          // go right
          const nextIndex = walk.nextBackIndex();
          F[nextIndex] = nz + (first - first0); // = nz + first1
          L[nextIndex] = nz + (last - last0); // = nz + last1
          S[nextIndex] = symbol | levelBitMask;
          C[nextIndex] = leftChildCount < c ? c - leftChildCount : 0; // = max(c - leftChildCount, 0);
          C2[nextIndex] = c2 - leftChildCount;
        }

        if (c < leftChildCount) {
          // go left
          const nextIndex = walk.nextBackIndex();
          F[nextIndex] = first0;
          L[nextIndex] = last0;
          S[nextIndex] = symbol;
          C[nextIndex] = c;
          C2[nextIndex] = leftChildCount < c2 ? leftChildCount : c2; // = min(leftChildCount, c2);
        }
      }
      walk.reset(false, F, L, S, C, C2);
    }

    for (let i = 0; i < walk.length; i++) L[i] -= F[i];
    const counts = L.subarray(0, walk.length).slice();
    const symbols = S.subarray(0, walk.length).slice();
    return { symbols, counts, nRankCalls };
  }

  simpleMajority({ first = 0, last = this.length } = {}) {
    const index = first + Math.trunc((last - first) / 2);
    const q = this.quantile(index, { first, last });
    const total = last - first;
    const half = Math.trunc(total / 2);
    if (q.count > half) return q;
    return null;
  }

  majority(k, { first = 0, last = this.length } = {}) {
    // Returns the 1/k-majority. Ie. for k = 4, return the elements (if any) with
    // frequency larger than 1/4th (25%) of the specified index range.
    // `k` is the number of evenly-spaced samples that will be made in the range [first, last).
    if (k < 1 || !Number.isInteger(k)) throw new Error('k must be a positive integer');
    const total = last - first; // todo: change if inclusive
    // note: can oversample, eg. if k > total.
    // note that total includes multiplicity, so the numbers can be quite large...
    // if k === 1, we sample the first element at index 0
    const indices = new Uint32Array(Math.max(1, k - 1));
    for (let i = 1; i < k; i++) {
      const pc = i / k;
      // implicit floor; we consistently round down.
      // quantileBatch indices are in the range [0, last - first).
      indices[i - 1] = total * pc;
    }
    // Filter out results that do not have a count high enough to be a k-majority.
    const res = this.quantileBatch(indices, { first, last });
    const targetCount = Math.floor((last - first) / k);
    let n = 0;
    for (let i = 0; i < res.symbols.length; i++) {
      if (res.counts[i] > targetCount) {
        res.counts[n] = res.counts[i];
        res.symbols[n] = res.symbols[i];
        n += 1;
      }
    }
    return { symbols: res.symbols.subarray(0, n), counts: res.counts.subarray(0, n) };
  }

  // set-union the symbols in two index ranges
  union(opts) {
    return this.setOp({
      ...opts,
      // recuse if the desired symbols appear in either index range
      op: (a, b) => a > 0 || b > 0,
    });
  }

  // set-intersect the symbols in two index ranges
  intersect(opts) {
    return this.setOp({
      ...opts,
      // recuse if the desired symbols appear in both index ranges
      op: (a, b) => a > 0 && b > 0,
    });
  }

  // Range intersection. Supports the same arguments as `count`, plus first2 and last2. Draft implementation.
  // The idea is to traverse the tree as usual, but only recurse if both of the intervals are non-empty.
  // it would be useful to generalize this to intersect more ranges, accepting `firsts` and `lasts`.
  // i think we would need even more scratch space for the intermediate rank computations, though...
  // We could also support finding elements that appear only in the first range, by recursing only when
  // leftCount > 0 && leftCount2 == 0 for left, and rightCount > 0 && rightCount2 == 0 for right
  setOp({
    // todo: ignorebits
    // todo: for union, we currently keep recursing down nodes that are zero count,
    //       and doing rank calls for them. Is there a way to avoid that without separately implementing
    //       each set op (only compute rank if first!=last or something)
    //       [update: implemented a prototype via the empty check.]
    // operation taking (count, count2) node counts and returning true or false – whether to recurse
    // into subnodes with those counts.
    op,
    first,
    last,
    first2,
    last2,
    lower = 0,
    upper = this.maxSymbol,
    subcodeSeparator = 0,
    sort = false,
  } = {}) {
    if (first === undefined || last === undefined || first2 === undefined || last2 === undefined) {
      throw new Error('first, last, first2, and last2 must all be specified');
    }
    if (lower > upper) throw new Error('lower must be < upper');
    if (last > this.length) throw new Error('last must be <= wavelet matrix length');
    if (first >= last) return Object.assign(this.emptyResult(), { counts2: new Uint32Array() });
    if (last2 > this.length) throw new Error('last2 must be < wavelet matrix length');
    if (first2 > last2) throw new Error('first2 must be <= last2');

    // for intersect
    // const scratchLength = Math.min(upper - lower + 1, last - first, last2 - first2);
    // could pass this in to be smaller for different kinds of set ops, eg. when constrained
    // by the smallest set like in union
    const scratchLength = Math.min(upper - lower + 1);
    this.scratch.reset();
    const F = this.scratch.allocU32(scratchLength); // firsts
    const L = this.scratch.allocU32(scratchLength); // lasts
    const F2 = this.scratch.allocU32(scratchLength); // firsts2
    const L2 = this.scratch.allocU32(scratchLength); // lasts2
    const S = this.scratch.allocU32(scratchLength); // symbols
    const walk = new ArrayWalker(1, scratchLength);
    const reverse = !sort; // for walk.reset(reverse, ...)
    const nextIndex = walk.nextFrontIndex();
    F[nextIndex] = first;
    L[nextIndex] = last;
    F2[nextIndex] = first2;
    L2[nextIndex] = last2;
    S[nextIndex] = 0;
    walk.reset(reverse, F, L, F2, L2, S);
    let nRankCalls = 0;

    const numLevels = this.numLevels;
    // see `counts` for subcode-related comments.
    let subcodeMask = 0xffffffff << numLevels;

    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBit = this.maxLevel - l;
      const levelBitMask = 1 << levelBit;

      if ((subcodeSeparator & levelBitMask) === 0) subcodeMask |= levelBitMask;
      else subcodeMask = levelBitMask;
      const subcodeLower = lower & subcodeMask;
      const subcodeUpper = upper & subcodeMask;

      // if we want sorted outputs, iterate in reverse to ensure that we don't
      // overwrite unprocessed elements when writing from the back of the array.
      const start = sort ? walk.length - 1 : 0;
      const step = sort ? -1 : 1;
      const end = sort ? -1 : walk.length;
      for (let i = start; i != end; i += step) {
        const first = F[i];
        const last = L[i];
        const empty = first === last; // node count == 0

        const first1 = empty ? first : level.rank1(first - 1);
        const first0 = first - first1;

        const last1 = empty ? last : level.rank1(last - 1);
        const last0 = last - last1;

        const rightCount = last1 - first1;
        const leftCount = last0 - first0;

        const first2 = F2[i];
        const last2 = L2[i];
        const empty2 = first2 === last2; // node count == 0

        const first1_2 = empty ? first1 : level.rank1(first2 - 1);
        const first0_2 = first2 - first1_2;

        const last1_2 = empty ? last1 : level.rank1(last2 - 1);
        const last0_2 = last2 - last1_2;

        const rightCount2 = last1_2 - first1_2;
        const leftCount2 = last0_2 - first0_2;

        const symbol = S[i];

        nRankCalls += 4;

        if (op(rightCount, rightCount2)) {
          // go right [if symbol ranges overlap each other and the target range]
          const a = (symbol | levelBitMask) & subcodeMask;
          const b = (a | (levelBitMask - 1)) & subcodeMask;
          if (intervalsOverlapInclusive(a, b, subcodeLower, subcodeUpper)) {
            const nextIndex = walk.nextBackIndex();
            F[nextIndex] = nz + first1;
            L[nextIndex] = nz + last1;
            F2[nextIndex] = nz + first1_2;
            L2[nextIndex] = nz + last1_2;
            S[nextIndex] = symbol | levelBitMask;
          }
        }

        if (op(leftCount, leftCount2)) {
          // go left
          const a = symbol & subcodeMask;
          const b = (a | (levelBitMask - 1)) & subcodeMask;
          if (intervalsOverlapInclusive(a, b, subcodeLower, subcodeUpper)) {
            const nextIndex = sort ? walk.nextBackIndex() : walk.nextFrontIndex();
            F[nextIndex] = first0;
            L[nextIndex] = last0;
            F2[nextIndex] = first0_2;
            L2[nextIndex] = last0_2;
            S[nextIndex] = symbol;
          }
        }
      }
      walk.reset(reverse, F, L, F2, L2, S);
    }

    for (let i = 0; i < walk.length; i++) {
      L[i] -= F[i];
      L2[i] -= F2[i];
    }
    const counts = L.subarray(0, walk.length).slice();
    const counts2 = L2.subarray(0, walk.length).slice();
    const symbols = S.subarray(0, walk.length).slice();
    return { symbols, counts, counts2, nRankCalls };
  }

  // status: initial draft.
  // note: will usually return more than required, since that is what was
  // retrieved.
  // note: does not return them in sorted order. we have a 'topK' function for that instead.
  // figure out how to make this less confusing.
  mostFrequentUnsorted(k, { first = 0, last = this.length } = {}) {
    // todo: have some heuristic about when it might make sense to jump straight to querying counts,
    // or having a different strategy than doubing (when k is small, eg. 10).
    let ret = this.majority(k, { first, last });
    let iters = 0;
    while (ret.symbols.length < k) {
      iters++;
      const n = k * 2 ** iters;
      // console.log(iters, n, ret, ret.symbols.length);
      // note:
      // we may want a different strategy than when n > alphabetSize for breaking out.
      // also, the iters breakpoint is a magic number. should it be configurable?
      if (n > this.alphabetSize || iters > 6) break;
      ret = this.majority(n, { first, last });
    }
    if (ret.symbols.length < k) {
      ret = this.counts({ first, last, sort: false });
    }
    return ret;
  }

  // this is all hard to use... maybe encodeSubcodes should take a list of subcode sizes
  // as integers

  subcodeSeparator(subcodeSizesInBits) {
    let separator = 0;
    let offset = 0;
    for (const sz of subcodeSizesInBits) {
      if (sz === 0) throw 'cannot have zero-sized field';
      separator |= 1 << (sz - 1 + offset);
      offset += sz;
    }
    return separator >>> 0;
  }

  encodeSubcodes(subcodeSeparator, values) {
    if (subcodeSeparator === 0) {
      if (values.length !== 1) {
        throw new Error('number of values must be one if the subcodeSeparator is zero');
      }
      return values[0];
    }
    if (popcount(subcodeSeparator) !== values.length) {
      throw new Error('number of values must be equal to the number of 1 bits in the subcodeSeparator');
    }
    let code = 0;
    let offset = 0;
    let i = 0;
    while (subcodeSeparator > 0) {
      // todo: validate that values[i] is 0 <= v < 2^subcodeSize
      const subcodeSize = trailing0(subcodeSeparator) + 1;
      code |= values[i] << offset;
      i += 1;
      offset += subcodeSize;
      subcodeSeparator >>>= subcodeSize; // shift off this subcode
    }
    return code >>> 0;
  }

  approxSizeInBits() {
    // ignores fixed-size fields
    let size = 0;
    for (const level of this.levels) {
      size += level.approxSizeInBits();
    }
    return size;
  }
}

// note:
// we can express rank1 in terms of rank0
// when i and j are both in bounds:
//    rank0(i)   = i - rank1(i) + 1
// => rank0(i-1) = i-1 - rank1(i-1) + 1
// => rank0(i-1) = i - rank1(i-1)
// and vice versa, so i0 = i - i1; (see impl. of rank0)

// Test two intervals for inclusive overlap.
function intervalsOverlapInclusive(aLo, aHi, bLo, bHi) {
  return aLo <= bHi && bLo <= aHi;
}
