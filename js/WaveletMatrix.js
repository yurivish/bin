import { BitVector } from './BitVector.js';
import { ZeroCompressedBitVector } from './ZeroCompressedBitVector.js';
import { reverseBits, reverseBits32, clamp, trailing0, popcount, isObjectLiteral } from './util.js';

// Implements a binary wavelet matrix that splits on power-of-two alphabet
// boundaries, rather than splitting based on the true alphabet midpoint.
export class WaveletMatrix {
  // This implements Algorithm 1 (seq.pc) from the paper "Practical Wavelet Tree Construction".
  // The return value is the array of wavelet tree levels. Adapting the algorithm to construct
  // a wavelet matrix instead requires changing the borders computation (see section 5.3).
  // todo: check that all symbols are < alphabetSize
  // todo: pass in maxSymbol, with alphabetSize = maxSymbol + 1?
  constructor(data, alphabetSize, opts = {}) {
    // As a simple heuristic, by default use the large alphabet constructor when
    // the alphabet sides exceeds the number of data points.
    const { largeAlphabet = alphabetSize > data.length } = opts;
    if (largeAlphabet) return this.constructLargeAlphabet(data, alphabetSize, opts);
    // data is an array of integer values in [0, alphabetSize)
    const numLevels = Math.ceil(Math.log2(alphabetSize));
    const maxLevel = numLevels - 1;
    // todo: can we get away with non-pow2, storing just one entry per symbol?
    const hist = new Uint32Array(2 ** numLevels);
    const borders = new Uint32Array(2 ** numLevels);
    // note: if data is sorted, could we use a compressed set data structure for hist/borders?
    // todo: can we perform better for sparse code sets, eg. [0, 2^32)?

    const levels = new Array(numLevels);
    // Initialize the level bit vectors
    for (let i = 0; i < numLevels; i++) {
      levels[i] = new BitVector(data.length);
      // try this once bits can be added to it out-of-order
      // levels[i] = new ZeroCompressedBitVector(data.length, { rank: true });
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
        // Update the positions in-place (wavelet matrix)
        const prevIndex = reverseBits(i - 1, l);
        borders[reverseBits(i, l)] = borders[prevIndex] + hist[prevIndex];
        // To compute a wavelet tree, do this instead:
        // borders[i] = borders[i - 1] + hist[i - 1];
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
        const nodeIndex = (d & bitPrefixMask) >>> bitPrefixShift;
        const p = borders[nodeIndex];
        borders[nodeIndex] += 1;
        // Set the bit in the bitvector
        if (d & levelBitMask) level.one(p);
      }
    }

    // Compute the number of zeros at each level
    const numZeros = new Uint32Array(numLevels);
    for (let i = 0; i < numLevels; i++) {
      levels[i].finish();
      numZeros[i] = levels[i].rank0(levels[i].length);
    }

    this.levels = levels;
    this.alphabetSize = alphabetSize;
    this.numZeros = numZeros;
    this.numLevels = numLevels;
    this.maxLevel = maxLevel;
    this.length = data.length;
    this.scratch = new ScratchSpace();
  }

  // Alternative construction algorithm for the 'sparse' case when the alphabet size
  // is significantly larger than the number of symbols that actually occur in the data.
  constructLargeAlphabet(data, alphabetSize, opts = {}) {
    data = new Uint32Array(data); // copy data because will be mutated
    let next = new Uint32Array(data.length);

    // data is an array of integer values in [0, alphabetSize)
    const numLevels = Math.ceil(Math.log2(alphabetSize));
    const maxLevel = numLevels - 1;

    // Initialize the level bit vectors
    const levels = new Array(numLevels);
    for (let i = 0; i < numLevels; i++) {
      levels[i] = new BitVector(data.length);
    }

    const numZeros = new Uint32Array(numLevels);
    const walk = new ArrayWalker(0, data.length);

    // For each level, sort the data point by its bit value at that level.
    // Zero bits get sorted left, one bits get sorted right. This amounts
    // to a bucket sort with two buckets.
    // We sort into `next`, then swap `next` and `data`.
    for (let l = 0; l < maxLevel; l++) {
      const level = levels[l];
      const levelBit = maxLevel - l;
      const levelBitMask = 1 << levelBit;
      for (let i = 0; i < data.length; i++) {
        const d = data[i];
        if (d & levelBitMask) {
          next[walk.nextBackIndex()] = d;
          level.one(i);
        } else {
          next[walk.nextFrontIndex()] = d;
        }
        numZeros[l] = walk.frontIndex;
      }
      walk.reset(true, next);
      const tmp = data;
      data = next;
      next = tmp;
    }

    // For the last level we don't need to build anything but the bitvector
    const level = levels[maxLevel];
    const levelBitMask = 1 << 0;
    for (let i = 0; i < data.length; i++) {
      if (data[i] & levelBitMask) level.one(i);
    }
    numZeros[maxLevel] = level.rank0(level.length);

    // Mark the level bitvectors as finished
    for (let l = 0; l < numLevels; l++) levels[l].finish();

    this.levels = levels;
    this.alphabetSize = alphabetSize;
    this.numZeros = numZeros;
    this.numLevels = numLevels;
    this.maxLevel = maxLevel;
    this.length = data.length;
    this.scratch = new ScratchSpace();
  }

  access(index) {
    if (index < 0 || index > this.length) throw new Error('access: out of bounds');
    let symbol = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const index1 = level.rank1(index - 1); // todo: combine rank and access queries since they symbol the same block
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

  // Returns the number of occurrences of `symbol` in the range [first, last).
  countSymbol(first, last, symbol, opts) {
    const indices = this.symbolIndices(first, last, symbol, opts);
    return indices.last - indices.first;
  }

  // Returns the index of the kth occurrence of `symbol` in the range [first, last).
  select(first, last, symbol, k) {
    if (symbol < 0 || symbol >= this.alphabetSize) return -1;
    if (k < 1 || symbol > this.length) return -1;
    const indices = this.symbolIndices(first, last, symbol, 0);
    if (indices.last - indices.first < k) return -1; // in analogy with select
    let index = indices.first + k - 1;
    for (let l = this.numLevels; l > 0; ) {
      l -= 1;
      const level = this.levels[l];
      const nz = this.numZeros[l];
      if (index < nz) {
        // this position was mapped from a zero at the previous level
        const k = index + 1; // this was the nth zero on the that level
        index = level.select0(k); // locate the corresponding index on the preceding level
      } else {
        // this position was mapped from a one at the previous level
        const k = index - nz + 1; // this was the nth one on the that level
        index = level.select1(k); // locate the corresponding index on the preceding level
      }
    }
    return index;
  }

  // Internal function returning the (first, last] index range covered by this symbol on the virtual bottom level,
  // or on a higher level if groupBits > 0. `last - first` gives the symbol count within the provided range.
  symbolIndices(first, last, symbol, { groupBits = 0 } = {}) {
    const symbolGroupSize = 1 << groupBits;
    if (symbol % symbolGroupSize !== 0) {
      // note: could be done with bit math (check that low bits are zero)
      throw new Error('symbol must evenly divide the block size implied by groupBits');
    }
    if (symbol >= this.alphabetSize) throw new Error('symbol must be < alphabetSize');
    if (first > last) throw new Error('last must be <= first');
    if (first === last) return 0;
    if (first > this.length) throw new Error('first must be < wavelet matrix length');
    if (this.numLevels === 0) return 0;
    const numLevels = this.numLevels - groupBits;
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
  countLessThan(first, last, symbol) {
    if (first < 0) throw new Error('first must be >= 0');
    if (first > last) throw new Error('first must be <= last');
    if (last > this.length) throw new Error('last must be < wavelet matrix length');
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
  // of rank calls, but this comes at a cost of implementation complexity.
  // See the paper "New algorithms on wavelet trees and applications to
  // information retrieval" for details. Another approach is to modify the
  // implementation of countLessThan to perform two interleaved calls. The
  // subtlety there is that, as written, the algorithm does not work when
  // symbol >= alphabetSize (the one-symbol impl. can return early in this case).
  count(first, last, lower, upper) {
    return this.countLessThan(first, last, upper) - this.countLessThan(first, last, lower);
  }

  countSymbolBatch(first, last, sortedSymbols, { groupBits = 0 } = {}) {
    // splitByMsb requires the same sortedSymbol to be searched for in each of the split paths.
    // for now, we'll go with the relatively inefficient route of asking that this be done by
    // supplying a larger set of sortedSymbols, enumerating all MSB variations in the high bits.
    const symbolGroupSize = 1 << groupBits;
    for (const symbol of sortedSymbols) {
      if (symbol % symbolGroupSize !== 0)
        // note: could be done with bit math (check that low bits are zero)
        throw new Error('symbol must evenly divide the block size implied by groupBits');
    }
    if (first > last) throw new Error('first must be <= last');
    if (last > this.length) throw new Error('last must be < wavelet matrix length');

    // account for duplicate sorted symbols and the fact that there cannot be more outputs than elements
    const scratchSize = Math.min(this.alphabetSize, sortedSymbols.length, last - first);
    this.scratch.reset();
    const F = this.scratch.alloc(scratchSize); // firsts
    const L = this.scratch.alloc(scratchSize); // lasts
    const S = this.scratch.alloc(scratchSize); // symbols
    const C = this.scratch.alloc(scratchSize); // counts
    const walk = new ArrayWalker(sortedSymbols.length === 0 ? 0 : 1, scratchSize);
    const nextIndex = walk.nextFrontIndex();
    F[nextIndex] = first;
    L[nextIndex] = last;
    C[nextIndex] = sortedSymbols.length;
    S[nextIndex] = 0;
    walk.reset(false, F, L, C, S);
    let nRankCalls = 0;

    const numLevels = this.numLevels - groupBits;
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

        // perform two binary searches over the sorted symbols for this node to deterine
        // the number in the left child node. We do a single linear pass across sortedSymbols
        // over the course of the walk.
        // Note that this relies on the fact that the symbols in `S` are in sorted order,
        // since we're marching `k` across the sorted symbols array to reduce the search
        // space for the binary search calls to those symbols corresponding to this node.
        // This means that we have to use `walk.nextBackIndex()` for both left and right children,
        // since that's what gives rise to the sorted order of `S` its sorted order.
        // [lo, hi) is the range of sorted indices covered by this node
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
  // in groups of size 2^groupBits (symbols are grouped together when they
  // differ only in their lowest `groupBits` bits)
  // Each distinct group is labeled by its lowest element, which represents
  // the group containing symbols in the range [symbol, symbol + 2^groupBits).
  counts(first, last, lower, upper, { groupBits = 0, subcodeIndicator = 0, sort = true } = {}) {
    const symbolGroupSize = 1 << groupBits;
    // todo: handle lower === upper
    // these error messages could be improved, explaining that ignore bits tells us the power of two
    // that lower and upper need to be multiples of.
    if (lower % symbolGroupSize !== 0)
      throw new Error('lower must evenly divide the symbol block size implied by groupBits');
    if ((upper + 1) % symbolGroupSize !== 0)
      throw new Error('(upper + 1) must evenly divide the symbol block size implied by groupBits');
    // if (upper >= this.alphabetSize) throw new Error('upper must be < alphabetSize ([lower, upper] is inclusive)');
    // ^ we now allow this so that subcode stuff works without us having to be unrealistic about the true alphabet size
    // (eg. allow querying code consisting of all maximum subcodes)
    const numLevels = this.numLevels - groupBits;

    // todo: bound this more closely. slightly involved to upper-bound due to subcodes;
    // need to compute the product of the subcode ranges since that's the maximum possible
    // number of unique symbols.
    const scratchSize = Math.min(2 ** this.numLevels - groupBits, this.alphabetSize, this.length);
    this.scratch.reset();
    const F = this.scratch.alloc(scratchSize);
    const L = this.scratch.alloc(scratchSize);
    const S = this.scratch.alloc(scratchSize);
    const walk = new ArrayWalker(1, scratchSize);
    const reverse = !sort; // walk.reset(reverse, ...)
    const nextIndex = walk.nextFrontIndex();
    F[nextIndex] = first;
    L[nextIndex] = last;
    S[nextIndex] = 0;
    walk.reset(false, F, L, S);

    // count inner
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
      // It can be useful to instead treat the code as representing a concatenation of subcodeIndicator,
      // and the [lower, upper] values as representing a concatenation of the ranges of those
      // subcodeIndicator. This behavior can be specified by the `subcodeIndicator` argument, which is a
      // bitmask in which a 1 bit indicates the onset of a new subcode and a 0 implies the continuation
      // of the current subcode. All range comparisons are done within a subcode, and the default
      // subcodeIndicator of 0 gives us the default behavior in which the full code is treated as
      // a single subcode.
      if ((subcodeIndicator & levelBitMask) === 0) subcodeMask |= levelBitMask;
      else subcodeMask = levelBitMask;
      const subcodeLower = lower & subcodeMask;
      const subcodeUpper = upper & subcodeMask;

      // note: if we want sorted outputs, iterate in reverse to ensure that we don't
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
          // go right if the right node range [a, b] (inclusive) overlaps [lower, upper) (exclusive)
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
          // go left if the left node range [a, b] (inclusive) overlaps [lower, upper] (inclusive)
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

  // Adapted from https://github.com/noshi91/Library/blob/0db552066eaf8655e0f3a4ae523dbf8c9af5299a/data_structure/wavelet_matrix.cpp#L76
  // Range quantile query returning the kth largest symbol in A[i, j).
  quantile(first, last, index) {
    if (first > last) throw new Error('first must be <= last');
    if (last > this.length) throw new Error('last must be < wavelet matrix length');
    if (index < 0 || index >= last - first)
      throw new Error('index cannot be less than zero or exceed length of range [first, last)');
    let symbol = 0;
    let nRankCalls = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const first0 = level.rank0(first - 1);
      const last0 = level.rank0(last - 1);
      const count = last0 - first0;
      nRankCalls += 2;
      if (index < count) {
        // go left
        first = first0;
        last = last0;
      } else {
        // go right
        const nz = this.numZeros[l];
        first = nz + (first - first0); // = nz + first1
        last = nz + (last - last0); // = nz + last1
        // update symbol and new target sorted index in the child node
        const levelBitMask = 1 << (this.maxLevel - l);
        symbol |= levelBitMask;
        index -= count;
      }
    }
    return { symbol, count: last - first, nRankCalls };
  }

  quantileBatch(first, last, sortedIndices) {
    // these error messages could be improved, explaining that ignore bits tells us the power of two
    // that lower and upper need to be multiples of.
    if (first > last) throw new Error('first must be <= last');
    if (last > this.length) throw new Error('last must be < wavelet matrix length');
    for (let i = 1; i < sortedIndices.length; i++) {
      if (!(sortedIndices[i - 1] <= sortedIndices[i])) throw new Error('sorted indices must be sorted');
    }
    // todo: error if there are more sortedindices than last-first, since we copy them into a scratch space (or ensure the space can hold the size we need)
    if (sortedIndices[0] < 0 || sortedIndices[sortedIndices.length - 1] >= last - first)
      throw new Error('sortedIndex cannot be less than zero or exceed length of range [first, last)');

    // account for duplicate sorted indices and the fact that there cannot be more outputs than elements
    const scratchSize = Math.min(this.alphabetSize, sortedIndices.length, last - first);
    this.scratch.reset();
    const F = this.scratch.alloc(scratchSize); // firsts
    const L = this.scratch.alloc(scratchSize); // lasts
    const S = this.scratch.alloc(scratchSize); // symbols
    const C = this.scratch.alloc(scratchSize); // counts
    const I = this.scratch.alloc(sortedIndices.length); // sorted indices
    const walk = new ArrayWalker(sortedIndices.length === 0 ? 0 : 1, scratchSize);
    const nextIndex = walk.nextFrontIndex();
    F[nextIndex] = first;
    L[nextIndex] = last;
    S[nextIndex] = 0;
    C[nextIndex] = sortedIndices.length; // number of sortedIndices represented by node 0
    walk.reset(true, F, L, S, C);
    // copy sorted indices into a scratch space since they are mutated as we go
    I.set(sortedIndices);
    let nRankCalls = 0;

    // note: grouping by LSB computes approximate quantiles where
    // the count of symbols assigned to each range is given by I,
    // and I think the ranges are [symbol, symbol+2^groupBits).
    const numLevels = this.numLevels;
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBitMask = 1 << (this.maxLevel - l);
      let k = sortedIndices.length; // march k over the sorted indices array
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
        const sortedIndexCount = C[i]; // number of sorted indices inside this node
        nRankCalls += 2;

        // Determine the number of nodes that wants to be mapped to the right child of this node,
        // then subtract the count of left children from all of the nodes matched to the right child,
        // to account for the elements counted in the left counted.
        // [lo, hi) is the range of sorted symbols covered by this node
        const lo = k - sortedIndexCount;
        const hi = k;

        const splitIndex = binarySearchBefore(I, leftChildCount, lo, hi);
        const numGoLeft = splitIndex - lo;
        const numGoRight = sortedIndexCount - numGoLeft;
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
    // numSortedIndices indicates how many entries in sortedIndices are assigned
    // to each symbol. this is a more economical representation than a dense array of
    // length sortedIndices when multiple sortedIndices point to the same symbol.
    const numSortedIndices = C.subarray(0, walk.length).slice();
    return { symbols, counts, numSortedIndices, nRankCalls };
  }

  // The approach below is from by "New algorithms on wavelet trees and applications to information retrieval"
  quantiles(first, last, firstIndex, lastIndex) {
    if (first > last) throw new Error('first must be <= last');
    if (last > this.length) throw new Error('last must be < wavelet matrix length');
    if (firstIndex > lastIndex) throw new Error('firstIndex must be <= lastIndex');
    if (firstIndex < 0 || lastIndex > last - first)
      throw new Error('sortedIndex cannot be less than zero or exceed length of range [first, last)');
    const scratchSize = lastIndex - firstIndex;
    this.scratch.reset();
    const F = this.scratch.alloc(scratchSize); // firsts
    const L = this.scratch.alloc(scratchSize); // lasts
    const S = this.scratch.alloc(scratchSize); // symbols
    const C = this.scratch.alloc(scratchSize); // counts
    const C2 = this.scratch.alloc(scratchSize);
    const walk = new ArrayWalker(firstIndex === lastIndex ? 0 : 1, scratchSize);
    const nextIndex = walk.nextFrontIndex();
    F[nextIndex] = first;
    L[nextIndex] = last;
    S[nextIndex] = 0;
    C[nextIndex] = firstIndex;
    C2[nextIndex] = lastIndex;
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

  simpleMajority(first, last) {
    const index = (last + first) >>> 1;
    const q = this.quantile(first, last, index);
    const total = last - first;
    const half = total >>> 1;
    if (q.count > half) return q;
    return null;
  }

  majority(first, last, denominator = 2) {
    if (denominator < 1) throw new Error('denominator must be a positive integer');
    // if denominator === 1, we sample the first element at index 0
    const indices = new Uint32Array(Math.max(1, denominator - 1));
    const total = last - first; // todo: change if inclusive
    for (let i = 1; i < denominator; i++) {
      const pc = i / denominator;
      // implicit floor; we consistently round down.
      indices[i - 1] = first + total * pc;
    }
    const res = this.quantileBatch(first, last, indices);
    const count = Math.floor((last - first) / denominator);
    let n = 0;
    for (let i = 0; i < res.symbols.length; i++) {
      if (res.counts[i] > count) {
        res.counts[n] = res.counts[i];
        res.symbols[n] = res.symbols[i];
        n += 1;
      }
    }
    return { symbols: res.symbols.subarray(0, n), counts: res.counts.subarray(0, n) };
  }

  subcodeIndicator(subcodeSizesInBits) {
    let indicator = 0;
    let offset = 0;
    for (const sz of subcodeSizesInBits) {
      if (sz === 0) throw 'cannot have zero-sized field';
      indicator |= 1 << (sz - 1 + offset);
      offset += sz;
    }
    return indicator >>> 0;
  }

  encodeSubcodes(indicator, values) {
    if (indicator === 0) {
      if (values.length !== 1) {
        throw new Error('number of values must be one if the indicator is zero');
      }
      return values[0];
    }
    if (popcount(indicator) !== values.length) {
      throw new Error('number of values must be equal to the number of 1 bits in the indicator');
    }
    let code = 0;
    let offset = 0;
    let i = 0;
    while (indicator > 0) {
      // todo: validate that values[i] is 0 <= v < 2^subcodeSize
      const subcodeSize = trailing0(indicator) + 1;
      code |= values[i] << offset;
      i += 1;
      offset += subcodeSize;
      indicator >>>= subcodeSize; // shift off this subcode
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

// Find the rightmost insertion index in A for T
// in order to maintain A's sorted order.
// Searches the index range [L, R).
function binarySearchAfter(A, T, L, R) {
  while (L < R) {
    // This midpoint calculation will return incorrect results for large arrays (>2^30)
    // By that point we should switch to Zig. Correct alternative: Math.floor(L + (R - L) / 2);
    const m = (L + R) >>> 1;
    if (A[m] > T) R = m;
    else L = m + 1;
  }
  return R;
}

// Find the leftmost insertion index in A for T
// in order to maintain A's sorted order.
// Searches the index range [L, R).
function binarySearchBefore(A, T, L, R) {
  while (L < R) {
    // This midpoint calculation will return incorrect results for large arrays (>2^30)
    // By that point we should switch to Zig. Correct alternative: Math.floor(L + (R - L) / 2);
    const m = (L + R) >>> 1;
    if (A[m] < T) L = m + 1;
    else R = m;
  }
  return L;
}

// note:
// we can express rank1 in terms of rank0
// when i and j are both in bounds:
//    rank0(i)   = i - rank1(i) + 1
// => rank0(i-1) = i-1 - rank1(i-1) + 1
// => rank0(i-1) = i - rank1(i-1)
// and vice versa, so i0 = i - i1; (see impl. of rank0)

// helper for turning recursion into iteration, encapsulating the logic
// of walking an array zero to two outputs for each input.
// left elements are filled in from the front, while right elements are filled
// in from the back, and reset() will preserve all left children in order at
// the front of the array, and move the right children in order immediately
// after them. to interleave left and right in the order they are generated,
// use nextBackIndex() for both the left and right outputs (we cannot use nextFrontIndex()
// for both, since that would overwrite elements as they are being processed
// in left-to-right order).
class ArrayWalker {
  constructor(length, cap) {
    this.length = length; // length taken up by existing elements
    this.cap = cap; // capacity for additional elements
    this.frontIndex = 0;
    this.backIndex = cap; // nextIndex + 1
  }

  // Return the next index from the front (left side) of the array
  // note: ensure that the current value has been retrieved before
  // writing to the left, since it will overwrite the current value.
  nextFrontIndex() {
    const index = this.frontIndex;
    this.frontIndex += 1;
    return index;
  }
  // Return the next index from the back (right side) of the array
  nextBackIndex() {
    // return the next index at which we can append an element
    // (as we fill the array in backwards from arr[cap - 1])
    return (this.backIndex -= 1);
  }

  // `reverse`: Whether to reverse the [backIndex...cap] subarray.
  // In practice, we want to reverse it if we've been iterating
  // the source array from left-to-right, and want to not reverse
  // if we've been iterating right-to-left.
  reset(reverse, ...arrays) {
    if (this.backIndex < this.cap) {
      // move the filled-in elements from the end
      // to the front of the array
      for (let i = 0; i < arrays.length; i++) {
        const arr = arrays[i];
        const sub = arr.subarray(this.backIndex, this.cap);
        // reverse (if needed)
        if (reverse) sub.reverse();
        // move right elements to follow the left element directly
        arr.set(sub, this.frontIndex);
      }
    }
    // apply the same logical change to the last and length markers
    this.length = this.frontIndex + (this.cap - this.backIndex);
    this.frontIndex = 0;
    this.backIndex = this.cap;
  }
}

// start scratch with 10kb, then double each time we run out
// alloc, reset. reset at start, alloc whenever we need.
// note: one convenient aspect of this design is that the previous buffer subarrays
// can continue to be used if we are resized during an alloc.
class ScratchSpace {
  constructor(initialLength = 10 * 1024) {
    this.buf = new Uint32Array(initialLength);
    this.index = 0;
  }
  alloc(length) {
    if (this.index + length > this.buf.length) {
      this.buf = new Uint32Array(2 * this.buf.length);
      this.index = 0;
    }
    const sub = this.buf.subarray(this.index, this.index + length);
    this.index += length;
    return sub;
  }
  reset() {
    this.index = 0;
  }
}

// Test two intervals for inclusive overlap.
function intervalsOverlapInclusive(aLo, aHi, bLo, bHi) {
  return aLo <= bHi && bLo <= aHi;
}
