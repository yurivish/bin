import { RankBitVector } from './RankBitVector';
import { ZeroCompressedBitVector } from './ZeroCompressedBitVector';
import { reverseBits, reverseBits32, clamp } from './util';
// todo: range next value, range prev value (though these can be done using quantile), quantiles, majority,
// intersection and more from "New algorithms on wavelet trees and applications to information retrieval"

// wavelet tree is from 2003, wavelet matrix is from oct. 2012
// to do: hoist error objects to the top?

// note: implements a binary wavelet matrix that splits on power-of-two alphabet
// boundaries, rather than splitting based on the true alphabet midpoint.
export class WaveletMatrix {
  // This implements Algorithm 1 (seq.pc) from the paper "Practical Wavelet Tree Construction".
  // The return value is the array of wavelet tree levels. Adapting the algorithm to construct
  // a wavelet matrix instead requires changing the borders computation (see section 5.3).
  // todo: check that all symbols are < alphabetSize
  // todo: pass in maxSymbol, with alphabetSize = maxSymbol + 1?
  constructor(data, alphabetSize) {
    // data is an array of integer values in [0, alphabetSize)
    const numLevels = Math.ceil(Math.log2(alphabetSize));
    const maxLevel = numLevels - 1;

    // todo: can we get away with non-pow2, storing just one entry per symbol?
    const hist = new Uint32Array(2 ** numLevels);
    const borders = new Uint32Array(2 ** numLevels);
    const levels = new Array(numLevels);
    // Initialize the level bit vectors
    for (let i = 0; i < numLevels; i++) {
      levels[i] = new RankBitVector(data.length);
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

    // scratch spaces for intermediate processing.
    // todo: don't materialize these until needed - the alphabet might be big
    // and we can accept scratch space as an input parameter so we can reuse
    // the same space across wavelet trees
    const sz = Math.min(2 ** this.numLevels, this.alphabetSize);
    // todo: buffer pool of scratch spaces
    // todo: these names are also getting quite silly (and inaccurate)
    this.F = new Uint32Array(sz); // firsts
    this.L = new Uint32Array(sz); // lasts
    this.S = new Uint32Array(sz); // symbols
    this.C = new Uint32Array(sz);
    this.C2 = new Uint32Array(sz);
  }

  symbol(index) {
    if (index < 0 || index > this.length) throw new Error('symbol: out of bounds');
    let symbol = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const index1 = level.rank1(index - 1); // todo: combine rank and symbol queries since they symbol the same block
      if (level.symbol(index) === 0) {
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

  count(first, last, symbol, symbolBlockBits = 0) {
    const symbolBlockSize = 1 << symbolBlockBits;
    if (symbol % symbolBlockSize !== 0) {
      // note: could be done with bit math (check that low bits are zero)
      throw new Error('symbol must evenly divide the block size implied by symbolBlockBits');
    }
    if (symbol >= this.alphabetSize) throw new Error('symbol must be < alphabetSize');
    if (first > last) throw new Error('last must be <= first');
    if (first === last) return 0;
    if (first > this.length) throw new Error('first must be < wavelet matrix length');
    if (this.numLevels === 0) return 0;
    const numLevels = this.numLevels - symbolBlockBits;
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
    return last - first;
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
  countRange(first, last, lower, upper) {
    return this.countLessThan(first, last, upper) - this.countLessThan(first, last, lower);
  }

  // todo: think about the behavior wrt. repeated input syms; only returns it once, but unclear what symbol it refers to.
  // I think we should just note that the return value is one symbol per *unique* symbol in the input
  countBatch(first, last, sortedSymbols, symbolBlockBits = 0) {
    const symbolBlockSize = 1 << symbolBlockBits;
    for (const symbol of sortedSymbols) {
      if (symbol % symbolBlockSize !== 0)
        // note: could be done with bit math (check that low bits are zero)
        throw new Error('symbol must evenly divide the block size implied by symbolBlockBits');
    }
    if (first > last) throw new Error('first must be <= last');
    if (last > this.length) throw new Error('last must be < wavelet matrix length');
    const { F, L, C } = this;
    F[0] = first;
    L[0] = last;
    C[0] = sortedSymbols.length;
    const walk = new ReverseArrayWalker(sortedSymbols.length === 0 ? 0 : 1, F.length);
    const numLevels = this.numLevels - symbolBlockBits;
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBitMask = 1 << (this.maxLevel - l); // determines whether a symbol goes to the left or right
      const symbolsPerNode = (levelBitMask << 1) - 1;
      const symbolBitMask = 0xffffffff << (this.maxLevel - l); // clears the low bits from a symbol, giving the left edge of its node
      let k = sortedSymbols.length; // at every level, we sweep through the sorted indices in reverse
      for (let i = walk.len; i > 0; ) {
        i -= 1;

        const first = F[i];
        const first1 = level.rank1(first - 1);
        const first0 = first - first1;

        const last = L[i];
        const last1 = level.rank1(last - 1);
        const last0 = last - last1;

        // [hi, lo) is the index range of the sortedSymbols covered by this node
        const symbolCount = C[i];
        const lo = k - symbolCount;
        const hi = k;
        k -= symbolCount;

        const symbol = sortedSymbols[lo];
        const a = symbol & symbolBitMask; // leftmost symbol in this node
        const b = a + symbolsPerNode; // rightmost symbol in this node
        const m = (a + b) >>> 1;

        // todo: avoid the binary search when all nodes are mapped into the same child
        // by checking sortedSymbols[lo] & levelBitMask and sortedSymbols[hi] & levelBitMask
        const splitIndex = binarySearchAfter(sortedSymbols, m, lo, hi);
        const numGoLeft = splitIndex - lo;
        const numGoRight = symbolCount - numGoLeft;

        if (numGoRight > 0) {
          // go right
          const nextIndex = walk.next();
          F[nextIndex] = nz + first1;
          L[nextIndex] = nz + last1;
          C[nextIndex] = numGoRight;
        }

        if (numGoLeft > 0) {
          // go left
          const nextIndex = walk.next();
          F[nextIndex] = F[i] - first1; // = first0
          L[nextIndex] = L[i] - last1; // = last0
          C[nextIndex] = numGoLeft;
        }
      }

      // update the length and move processed elements back to the front of the list.
      const index = walk.moveToFront();
      F.set(F.subarray(index, walk.cap));
      L.set(L.subarray(index, walk.cap));
      C.set(C.subarray(index, walk.cap));
    }

    for (let i = 0; i < walk.len; i++) L[i] -= F[i];
    const counts = L.subarray(0, walk.len).slice();
    return counts;
  }

  // Returns all of the distinct symbols [lower, upper) in the range [first, last)
  // together with their number of occurrences. Symbols are grouped and processed
  // in groups of size 2^symbolBlockBits (symbols are grouped together when they
  // differ only in their lowest `symbolBlockBits` bits)
  // Each distinct group is labeled by its lowest element, which represents
  // the group containing symbols in the range [symbol, symbol + 2^symbolBlockBits).
  counts(first, last, lower, upper, symbolBlockBits = 0) {
    const symbolBlockSize = 1 << symbolBlockBits;
    // these error messages could be improved, explaining that ignore bits tells us the power of two
    // that lower and upper need to be multiples of.
    if (lower % symbolBlockSize !== 0)
      throw new Error('lower must evenly divide the symbol block size implied by symbolBlockBits');
    if (upper % symbolBlockSize !== 0)
      throw new Error('upper must evenly divide the symbol block size implied by symbolBlockBits');
    const numLevels = this.numLevels - symbolBlockBits;
    const { F, L, S } = this; // firsts, lasts, symbols
    F[0] = first;
    L[0] = last;
    S[0] = 0;
    const walk = new ReverseArrayWalker(lower === upper ? 0 : 1, F.length);
    let nRankCalls = 0;

    // In each iteration, we traverse F/L/S back to front and place the processed results at the end of the array,
    // filling it it in from right to left. This allows us to turn an individual element into more than one
    // processed element, for example if we want to recurse into both children of a node.
    // since we process from back to front, we also "go right" before we "go left" so that symbols retain their
    // left-to-right ordering when iterated in left-to-right order.
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBitMask = 1 << (this.maxLevel - l);
      for (let i = walk.len; i > 0; ) {
        i -= 1;

        const first = F[i];
        const first1 = level.rank1(first - 1);
        const first0 = first - first1;

        const last = L[i];
        const last1 = level.rank1(last - 1);
        const last0 = last - last1;

        const symbol = S[i];
        nRankCalls += 2;

        const num1 = last1 - first1; // count of right children
        if (num1 > 0) {
          // go right if the right node range [a, b] overlaps [lower, upper)
          const a = symbol | levelBitMask;
          const b = a | (levelBitMask - 1);
          const intervalsOverlap = lower <= b && a < upper;
          if (intervalsOverlap) {
            const nextIndex = walk.next();
            F[nextIndex] = nz + first1;
            L[nextIndex] = nz + last1;
            S[nextIndex] = a;
          }
        }

        const num0 = last0 - first0; // count of left children
        if (num0 > 0) {
          // go left if the left node range [a, b] overlaps [lower, upper)
          const a = symbol;
          const b = a | (levelBitMask - 1);
          const intervalsOverlap = lower <= b && a < upper;
          if (intervalsOverlap) {
            const nextIndex = walk.next();
            F[nextIndex] = first0;
            L[nextIndex] = last0;
            S[nextIndex] = symbol;
          }
        }
      }

      // update the length and move processed elements back to the front of the list.
      const index = walk.moveToFront();
      F.set(F.subarray(index, walk.cap));
      L.set(L.subarray(index, walk.cap));
      S.set(S.subarray(index, walk.cap));
    }
    for (let i = 0; i < walk.len; i++) L[i] -= F[i];
    const counts = L.subarray(0, walk.len).slice();
    const symbols = S.subarray(0, walk.len).slice();
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
    if (first > last) throw new Error('first must be <= last');
    if (last > this.length) throw new Error('last must be < wavelet matrix length');
    for (let i = 1; i < sortedIndices.length; i++) {
      if (!(sortedIndices[i - 1] <= sortedIndices[i])) throw new Error('sorted indices must be sorted');
    }
    // todo: error if there are more sortedindices than last-first, since we copy them into a scratch space (or ensure the space can hold the size we need)
    if (sortedIndices[0] < 0 || sortedIndices[sortedIndices.length - 1] >= last - first)
      throw new Error('sortedIndex cannot be less than zero or exceed length of range [first, last)');

    const { F, L, S, C } = this; // firsts, lasts, symbols, counts
    F[0] = first;
    L[0] = last;
    S[0] = 0;
    C[0] = sortedIndices.length; // number of sortedIndices represented by node 0
    // copy sorted indices into a scratch space since they are mutated as we go
    const I = this.C2.subarray(0, sortedIndices.length);
    I.set(sortedIndices);
    let nRankCalls = 0;

    const walk = new ReverseArrayWalker(sortedIndices.length === 0 ? 0 : 1, F.length);
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBitMask = 1 << (this.maxLevel - l);
      let k = I.length; // at every level, we sweep through the sorted indices in reverse
      for (let i = walk.len; i > 0; ) {
        i -= 1;

        const first = F[i];
        const first1 = level.rank1(first - 1);
        const first0 = first - first1;

        const last = L[i];
        const last1 = level.rank1(last - 1);
        const last0 = last - last1;

        let symbol = S[i];
        nRankCalls += 2;

        // Determine the number of nodes that wants to be mapped to the right child of this node,
        // then subtract the count of left children from all of the nodes matched to the right child,
        // to account for the elements counted in the left counted.
        const sortedIndexCount = C[i]; // number of sorted indices inside this node
        // [hi, lo) is the range of sorted indices covered by this node
        const lo = k - sortedIndexCount;
        const hi = k;

        const leftChildCount = last0 - first0; // left child count
        // index of the first right child
        // todo: avoid the binary search when all elems are mapped into the same child
        const splitIndex = binarySearchBefore(I, leftChildCount, lo, hi);
        const numGoLeft = splitIndex - lo;
        const numGoRight = sortedIndexCount - numGoLeft;

        // adjust count for quantiles mapped to the right child, taking into account the left count,
        // so that we look only for the remaining count of elements in the child node.
        for (let n = splitIndex; n < hi; n++) I[n] -= leftChildCount;

        k -= sortedIndexCount;

        if (numGoRight > 0) {
          // go right
          const nextIndex = walk.next();
          F[nextIndex] = nz + (first - first0); // = nz + first1
          L[nextIndex] = nz + (last - last0); // = nz + last1
          S[nextIndex] = symbol | levelBitMask;
          C[nextIndex] = numGoRight;
        }

        if (numGoLeft > 0) {
          // go left
          const nextIndex = walk.next();
          F[nextIndex] = first0;
          L[nextIndex] = last0;
          S[nextIndex] = symbol;
          C[nextIndex] = numGoLeft;
        }
      }

      // update the length and move processed elements back to the front of the list.
      const index = walk.moveToFront();
      F.set(F.subarray(index, walk.cap));
      L.set(L.subarray(index, walk.cap));
      S.set(S.subarray(index, walk.cap));
      C.set(C.subarray(index, walk.cap));
      // note: terminating this loop early computes approximate quantiles,
      // recursively dividing the alphabet in two each iteration.
      // the count of symbols assigned to each range is given by I,
      // and I think the ranges are [symbol, symbol+2^(maxLevel-l)).
    }

    for (let i = 0; i < walk.len; i++) L[i] -= F[i];
    const counts = L.subarray(0, walk.len).slice();
    const symbols = S.subarray(0, walk.len).slice();
    // assignments indicates how many entries in sortedIndices are assigned
    // to each symbol. this is a more economical representation than a dense array of
    // length sortedIndices when multiple sortedIndices point to the same symbol.
    const assignments = C.subarray(0, walk.len).slice();
    return { symbols, counts, assignments, nRankCalls };
  }

  quantiles(first, last, firstIndex, lastIndex) {
    // todo: for some reason  quantiles(first, last, 0, 0) returns a single value rather than nothing.
    if (first > last) throw new Error('first must be <= last');
    if (last > this.length) throw new Error('last must be < wavelet matrix length');
    if (firstIndex > lastIndex) throw new Error('firstIndex must be <= lastIndex');
    if (firstIndex < 0 || lastIndex > last - first)
      throw new Error('sortedIndex cannot be less than zero or exceed length of range [first, last)');
    const { F, L, S, C, C2 } = this; // firsts, lasts, symbols, counts
    F[0] = first;
    L[0] = last;
    S[0] = 0;
    C[0] = firstIndex;
    C2[0] = lastIndex;
    let nRankCalls = 0;

    const walk = new ReverseArrayWalker(firstIndex === lastIndex ? 0 : 1, F.length);
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBitMask = 1 << (this.maxLevel - l);
      for (let i = walk.len; i > 0; ) {
        i -= 1;

        const first = F[i];
        const first1 = level.rank1(first - 1);
        const first0 = first - first1;

        const last = L[i];
        const last1 = level.rank1(last - 1);
        const last0 = last - last1;

        let symbol = S[i];
        nRankCalls += 2;

        const leftChildCount = last0 - first0;
        if (C2[i] > leftChildCount) {
          // go right
          const nextIndex = walk.next();
          F[nextIndex] = nz + (first - first0); // = nz + first1
          L[nextIndex] = nz + (last - last0); // = nz + last1
          S[nextIndex] = symbol | levelBitMask;
          C[nextIndex] = leftChildCount < C[i] ? C[i] - leftChildCount : 0; // = max(C[i] - leftChildCount, 0);
          C2[nextIndex] = C2[i] - leftChildCount;
        }

        if (C[i] < leftChildCount) {
          // go left
          const nextIndex = walk.next();
          F[nextIndex] = first0;
          L[nextIndex] = last0;
          S[nextIndex] = symbol;
          C[nextIndex] = C[i];
          C2[nextIndex] = leftChildCount < C2[i] ? leftChildCount : C2[i]; // = min(leftChildCount, C2[i]);
        }
      }

      // update the length and move processed elements back to the front of the list.
      const index = walk.moveToFront();
      F.set(F.subarray(index, walk.cap));
      L.set(L.subarray(index, walk.cap));
      S.set(S.subarray(index, walk.cap));
      C.set(C.subarray(index, walk.cap));
      C2.set(C2.subarray(index, walk.cap));
    }

    for (let i = 0; i < walk.len; i++) L[i] -= F[i];
    const counts = L.subarray(0, walk.len).slice();
    const symbols = S.subarray(0, walk.len).slice();
    return { symbols, counts, nRankCalls };
  }

  // simple majority. todo: make configurable.
  majority(first, last) {
    const index = (last + first) >>> 1;
    const q = this.quantile(first, last, index);
    const count = last - first;
    const half = count >>> 1;
    if (q.count > half) return q.symbol;
    return null;
  }

  approxSizeInBits() {
    // ignores fixed-size fields
    let size = 0;
    for (const level of this.levels) {
      size += level.approxSizeInBits();
    }
    return size;
  }

  // For visualization
  symbolPath(i) {
    if (i < 0 || i > this.length) throw new Error('symbol: out of bounds');
    let l = 0; // level index
    let a = 0; // left symbol index
    let b = (1 << this.numLevels) - 1; // right symbol index
    const path = [];
    while (a !== b) {
      const level = this.levels[l];
      const m = (a + b) >>> 1;
      path.push({ index: i, bit: level.symbol(i) });
      if (level.symbol(i) === 0) {
        // go left
        i = level.rank0(i - 1);
        b = m;
      } else {
        // go right
        const nz = this.numZeros[l];
        i = nz + level.rank1(i - 1);
        a = m + 1;
      }
      l += 1;
    }
    return { symbol: a, path, virtualLeafIndex: i };
  }
}

// helper for turning recursion into iteration, encapsulating the logic
// of walking an array in reverse and generating zero or more outputs
// for each input.
// This is most useful when every element expands into at most k elements
// because then we can pre-allocate an array of size K times max elements
// and be sure that we will never end up overriding existing elements as
// we expand.
// new elements are filled in from the end of the array, using the space
// beyond the last element as scratch space.
// todo: next() could be nextWrite(), and we can add a nextRead(), and allow
// flipping the elements from back to front incrementally, rather than
// memcpying to the front after every iteration (with a final pass to
// copy to the front if needed).
class ReverseArrayWalker {
  constructor(len, cap) {
    this.len = len; // length taken by elements
    this.cap = cap; // capacity for additional elements
    this.index = cap; // nextIndex + 1
  }
  next() {
    // return the next index at which we can append an element
    // (as we fill the array in backwards from arr[cap - 1])
    return (this.index -= 1);
  }
  moveToFront() {
    // logic for moving the filled-in elements from the end
    // to the front of the array, returning the start index
    // that can be used to copy elements from the back to front:
    // > const index = walk.moveToFront()
    // > typedArray.set(typedArray.subarray(index), walk.cap);
    const index = this.index;
    this.len = this.cap - this.index;
    this.index = this.cap;
    return index;
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
