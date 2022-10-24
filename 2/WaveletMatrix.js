import { RankBitVector } from './RankBitVector';
import { ZeroCompressedBitVector } from './ZeroCompressedBitVector';
import { reverseBits, reverseBits32, clamp } from './util';
// todo: range next value, range prev value (though these can be done using quantile), quantiles, majority,
// intersection and more from "New algorithms on wavelet trees and applications to information retrieval"

// wavelet tree is from 2003, wavelet matrix is from oct. 2012
// to do: hoist error objects to the top?

// note: implements a binary wavelet matrix that always splits on power-of-two
// alphabet boundaries, rather than splitting based on the true alphabet midpoint.
export class WaveletMatrix {
  // This implements Algorithm 1 (seq.pc) from the paper "Practical Wavelet Tree Construction".
  // The return value is the array of wavelet tree levels. Adapting the algorithm to construct
  // a wavelet matrix instead requires changing the borders computation (see section 5.3).
  // todo: check that all symbols are < alphabetSize
  constructor(data, alphabetSize) {
    // console.clear();
    // data is an array of integer values in [0, alphabetSize)
    const n = data.length;
    const numLevels = Math.ceil(Math.log2(alphabetSize));
    const maxLevel = numLevels - 1;

    // for each level, store the number of nodes at that level..
    // for a balanced binary tree, there are twice as many nodes
    // at each child level as its parent level, other than the
    // final level, which is limited by the alphabet size.
    // (the previous levels are also limited by the alphabet size,
    // but the way that we determine the number of levels means
    // at the only power of 2 that may be greater than alphabetSize
    // is the last one.
    // This will look different with Huffman-shaped matrices!
    // const numLeafNodes = Math.ceil(this.alphabetSize / 2); // ! there are alphabetSize 'virtual' leaves
    // const numNodes = new Uint32Array(numLevels);
    // numNodes[0] = 1;
    // for (let i=1;i<maxLevel;i++) numNodes[i] = 2 * numNodes[i-1];
    // numNodes[maxLevel] = numLeafNodes;

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
    const level0 = levels[0];
    const highBitMask = 1 << maxLevel;
    for (let i = 0; i < n; i++) {
      const d = data[i];
      // Compute the histogram of the data
      hist[d] += 1;
      // Fill the first level's bit vector (MSBs in data order)
      if (d & highBitMask) level0.one(i);
    }
    // Construct the other levels bottom-up
    for (let l = maxLevel; l > 0; l--) {
      // console.log('');
      // console.log('level', l);
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
      // console.log("levelBit:", levelBit);
      // console.log("levelBitMask:", (levelBitMask >>> 0).toString(2));
      // console.log("bitPrefixMask:", (bitPrefixMask >>> 0).toString(2));
      for (let i = 0; i < n; i++) {
        const d = data[i];
        // Get and update position for bit by computing its bit prefix,
        // which encodes the path from the root to the node at level l
        // containing this bit
        const nodeIndex = (d & bitPrefixMask) >>> bitPrefixShift;
        // console.assert(nodeIndex < m);
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
    // this.numNodes = numNodes;
    this.maxLevel = maxLevel;
    this.length = data.length;

    this.symbols = new Uint32Array(this.alphabetSize);
    for (let s = 0; s < this.symbols.length; s++) this.symbols[s] = s;

    // todo: do not materialize these until needed - alphabet might be big.
    // accept scratch space as an input parameter so we can reuse the same
    // space across wavelet trees
    const sz = Math.min(2 ** this.numLevels, this.alphabetSize);
    this.F = new Uint32Array(sz); // firsts
    this.L = new Uint32Array(sz); // lasts
    this.S = new Uint32Array(sz); // symbols
    this.I = new Uint32Array(sz); // indices
  }

  access(index) {
    if (index < 0 || index > this.length) throw new Error('access: out of bounds');
    let symbol = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const index1 = level.rank1(index - 1); // todo: combine rank and access queries since they access the same block
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

  countSymbol(first, last, symbol) {
    if (symbol >= this.alphabetSize) throw new Error('symbol must be < alphabetSize');
    if (first > last) throw new Error('last must be <= first');
    if (first === last) return 0;
    if (first > this.length) throw new Error('first must be < wavelet matrix length');
    if (this.numLevels === 0) return 0;
    for (let l = 0; l < this.numLevels; l++) {
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
  countLess(first, last, symbol) {
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
  // implementation of countLess to perform two interleaved calls. The
  // subtlety there is that, as written, the algorithm does not work when
  // symbol >= alphabetSize (the one-symbol impl. can return early in this case).
  count(first, last, lower, upper) {
    return this.countLess(first, last, upper) - this.countLess(first, last, lower);
  }

  // Returns all of the distinct symbols [lower, upper) in the range [first, last)
  // together with their number of occurrences. If symbol block bits is specified,
  // then symbols are grouped together when they differ only in their `groupByLowBits`
  // lowest bits. Each distinct group is labeled by its lowest element, which represents
  // the group containing symbols in the range [symbol, symbol + 2^groupByLowBits).
  counts(first, last, lower, upper, groupByLowBits = 0) {
    const symbolBlockSize = 1 << groupByLowBits;
    // these error messages could be improved, explaining that ignore bits tells us the power of two
    // that lower and upper need to be multiples of.
    if (lower % symbolBlockSize !== 0)
      throw new Error('lower must evenly divide the symbol block size implied by groupByLowBits');
    if (upper % symbolBlockSize !== 0)
      throw new Error('upper must evenly divide the symbol block size implied by groupByLowBits');
    const numLevels = this.numLevels - groupByLowBits;
    // if (upper - lower < ) throw new Error('step size implied by groupByLowBits is greater than the specified symbol range (results would be misleading)')
    const { F, L, S } = this; // firsts, lasts, symbols
    // F.fill(123); // for debugging
    // L.fill(123);
    // S.fill(123);
    F[0] = first;
    L[0] = last;
    S[0] = 0;
    const walk = new ReverseArrayWalker(1, F.length);
    let nRankCalls = 0;

    // In each iteration, we traverse F/L/S back to front and place the processed results at the end of the array,
    // filling it it in from right to left. This allows us to turn an individual element into more than one
    // processed element, for example if we want to recurse into both children of a node.
    // since we process from back to front, we also "go right" before we "go left" so that symbols retain their
    // left-to-right ordering when iterated in left-to-right order.
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
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
            const nz = this.numZeros[l];
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
  quantile(first, last, sortedIndex) {
    if (first > last) throw new Error('first must be <= last');
    if (last > this.length) throw new Error('last must be < wavelet matrix length');
    if (sortedIndex < 0 || sortedIndex >= last - first)
      throw new Error('sortedIndex cannot be less than zero or exceed length of range [first, last)');
    let symbol = 0;
    let nRankCalls = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const first0 = level.rank0(first - 1);
      const last0 = level.rank0(last - 1);
      const count = last0 - first0;
      nRankCalls += 2;
      if (sortedIndex < count) {
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
        sortedIndex -= count;
      }
    }
    return { symbol, count: last - first, nRankCalls };
  }

  quantiles(first, last, sortedIndices) {
    if (first > last) throw new Error('first must be <= last');
    if (last > this.length) throw new Error('last must be < wavelet matrix length');
    for (let i = 1; i < sortedIndices.length; i++) {
      if (!(sortedIndices[i - 1] <= sortedIndices[i])) throw new Error('sorted indices must be sorted');
    }
    if (sortedIndices[0] < 0 || sortedIndices[sortedIndices.length - 1] >= last - first)
      throw new Error('sortedIndex cannot be less than zero or exceed length of range [first, last)');

    const { F, L, S, I } = this; // firsts, lasts, symbols, counts
    F[0] = first;
    L[0] = last;
    S[0] = 0;
    I[0] = sortedIndices.length;
    let nRankCalls = 0;

    const walk = new ReverseArrayWalker(1, F.length);
    let symbol = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const levelBitMask = 1 << (this.maxLevel - l);
      let k = sortedIndices.length;
      for (let i = walk.len; i > 0; ) {
        i -= 1;

        const first = F[i];
        const first1 = level.rank1(first - 1);
        const first0 = first - first1;

        const last = L[i];
        const last1 = level.rank1(last - 1);
        const last0 = last - last1;

        let symbol = S[i];
        const count = last0 - first0;
        nRankCalls += 2;

        let numGoLeft = 0;
        let numGoRight = 0;

        // assign target quantiles to the left or right child of this node
        k -= I[i];
        for (let n = I[i]; n > 0; ) {
          n -= 1; // iterates in reverse over the range [0, I[i]-1]
          const index = k + n;
          // sortedIndices[index] belongs to this tree node
          if (sortedIndices[index] < count) {
            numGoLeft++;
          } else {
            sortedIndices[index] -= count;
            numGoRight++;
          }
        }

        if (numGoRight > 0) {
          // go right
          const nz = this.numZeros[l];
          const nextIndex = walk.next();
          F[nextIndex] = nz + (first - first0); // = nz + first1
          L[nextIndex] = nz + (last - last0); // = nz + last1
          S[nextIndex] = symbol | levelBitMask;
          I[nextIndex] = numGoRight;
        }

        if (numGoLeft > 0) {
          // go left
          const nextIndex = walk.next();
          F[nextIndex] = first0;
          L[nextIndex] = last0;
          S[nextIndex] = symbol;
          I[nextIndex] = numGoLeft;
        }
      }

      // update the length and move processed elements back to the front of the list.
      const index = walk.moveToFront();
      F.set(F.subarray(index, walk.cap));
      L.set(L.subarray(index, walk.cap));
      S.set(S.subarray(index, walk.cap));
      I.set(I.subarray(index, walk.cap));
      // note: terminating this loop early computes approximate quantiles,
      // recursively dividing the alphabet in two each iteration.
      // the count of symbols assigned to each range is given by I,
      // and I think the ranges are [symbol, symbol+2^(maxLevel-l)).
    }

    for (let i = 0; i < walk.len; i++) L[i] -= F[i];
    const counts = L.subarray(0, walk.len).slice();
    const symbols = S.subarray(0, walk.len).slice();
    // assignments indicates how many quantiles are assigned to each symbol.
    // this is a more concise representation than a dense array with one entry
    // per sortedIndex in those cases where the same symbol is at multiple input indices.
    const assignments = I.subarray(0, walk.len).slice();
    return { symbols, counts, assignments, nRankCalls };
  }

  // note:
  // we can express rank1 in terms of rank0
  // when i and j are both in bounds:
  //    rank0(i)   = i - rank1(i) + 1
  // => rank0(i-1) = i-1 - rank1(i-1) + 1
  // => rank0(i-1) = i - rank1(i-1)
  // and vice versa, so i0 = i - i1; (see impl. of rank0)

  approxSizeInBits() {
    // ignores fixed-size fields
    let size = 0;
    for (const level of this.levels) {
      size += level.approxSizeInBits();
    }
    return size;
  }

  // For visualization
  accessPath(i) {
    if (i < 0 || i > this.length) throw new Error('access: out of bounds');
    let l = 0; // level index
    let a = 0; // left symbol index
    let b = (1 << this.numLevels) - 1; // right symbol index
    const path = [];
    while (a !== b) {
      const level = this.levels[l];
      path.push({ index: i, bit: level.access(i) });
      if (level.access(i) === 0) {
        // go left
        i = level.rank0(i - 1);
        b = (a + b) >>> 1;
      } else {
        // go right
        const nz = this.numZeros[l];
        i = nz + level.rank1(i - 1);
        a = ((a + b) >>> 1) + 1;
      }
      l += 1;
    }
    return { symbol: a, path, virtualLeafIndex: i };
  }
}

// helper for turning recursion into iteration, encapsulating the logic
// of walking an array in reverse and generating zero or more outputs
// for each input.
// new elements are filled in from the end of the array, using the space
// beyond then last element as scratch space.
// in our case we never generate more than two outputs per input, which
// upper-bounds the size of the scratch space required to 2x the number
// of input elements.
class ReverseArrayWalker {
  constructor(len, cap) {
    this.len = len; // length
    this.cap = cap; // capacity
    this.index = cap; // nextIndex + 1
  }
  next() {
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
