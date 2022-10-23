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
    // we add one because bit-selector based ranks might compute an extra symbol
    const sz = Math.min(2 ** this.numLevels, this.alphabetSize + 1);
    this.F = new Uint32Array(sz); // firsts
    this.L = new Uint32Array(sz); // lasts
    this.S = new Uint32Array(sz); // symbols
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

  count(first, last, symbol) {
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
    let n = 0
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const first1 = level.rank1(first - 1);
      const last1 = level.rank1(last - 1); // inclusive
      n += 2;
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
    console.log('less:', n,'calls to rank')
    return count
  }

  // Returns the number of occurrences of symbols [lower, upper)
  // in the index range [first, last). It is possible to implement
  // this function roughly twice as efficiently in terms of number 
  // of rank calls, but this comes at a cost of implementation complexity.
  // See the paper "New algorithms on wavelet trees and applications to
  // information retrieval" for details. Another approach is to modify the
  // implementation of rankLess to perform two interleaved calls. The 
  // subtlety there is that, as written, the algorithm does not work when
  // symbol >= alphabetSize (the one-symbol impl. can return early in this case).
  countRange(first, last, lower, upper) {
    return this.rankLess(first, last, upper) - this.rankLess(first, last, lower);
  }

  // Returns all of the distinct symbols [lower, upper) in the range [first, last)
  // together with their number of occurrences. If symbol block bits is specified,
  // then symbols are grouped together when they differ only in their `symbolBlockBits`
  // lowest bits. Each distinct group is labeled by its lowest element, which represents
  // the group containing symbols in the range [symbol, symbol + 2^symbolBlockBits).
  counts(first, last, lower, upper, symbolBlockBits = 0) {
    const symbolBlockSize = (1 << symbolBlockBits)
    // these error messages could be improved, explaining that ignore bits tells us the power of two
    // that lower and upper need to be multiples of.
    if (lower % symbolBlockSize !== 0) throw new Error('lower must evenly divide the symbol block size implied by symbolBlockBits') 
    if (upper  % symbolBlockSize !== 0) throw new Error('upper must evenly divide the symbol block size implied by symbolBlockBits') 
    const numLevels = this.numLevels - symbolBlockBits
    // if (upper - lower < ) throw new Error('step size implied by symbolBlockBits is greater than the specified symbol range (results would be misleading)')
    const { F, L, S } = this; // firsts, lasts, symbols
    // F.fill(123); // for debugging
    // L.fill(123);
    // S.fill(123);
    F[0] = first;
    L[0] = last;
    S[0] = 0;
    let len = 1;

    // In each iteration, we traverse F/L/S back to front and place the processed results at the end of the array,
    // filling it it in from right to left. This allows us to turn an individual element into more than one
    // processed element, for example if we want to recurse into both children of a node.
    // since we process from back to front, we also "go right" before we "go left" so that symbols retain their
    // left-to-right ordering when iterated in left-to-right order.
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const levelBitMask = 1 << (this.maxLevel - l);
      let nextIndex = F.length - 1;
      for (let i = len; i > 0; ) {
        i -= 1;

        const first = F[i];
        const first1 = level.rank1(first - 1);
        const first0 = first - first1;
        const last = L[i];
        const last1 = level.rank1(last - 1);
        const last0 = last - last1;
        const symbol = S[i];
        const num1 = last1 - first1; // count of right children
        if (num1 > 0) {
          // go right if the right node range (a, b) overlaps [lower, upper)
          const a = symbol | levelBitMask;
          const b = a | (levelBitMask - 1);
          const intervalsOverlap = lower <= b && a < upper;
          if (intervalsOverlap) {
            const nz = this.numZeros[l];
            F[nextIndex] = nz + first1;
            L[nextIndex] = nz + last1;
            S[nextIndex] = a;
            nextIndex -= 1;
          }
        }

        const num0 = last0 - first0; // count of left children
        if (num0 > 0) {
          // go left if the left node range (a, b) overlaps [lower, upper)
          const a = symbol;
          const b = a | (levelBitMask - 1);
          const intervalsOverlap = lower <= b && a < upper;
          if (intervalsOverlap) {
            F[nextIndex] = first0;
            L[nextIndex] = last0;
            S[nextIndex] = symbol;
            nextIndex -= 1;
          }
        }
      }

      // update the length and move processed elements back to the front of the list.
      len = F.length - (nextIndex + 1);
      F.set(F.subarray(nextIndex + 1));
      L.set(L.subarray(nextIndex + 1));
      S.set(S.subarray(nextIndex + 1));
    }
    for (let i = 0; i < len; i++) L[i] -= F[i];
    const counts = L.subarray(0, len).slice();
    const symbols = S.subarray(0, len).slice();
    return { symbols, counts };
  }

  // Adapted from https://github.com/noshi91/Library/blob/0db552066eaf8655e0f3a4ae523dbf8c9af5299a/data_structure/wavelet_matrix.cpp#L76
  // Range quantile query returning the kth largest symbol in A[i, j).
  quantile(first, last, sortedIndex) {
    if (first > last) throw new Error('first must be <= last');
    if (last > this.length) throw new Error('last must be < wavelet matrix length');
    if (sortedIndex < 0 || sortedIndex >= last - first)
      throw new Error('sortedIndex cannot be less than zero or exceed length of range [first, last)');
    let symbol = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const first0 = level.rank0(first - 1);
      const last0 = level.rank0(last - 1);
      const count = last0 - first0;
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
    return { symbol, count: last - first };
  }

  // note: there is an extended version of the quantile algorithm in the paper
  // "New algorithms on wavelet trees and applications to information retrieval"
  // It calls the alg above rqq (Algorithm 3) and the extended version mrqq (Algorithm 5).
  // it's structurally the same, except there are now potentially multiple return values
  // so we would need to use the same reverse-the-list trick as in ranksRange to manage
  // an unpredictable number of nodes.

  // i think this is the bigger performance analysis for the rank of all individual symbols with bit selectors:
  // σ        = numLevels
  // k <= 2^σ = alphabetSize
  // k - 1    = numNodes = number of rank queries for all symbols rank via ranks()
  // 2σ * k   = number of rank queries for all symbols rank via rank();


  // todo: include symbolset functionality to aid constructing ranks queries:
  // https://observablehq.com/d/05a4c693328c3c34

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
