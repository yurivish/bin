import { RankBitVector } from './RankBitVector';
import { ZeroCompressedBitVector } from './ZeroCompressedBitVector';
import { reverseBits, reverseBits32, clamp } from './util';
// todo: range next value, range prev value (though these can be done using quantile), quantiles, majority,
// intersection and more from "New algorithms on wavelet trees and applications to information retrieval"

// wavelet tree is from 2003, wavelet matrix is from oct. 2012

// bit selectors
export const LEFT = 0; // select left bit
export const RIGHT = 1; // select right bit
export const BOTH = 2; // select both bits (eg. rank: return ranks for both subtrees)
export const STOP = 3; // stop at this level (eg. rank: return sums so far)

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
    this.P = new Uint32Array(this.alphabetSize + 1); // extra space for the full 'both' path
    this.I = new Uint32Array(this.alphabetSize + 1); // case and an odd alphabet size
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

  rank(first, last, symbol) {
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
        first = nz + first - first0; // = nz + first1
        last = nz + last - last0; // = nz + last1
        // update symbol and new target sorted index in the child node
        const levelBitMask = 1 << (this.maxLevel - l);
        symbol |= levelBitMask;
        sortedIndex -= count;
      }
    }
    return { symbol, count: last - first };
  }

  // note: less and rank look *very* similar - can we compute the less values 'for (almost) free' here,
  // by also adding to each element's count when we go right? probably, but is it useful...
  // and i think when we STOP, we want to add to each symbol the count
  // as if we went right for all subsequent levels.
  // is this related to numZeros at those lower levels? need to be careful...
  // could make a lesses function, and another that does both, if we ever need both...
  // is there a way to do range less in one pass? maybe just track 2 counts or something
  ranks(selectors, first, last) {
    if (first === undefined || last === undefined) throw 'wat';
    if (first > last) throw new Error('last must be <= first');
    // if (first === last) return 0; // todo: return a 0 array of the number of returned symbols
    if (first > this.length) throw new Error('first must be < wavelet matrix length');
    // note: bit selectors could be a u64 for up to 32 levels (2 bits per selector)
    if (selectors.length != this.numLevels) throw new Error('selectors.length must be equal to numLevels');
    let len = 1; // number of symbols we're currently tracking / updating
    // important: round up. This means that for odd alphabet sizes,
    // we will computer an extra element if we went 'both' directions,
    // which will be omitted from the return value with `subarray`.
    const numSymbols = this.alphabetSize;
    const halfLimit = Math.ceil(numSymbols / 2);

    const { P, I } = this;
    // P.fill(123);
    // I.fill(123); // clear for easier debugging

    I[0] = last; // clamp(last + 1, 1, this.length);
    P[0] = first; // todo: would starting this off at nonzero allow us to do a 'range count'? or do I need to do 2 'ranks' calls after all

    let levelBitMask = 1 << this.maxLevel;
    loop: for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const selector = selectors[l];
      switch (selector) {
        case BOTH:
          for (let n = Math.min(len, halfLimit); n > 0; ) {
            n -= 1;
            const last = I[n];
            const p = P[n];

            const last1 = level.rank1(last - 1);
            I[2 * n + 1] = nz + last1; // go right
            I[2 * n] = last - last1; // go left (=== level.rank0(last - 1))

            const p1 = level.rank1(p - 1);
            P[2 * n + 1] = nz + p1; // go right
            P[2 * n] = p - p1; // go left (=== level.rank0(p - 1))
          }
          len = Math.min(2 * len, numSymbols);
          break;
        case LEFT:
          for (let n = 0; n < len; n++) {
            // go left
            I[n] = level.rank0(I[n] - 1);
            P[n] = level.rank0(P[n] - 1);
          }
          break;
        case RIGHT:
          for (let n = 0; n < len; n++) {
            // go right
            I[n] = nz + level.rank1(I[n] - 1);
            P[n] = nz + level.rank1(P[n] - 1);
          }
          break;
        case STOP:
          break loop;
      }
      levelBitMask >>>= 1;
    }
    for (let i = 0; i < len; i++) I[i] -= P[i];
    return I.subarray(0, len).slice();
  }

  // i think this is the bigger performance analysis for the rank of all individual symbols:
  // σ        = numLevels
  // k <= 2^σ = alphabetSize
  // k - 1    = numNodes = number of rank queries for all symbols rank via ranks()
  // 2σ * k   = number of rank queries for all symbols rank via rank();


  // Returns the number of values with symbol strictly less than the given symbol.
  // This is kind of like a ranged rank operation over a symbol range.
  // Adapted from https://github.com/noshi91/Library/blob/0db552066eaf8655e0f3a4ae523dbf8c9af5299a/data_structure/wavelet_matrix.cpp#L25
  rankLess(first, last, symbol) {
    if (first < 0) throw new Error('first must be >= 0');
    if (first > last) throw new Error('first must be <= last');
    if (last > this.length) throw new Error('last must be < wavelet matrix length');
    if (symbol <= 0) return 0;
    if (symbol >= this.alphabetSize) return this.length;
    let count = 0;
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
        // update count before going right
        count += (last - last1) - (first - first1); // = last0 - first0
        // go right
        const nz = this.numZeros[l];
        first = nz + first1;
        last = nz + last1;
      }
    }
    return count;
  }


  // todo: include symbolset functionality to aid constructing ranks queries:
  // https://observablehq.com/d/05a4c693328c3c34

  // note:
  // we can express rank1 in terms of rank0
  // when i and j are both in bounds:
  //    rank0(i)   = i - rank1(i) + 1
  // => rank0(i-1) = i-1 - rank1(i-1) + 1
  // => rank0(i-1) = i - rank1(i-1)
  // and vice versa, so i0 = i - i1; (see impl. of rank0)

  // Returns the number of occurrences of symbols [lower, upper)
  // in the index range [i, j).
  count(i, j, lower, upper) {
    if (lower > upper) throw new Error('lower must be <= upper');
    return this.less(i, j, upper) - this.less(i, j, lower);
  }

  allSelector() {
    return new Uint8Array(this.numLevels).fill(BOTH);
  }

  selector(prefix) {
    const ret = new Uint8Array(this.numLevels).fill(STOP);
    ret.set(prefix);
    return ret;
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
