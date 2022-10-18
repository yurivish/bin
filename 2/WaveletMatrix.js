import { RankBitVector } from './RankBitVector';

// note: implements a binary wavelet matrix that always splits on power-of-two
// alphabet boundaries, rather than splitting based on the true alphabet midpoint.
export class WaveletMatrix {
  // This implements Algorithm 1 (seq.pc) from the paper "Practical Wavelet Tree Construction".
  // The return value is the array of wavelet tree levels. Adapting the algorithm to construct
  // a wavelet matrix instead requires changing the borders computation (see section 5.3).
  constructor(data, alphabetSize) {
    // console.clear();
    // data is an array of integer values in [0, alphabetSize)
    const n = data.length;
    const numLevels = Math.ceil(Math.log2(alphabetSize));
    const maxLevel = numLevels - 1;
    // todo: can we get away with non-pow2, storing just one entry per symbol?
    const hist = new Uint32Array(2 ** numLevels);
    const borders = new Uint32Array(2 ** numLevels);
    const levels = new Array(numLevels);
    // Initialize the level bit vectors
    for (let i = 0; i < numLevels; i++) {
      levels[i] = new RankBitVector(data.length);
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
        const prevIndex = bitReverse(i - 1, l);
        borders[bitReverse(i, l)] = borders[prevIndex] + hist[prevIndex];
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
    this.maxLevel = maxLevel;
    this.length = data.length;
  }

  access(i) {
    if (i < 0 || i > this.length) throw new Error('access: out of bounds');
    let l = 0; // level index
    let a = 0;
    let b = (1 << this.numLevels) - 1;
    while (a !== b) {
      const level = this.levels[l];
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
    return a;
  }

  // Adapted from Compact Data Structures: A Practical Approach (Algorithm 6.6)
  // This implements the 'strict' version, using only this.levels and this.numZeros.
  rank(symbol, i) {
    if (this.numLevels === 0) return 0;
    i += 1;
    let p = 0; // index of the start of the current node
    let l = 0; // level index
    let a = 0;
    let b = 2 ** this.numLevels - 1;
    let levelBitMask = 1 << this.maxLevel;
    while (a !== b) {
      const level = this.levels[l];
      if ((symbol & levelBitMask) === 0) {
        // go left
        i = level.rank0(i - 1);
        p = level.rank0(p - 1);
        b = (a + b) >>> 1;
      } else {
        // go right
        const nz = this.numZeros[l];
        i = nz + level.rank1(i - 1);
        p = nz + level.rank1(p - 1);
        a = ((a + b) >>> 1) + 1;
      }
      l += 1;
      levelBitMask >>>= 1;
    }
    return i - p;
  }

  // Adapted from https://github.com/noshi91/Library/blob/0db552066eaf8655e0f3a4ae523dbf8c9af5299a/data_structure/wavelet_matrix.cpp#L76
  // Range quantile query returning the kth largest symbol in A[i, j).
  // I wonder if there's a way to early-out in this implementation; as
  // written, it always looks through all levels. Does this have implications
  // for e.g. a Huffman-shaped wavelet matrix?
  // Note: Compact Data Structures: A Practical Approach contains a
  // pseudocode implementation for Wavelet Trees as Algorithm 6.16.
  quantile(i, j, k) {
    if (i > j) throw new Error('i must be <= j');
    if (j > this.length) throw new Error('j must be < wavelet matrix length');
    if (k < 0 || k >= j - i) throw new Error('k cannot be less than zero or exceed length of range [i, j)');
    const msbMask = 1 << this.maxLevel;
    let symbol = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const i0 = level.rank0(i - 1);
      const j0 = level.rank0(j - 1);
      const count = j0 - i0;
      if (k < count) {
        i = i0;
        j = j0;
      } else {
        symbol |= msbMask >>> l;
        k -= count;
        const nz = this.numZeros[l];
        i = nz + (i - i0); // === nz + level.rank1(i - 1);
        j = nz + (j - j0); // === nz + level.rank1(j - 1);
      }
    }
    return { symbol, frequency: j - i };
  }

  // Returns the number of values with symbol strictly less than the given symbol.
  // This is kind of like a ranged rank operation over a symbol range.
  // Adapted from https://github.com/noshi91/Library/blob/0db552066eaf8655e0f3a4ae523dbf8c9af5299a/data_structure/wavelet_matrix.cpp#L25
  less(i, j, symbol) {
    if (i < 0) throw new Error('i must be >= 0');
    if (i > j) throw new Error('i must be <= j');
    if (j > this.length) throw new Error('j must be < wavelet matrix length');
    if (symbol <= 0) return 0;
    if (symbol >= this.alphabetSize) return this.length;
    let levelBitMask = 1 << this.maxLevel;
    let count = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const i0 = level.rank0(i - 1);
      const j0 = level.rank0(j - 1);
      if ((symbol & levelBitMask) === 0) {
        i = i0;
        j = j0;
      } else {
        count += j0 - i0;
        const nz = this.numZeros[l];
        // we can simplifiy this way since i and j are both
        // in range and rank0(i) === i - this.rank1(i) + 1
        i = nz + (i - i0); // === nz + level.rank1(i - 1);
        j = nz + (j - j0); // === nz + level.rank1(j - 1);
      }
      levelBitMask >>>= 1;
    }
    return count;
  }

  // Returns the number of occurrences of symbols [lower, upper)
  // in the index range [i, j).
  count(i, j, lower, upper) {
    if (lower > upper) throw new Error('lower must be <= upper');
    return this.less(i, j, upper) - this.less(i, j, lower);
  }
}

// Reverse the lowest `numBits` bits of `v`.
// E.g. bitReverse(0b0000100100, 6) === 0b0000001001
//                       ^^^^^^               ^^^^^^
function bitReverse(v, numBits) {
  return bitReverse32(v) >>> (32 - numBits);
}

// Adapted from https://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel
// Follow this operation by a >>> shift to
function bitReverse32(v) {
  // unsigned int v; // 32-bit word to reverse bit order
  // swap odd and even bits
  v = ((v >>> 1) & 0x55555555) | ((v & 0x55555555) << 1);
  // swap consecutive pairs
  v = ((v >>> 2) & 0x33333333) | ((v & 0x33333333) << 2);
  // swap nibbles ...
  v = ((v >>> 4) & 0x0f0f0f0f) | ((v & 0x0f0f0f0f) << 4);
  // swap bytes
  v = ((v >>> 8) & 0x00ff00ff) | ((v & 0x00ff00ff) << 8);
  // swap 2-byte long pairs
  v = (v >>> 16) | (v << 16);
  return v;
}