import { RankBitVector } from './RankBitVector';

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
    // this.numNodes = numNodes;
    this.maxLevel = maxLevel;
    this.length = data.length;

    this.symbols = new Uint32Array(this.alphabetSize);
    for (let s = 0; s < this.symbols.length; s++) this.symbols[s] = s;
    this.P = new Uint32Array(this.alphabetSize + 1); // extra space for the full 'both' path
    this.I = new Uint32Array(this.alphabetSize + 1); // case and an odd alphabet size
  }

  access(i) {
    if (i < 0 || i > this.length) throw new Error('access: out of bounds');
    let l = 0; // level index
    let a = 0; // left symbol index
    let b = (1 << this.numLevels) - 1; // right symbol index
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
    if (symbol >= this.alphabetSize) throw new Error('symbol must be < alphabetSize');
    let p = 0; // index of the start of the current node
    let l = 0; // level index
    let a = 0; // left symbol index
    let b = (1 << this.numLevels) - 1; // right symbol index
    let levelBitMask = 1 << this.maxLevel;
    i = clamp(i + 1, 1, this.length);
    while (a !== b) {
      const level = this.levels[l];
      const m = (a + b) >>> 1;
      const i0 = level.rank0(i - 1);
      const p0 = level.rank0(p - 1);
      if ((symbol & levelBitMask) === 0) {
        // go left
        i = i0;
        p = p0;
        b = m;
      } else {
        // go right
        const nz = this.numZeros[l];
        i = nz + (i - i0); // === nz + level.rank1(i - 1);
        p = nz + (p - p0); // === nz + level.rank1(p - 1);
        a = m + 1;
      }
      l += 1;
      levelBitMask >>>= 1;
    }
    return i - p;
  }

  // todo: remove, but not before adding more comments to `ranks`.
  // Batched rank: https://www.sciencedirect.com/science/article/pii/S0890540112001526#se0100
  // Returns the rank of all the symbols. I think this is now superseded by the "ranks" function.
  // This version rearranges the computation to make level bitvector accesses contiguous in time.
  // Since we proceed from high bits to low bits, we save on rank queries when a symbol has the
  // same level bit as its predecessor. Currently done with an `if`; could instead be a nested
  // loop where the inner loop iterates all of the contiguous symbols with the same bit.
  // note: this comment is out of date, but the strategy it describes could be useful for arbitrary symbol sets.
  // we should also try implementing a batched rank over a contiguous symbol range that returns individual counts.
  batchedRank(i) {
    // todo: call this allRanks? ranks?
    const { symbols, P, I } = this;
    // P.fill(123);
    // I.fill(123); // clear for easier debugging
    P[0] = 0;
    I[0] = clamp(i + 1, 1, this.length);
    let levelBitMask = 1 << this.maxLevel;
    let numLevelNodes = 1; // tracks the number of branching paths; 2 * numLevelNodes paths at the current level (2 * leaves at this level).
    const numLeafNodes = this.alphabetSize >>> 1; // don't go beyond the last symbol when len(symbols) is not a power of 2
    let prevI = 0; // we want to store I[numLeafNodes] for the level maxLevels - 1 so that
    let prevP = 0; // we have it around to compute the final rank value if needed.

    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      prevI = I[numLeafNodes];
      prevP = P[numLeafNodes];
      // Perform all left and right mappings to the next level of the tree.
      // reach level, we go both left and right for each existing entry in
      // the array (when we come to a fork in the road, we take it both directions).
      // We started out with a single index i and p; at level zero we want
      // to expand go both left and right; then at level one we want to go
      // both left and right for the level zero "lefts", and same for the rights;
      // this way we don't have to perform O(len(symbols)) rank ops at each level
      // and can instead perform O(level) rank ops where level is O(log2(len(symbols))).
      // In the end, I think we perform just under 2 * len(symbols) - 1 rank operations
      // since we do two per tree node and eg. an 8-symbol tree has 4 + 2 + 1 = 7 nodes.
      // additionally, all of the rank operations at a level are done in a row.
      // todo: better explanation, docs, and cleaner + commented code.
      // todo: can we generalize this to arbitrary symbol sets? is it already general
      // (based on this.symbols)? I think we need to adjust numLeafNodes to the size of the symbol set.
      // all levels other than the last are completely filled if we assume that the number of levels is
      // the closest power of two to the number of symbols, but otherwise any level may benefit from
      // the reduction in the number of iterations.
      for (let n = Math.min(numLeafNodes, numLevelNodes); n > 0; ) {
        n -= 1;

        const i = I[n];
        const i0 = level.rank0(i - 1);
        I[2 * n] = i0;
        I[2 * n + 1] = nz + (i - i0);

        const p = P[n];
        const p0 = level.rank0(p - 1);
        P[2 * n] = p0;
        P[2 * n + 1] = nz + (p - p0);

        // question: how can we bail out early in the case of an unbalanced tree?
        // question: how can we use similar logic to query arbitrary symbol ranges?
        // question: can we use similar logic to query arbitrary symbol sets?
      }
      levelBitMask >>>= 1;
      numLevelNodes <<= 1;
    }

    // If the alphabet size is odd, fill in the rank of the final symbol since
    // it was not filled in the previous loop, which operates two symbols at a time.
    if (this.alphabetSize % 2 === 1) {
      const level = this.levels[this.maxLevel];
      I[2 * numLeafNodes] = level.rank0(prevI - 1);
      P[2 * numLeafNodes] = level.rank0(prevP - 1);
    }

    for (let i = 0; i < I.length; i++) I[i] -= P[i];
    return I.slice(0, this.alphabetSize);
  }

  ranks(i, bitSelectors) {
    // note: bitSelectors could be a u64 for up to 32 levels (2 bits per selector)
    if (bitSelectors.length != this.numLevels) throw new Error('bitSelectors.length must be equal to numLevels');
    let len = 1;
    // important: round up. This means that for odd alphabet sizes,
    // we will computer an extra element if we went 'both' directions,
    // which will be omitted from the return value with `subarray`.
    const numSymbols = this.alphabetSize;
    const halfLimit = Math.ceil(numSymbols / 2);

    const { P, I } = this;
    // P.fill(123);
    // I.fill(123); // clear for easier debugging

    I[0] = clamp(i + 1, 1, this.length);
    P[0] = 0;

    let levelBitMask = 1 << this.maxLevel;
    loop: for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const selector = bitSelectors[l];
      switch (selector) {
        case BOTH:
          for (let n = Math.min(len, halfLimit); n > 0; ) {
            n -= 1;
            const i = I[n];
            const p = P[n];

            const i1 = level.rank1(i - 1);
            I[2 * n + 1] = nz + i1; // go right
            I[2 * n] = i - i1; // go left (=== level.rank0(i - 1))

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

  // Adapted from https://github.com/noshi91/Library/blob/0db552066eaf8655e0f3a4ae523dbf8c9af5299a/data_structure/wavelet_matrix.cpp#L76
  // Range quantile query returning the kth largest symbol in A[i, j).
  // I wonder if there's a way to early-out in this implementation; as
  // written, it always looks through all levels. Does this have implications
  // for e.g. a Huffman-shaped wavelet matrix?
  // Note: this  may be noticeably slower (needs rigorous testing) than the previous
  // version that looped over all levels rather than bisecting a symbol interval
  // when the tree is balanced. Might be worth keeping both implementations around.
  quantile(i, j, k) {
    if (i > j) throw new Error('i must be <= j');
    if (j > this.length) throw new Error('j must be < wavelet matrix length');
    if (k < 0 || k >= j - i) throw new Error('k cannot be less than zero or exceed length of range [i, j)');
    let symbol = 0;
    let l = 0;
    let a = 0; // left symbol index
    let b = (1 << this.numLevels) - 1; // right symbol index
    let levelBitMask = 1 << this.maxLevel;
    while (a !== b) {
      const level = this.levels[l];
      const m = (a + b) >>> 1;
      const i0 = level.rank0(i - 1);
      const j0 = level.rank0(j - 1);
      const count = j0 - i0;
      if (k < count) {
        // go left
        i = i0;
        j = j0;
        b = m;
      } else {
        symbol |= levelBitMask;
        k -= count;
        const nz = this.numZeros[l];
        i = nz + (i - i0); // === nz + level.rank1(i - 1);
        j = nz + (j - j0); // === nz + level.rank1(j - 1);
        a = m + 1;
      }
      l += 1;
      levelBitMask >>>= 1;
    }
    return { symbol, frequency: j - i };
  }

  // todo: describe
  less(i, j, symbol) {
    if (i < 0) throw new Error('i must be >= 0');
    if (i > j) throw new Error('i must be <= j');
    if (j > this.length) throw new Error('j must be < wavelet matrix length');
    if (symbol <= 0) return 0;
    if (symbol >= this.alphabetSize) return this.length;
    let l = 0; // level index
    let a = 0; // left symbol index
    let b = (1 << this.numLevels) - 1; // right symbol index // right symbol
    let count = 0;
    let levelBitMask = 1 << this.maxLevel;
    while (a !== b) {
      const level = this.levels[l];
      const m = (a + b) >>> 1;
      const i0 = level.rank0(i - 1);
      const j0 = level.rank0(j - 1);
      if ((symbol & levelBitMask) === 0) {
        i = i0;
        j = j0;
        b = m;
      } else {
        count += j0 - i0;
        const nz = this.numZeros[l];
        i = nz + (i - i0); // === nz + level.rank1(i - 1);
        j = nz + (j - j0); // === nz + level.rank1(j - 1);
        a = m + 1;
      }
      l += 1;
      levelBitMask >>>= 1;
    }
    return count;
  }

  __quantile(i, j, k) {
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

  // note:
  // we can express rank1 in terms of rank0
  // when i and j are both in bounds:
  //    rank0(i)   = i - rank1(i) + 1
  // => rank0(i-1) = i-1 - rank1(i-1) + 1
  // => rank0(i-1) = i - rank1(i-1)
  // and vice versa (see impl. of rank0)
  
  // Returns the number of values with symbol strictly less than the given symbol.
  // This is kind of like a ranged rank operation over a symbol range.
  // Adapted from https://github.com/noshi91/Library/blob/0db552066eaf8655e0f3a4ae523dbf8c9af5299a/data_structure/wavelet_matrix.cpp#L25
  __less(i, j, symbol) {
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

function clamp(x, lo, hi) {
  return x < lo ? lo : x > hi ? hi : x;
}
