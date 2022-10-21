import { RankBitVector } from './RankBitVector';
import { ZeroCompressedBitVector } from './ZeroCompressedBitVector';

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

  // todo: consistency w/ count of whether the symbol specifier comes first or last

  // Adapted from Compact Data Structures: A Practical Approach (Algorithm 6.6)
  // Searches [first, last).
  // This implements the 'strict' version, using only this.levels and this.numZeros.
  // I noticed that by setting the initial value of P, we can do rank range queries rather than requiring two calls to rank.
  rank(symbol, first, last) {
    if (this.numLevels === 0) return 0;
    if (symbol >= this.alphabetSize) throw new Error('symbol must be < alphabetSize');
    if (first > last) throw new Error('last must be <= first');
    if (first === last) return 0;
    if (first > this.length) throw new Error('first must be < wavelet matrix length');
    let l = 0; // level index
    let a = 0; // left symbol index
    let b = (1 << this.numLevels) - 1; // right symbol index
    let levelBitMask = 1 << this.maxLevel;
    while (a !== b) {
      const level = this.levels[l];
      const m = (a + b) >>> 1;
      const last0 = level.rank0(last - 1);
      const first0 = level.rank0(first - 1); // todo: so strange that we need to do rank of a negative number when first is 0...
      if ((symbol & levelBitMask) === 0) {
        // go left
        last = last0;
        first = first0;
        b = m;
      } else {
        // go right
        const nz = this.numZeros[l];
        last = nz + (last - last0); // === nz + level.rank1(last - 1);
        first = nz + (first - first0); // === nz + level.rank1(first - 1);
        a = m + 1;
      }
      l += 1;
      levelBitMask >>>= 1;
    }
    return last - first;
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
    // question: do we ever go less than numLevels iterations?
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

  // Returns the number of values with symbol strictly less than the given symbol.
  // This is kind of like a ranged rank operation over a symbol range.
  // Adapted from https://github.com/noshi91/Library/blob/0db552066eaf8655e0f3a4ae523dbf8c9af5299a/data_structure/wavelet_matrix.cpp#L25
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

  // is less related to the idea of the 'last' value in rank?
  // ie. on the last (virtual) level, we want the position of
  // the first instance of symbol S; that position IS the rank
  // if we have not restricted the time range.
  // so less is kinda the same as the ranks function (but slower)?
  // could we have a class of bit selectors that allow you to choose
  // successive prefixes of the code space, so that when we STOP we've
  // got all the LESS values that we want? do we already have this implemented,
  // in fact? if we calculate contiguous symbol range sums, then just do a 
  // cumulative sum on that... voila?
  // in short:
  // - `ranks` allows us to compute sums of contiguous power of 2 symbol ranges
  // - if we cumsum those sums, we get the less-than values for 
  //   evenly-spaced power of 2 symbols, eg. sum for symbols < 4,
  //   symbols < 8, symbols < 12, symbols < 16.
  // todo: does the halfRange stuff still work with STOP sum ranges?
  // seems to.
  // so it's not possible to get the kind of power of two code ranges 
  // from `ranks`that we use BEFORE range splitting. But I totally wonder
  // if the set of bitselector-based contiguous sequences isn't exactly
  // the same as what you get out of range-splitting...
  // Plus will need one less query to get all the symbols below the lowest
  // symbol covered by the `ranks` range 

  // todo: include symbolset functionality to aid constructing ranks queries:
  // https://observablehq.com/d/05a4c693328c3c34

  // Restore the cleaner version that iterates through all levels each time.
  // I think with balance trees this is always the case, at least for certain operations.
  // if we ever have unbalanced trees using prefix free codes, I think this means that
  // the codes will be strictly ascending in terms of the number of bits they use.
  // which means that we can use a look up table for the cumulative number of codes below
  // a certain bit length that is of length #bits (store #bits => code offset) or such. so
  // based on the symbol code, we can tell in O(numLevels) the number of levels to traverse.
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

  // note:
  // we can express rank1 in terms of rank0
  // when i and j are both in bounds:
  //    rank0(i)   = i - rank1(i) + 1
  // => rank0(i-1) = i-1 - rank1(i-1) + 1
  // => rank0(i-1) = i - rank1(i-1)
  // and vice versa (see impl. of rank0)

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
}

// Reverse the lowest `numBits` bits of `v`.
// E.g. bitReverse(0b0000100100, 6) === 0b0000001001
//                       ^^^^^^               ^^^^^^
function bitReverse(v, numBits) {
  return bitReverse32(v) >>> (32 - numBits);
}

// Adapted from https://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel
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
