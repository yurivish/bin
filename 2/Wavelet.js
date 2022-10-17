import { RankBitVector } from './RankBitVector';

// Adapted from Compact Data Structures: A Practical Approach (Algorithm 6.3)
function waveletTreeRank(tree, c, i) {
  if (tree.length === 0) return 0;
  let l = 0; // level index
  let a = 0;
  let b = tree.alphabetSize - 1;
  let left = 0;
  let right = tree.levels[0].length;
  let levelBitMask = 1 << (levels.length - 1);
  while (a !== b) {
    const level = tree[l];
    if ((c & levelBitMask) === 0) {
      // go left
      i = level.rank0(i);
      // v = v.l
      // right = level.rank1(
      b = (a + b) >>> 1;
    } else {
      // go right
      i = z + level.rank1(i);
      // v = v.r
      a = ((a + b) >>> 1) + 1;
    }
    l += 1;
    levelBitMask >>>= 1;
  }
  return i - p;
}

// Adapted from Compact Data Structures: A Practical Approach (Algorithm 6.6)
// This implementation has not been tested and is here for possible future use.
function waveletMatrixRank(tree, c, i) {
  if (tree.length === 0) return 0;
  let p = 0;
  let l = 0; // level index
  let a = 0;
  let b = tree.alphabetSize - 1;
  let levelBitMask = 1 << (levels.length - 1);
  while (a !== b) {
    const level = tree[l];
    if ((c & levelBitMask) === 0) {
      // go left
      i = level.rank0(i);
      p = level.rank0(p);
      b = (a + b) >>> 1;
    } else {
      // go right
      const z = tree.z[l];
      i = z + level.rank1(i);
      p = z + level.rank1(p);
      a = ((a + b) >>> 1) + 1;
    }
    l += 1;
    levelBitMask >>>= 1;
  }
  return i - p;
}

// todo: can we include a signalling bitvector for the start of each symbol,
// to allow us to effectively update the final layer indicating whether a symbol is filtered out?
//
// This implements Algorithm 1 (seq.pc) from the paper "Practical Wavelet Tree Construction".
// The return value is the array of wavelet tree levels. Adapting the algorithm to construct
// a wavelet matrix instead requires changing the borders computation (see section 5.3):
// "We compute the starting positions of the intervals using the bit-reversal permutation,
// i.e., the lines mentioned above are changed to Borders[ρl[i]] = Borders[ρl[i − 1]] + Hist[ρl[i − 1]].
// Then, the resulting starting positions of the intervals for bit prefixes are in bit reversal
// permutation order, i.e., the starting positions of the intervals for a wavelet matrix."
// todo: constructWaveletRepresentation [tree/matrix];
export function constructWaveletTree(data, alphabetSize) {
  // console.clear();
  // data is an integer array [0, alphabetSize)
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
  // return hist;
  // Construct the other levels of the wavelet tree bottom-up
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
    // todo: if statement to pick the tree/matrix for loop
    // todo: do we need to do this, or can we instead just sample
    // borders[2i] (or borders[2 << levelBit * i]) since it contains
    // the cumulative sum, meaning that this work is not needed?
    // (all the desired values already exist, just spaced out)
    // could be: 1. compute histogram 2. compute cumulative sum
    // then just use the cumulative sum every time, spacing out
    // accesses as needed (see slides and diagrams from Paul Dinklage's
    // "Translating Between Wavelet Tree and Wavelet Matrix Construction":
    // http://www.stringology.org/event/2019/p12.html
    for (let i = 1; i < m; i++) {
      // Update the positions in-place
      borders[i] = borders[i - 1] + hist[i - 1];
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
      const symbolIndex = (d & bitPrefixMask) >>> bitPrefixShift;
      const p = borders[symbolIndex];
      borders[symbolIndex] += 1;
      // Set the bit in the bitvector
      if (d & levelBitMask) level.one(p);
    }
  }
  for (let i = 0; i < numLevels; i++) {
    levels[i].finish();
  }
  levels.alphabetSize = alphabetSize;
  return levels;
}
