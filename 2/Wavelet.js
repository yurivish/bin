import { RankBitVector } from './RankBitVector';

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
  console.clear();
  // data is an integer array [0, alphabetSize)
  const n = data.length;
  const numLevels = Math.ceil(Math.log2(alphabetSize));
  const maxLevel = numLevels - 1;
  // todo: can we get away with non-pow2, storing just one entry per symbol?
  const hist = new Uint32Array(2 ** numLevels);
  const borders = new Uint32Array(2 ** numLevels);
  const adjustments = new Uint32Array(2 ** numLevels);
  const levels = new Array(numLevels);
  // Initialize the level bit vectors
  for (let i = 0; i < numLevels; i++) {
    levels[i] = new RankBitVector(data.length);
  }
  const level0 = levels[0];
  const highBitMask = 1 << maxLevel;
  console.log('numLevels', numLevels);

  for (let i = 0; i < n; i++) {
    const d = data[i];
    // Compute the histogram of the data
    hist[d] += 1;
    // hist[bitReverse(d, numLevels)] += 1;
    // Fill the first level's bit vector (MSBs in data order)
    if (d & highBitMask) level0.one(i);
  }

  for (let i = 1; i < hist.length; i++) {
    // Calculate cumulative sum
    borders[i] = borders[i - 1] + hist[i - 1];
    // console.log(i, bitReverse(i, numLevels))
    // borders[bitReverse(i, numLevels)] = borders[bitReverse(i - 1, numLevels)] + hist[bitReverse(i - 1, numLevels)];
  }
  console.log('hist', hist);
  console.log('borders', borders);

  // Construct the other levels of the wavelet tree bottom-up
  for (let l = maxLevel; l > 0; l--) {
    console.log('');
    console.log('level', l);

    const m = 2 ** l;
    // Fill the bit vector of the current level
    const level = levels[l];
    const levelBit = maxLevel - l;
    const levelBitMask = 1 << levelBit;
    const bitPrefixMask = 0xfffffffe << levelBit;
    const bitPrefixShift = levelBit + 1;
    console.log('levelBit:', levelBit);
    console.log('levelBitMask:', (levelBitMask >>> 0).toString(2));
    console.log('bitPrefixMask:', (bitPrefixMask >>> 0).toString(2));
    const numPriorNodes = 2 ** l - 1;
    console.log('numPriorNodes:', numPriorNodes);
    adjustments.subarray(0, 2 ** l).fill(0); // todo: can we avoid zeroing this an extra time?
    const hopSize = (2 ** numLevels) >>> l;
    console.log('hopSize:', hopSize);
    for (let i = 0; i < n; i++) {
      const d = data[i];
      // Get and update position for bit by computing its bit prefix,
      // which encodes the path from the root to the node at level l
      // containing this bit
      const levelNodeIndex = (d & bitPrefixMask) >>> bitPrefixShift; // index of this node along the nodes in this level
      const nodeIndex = hopSize * levelNodeIndex; // index of this node in breadth-first ordering across the full tree
      const p = borders[nodeIndex] + adjustments[levelNodeIndex]; // index of this symbol in its level bitvector
      adjustments[levelNodeIndex] += 1;
      // Set the bit in the bitvector
      if (d & levelBitMask) level.one(p);
    }
  }
  for (let i = 0; i < numLevels; i++) levels[i].finish();
  levels.alphabetSize = alphabetSize;
  levels.borders = borders;
  return levels;
}

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

// Adapted from Compact Data Structures: A Practical Approach (Algorithm 6.3)
// The implementation is incomplete.
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
      // original: v = v.l
      // [todo] ours: right = ...
      b = (a + b) >>> 1;
    } else {
      // go right
      i = z + level.rank1(i);
      // original: v = v.r
      // [todo] ours: left = ...
      a = ((a + b) >>> 1) + 1;
    }
    l += 1;
    levelBitMask >>>= 1;
  }
  return i;
}
