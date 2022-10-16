import { RankBitVector } from './RankBitVector';

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
export function constructWaveletTree(data, numSymbols) {
  console.clear();
  // data is an integer array [0, numSymbols)
  const n = data.length;
  const numLevels = Math.ceil(Math.log2(numSymbols));
  const maxLevel = numLevels - 1;
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
    console.log('');
    console.log('level', l);
    const m = 2 ** l;
    // Compute the histogram based on the previous level's one
    for (let i = 0; i < m; i++) {
      // Update the histogram in-place
      hist[i] = hist[2 * i] + hist[2 * i + 1];
    }
    // Get starting positions of intervals from the new histogram
    borders[0] = 0;
    // todo: if statement to pick the tree/matrix for loop
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
  return levels;
}
