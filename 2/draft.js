// note: incomplete
class __draft__WaveletTree {
  // This implements a specialized version of 'constructWaveletRepresentation' for wavelet trees that does less work per level.
  // It can keep around the symbol offsets in the 'virtual bottom row', which are the cumulative sum values in `borders`.
  // todo: construct borders in-place reusing the `hist` array, which we do not need once we have the borders.
  constructor(data, alphabetSize) {
    console.clear();
    // todo: return early if no data; otherwise maxLevel is -1
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

    this.alphabetSize = alphabetSize;
    this.maxLevel = maxLevel;
    this.borders = borders;
    this.levels = levels;
  }

  // Adapted from Compact Data Structures: A Practical Approach (Algorithm 6.3)
  // The implementation is incomplete.
  // function waveletTreeRank(tree, c, i) {
  //   if (tree.length === 0) return 0;
  //   let l = 0; // level index
  //   let a = 0;
  //   let b = tree.alphabetSize - 1;
  //   let left = 0;
  //   let right = tree.levels[0].length;
  //   let levelBitMask = 1 << (levels.length - 1);
  //   while (a !== b) {
  //     const level = tree[l];
  //     if ((c & levelBitMask) === 0) {
  //       // go left
  //       i = level.rank0(i);
  //       // original: v = v.l
  //       // [todo] ours: right = ...
  //       b = (a + b) >>> 1;
  //     } else {
  //       // go right
  //       i = z + level.rank1(i);
  //       // original: v = v.r
  //       // [todo] ours: left = ...
  //       a = ((a + b) >>> 1) + 1;
  //     }
  //     l += 1;
  //     levelBitMask >>>= 1;
  //   }
  //   return i;
  // }
}
