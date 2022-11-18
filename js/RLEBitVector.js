class RLEBitVector {
  constructor() {
    this.Z = [];
    this.ZO = [];
    this.length = 0;
  }

  // Encodes a run of `numZeros` zeros followed by `numOnes` ones
  run(numZeros, numOnes) {
    const len = this.Z.length;

    // Append the cumulative number of zeros to the Z array
    const prevZ = len > 0 ? this.Z[len - 1] : 0;
    this.Z.push(prevZ + numZeros);

    // Append the index of the first one of this run to the ZO array
    this.ZO.push(this.length + numZeros);

    // Update the bitvector length
    this.length += numZeros + numOnes;
  }

  // Encodes a run of ones starting at at index `i`
  oneRun(i, length) {
    if (i < this.length) throw new Error('oneRun cannot overlap any pre-existing runs');
    // Number of zeros preceding this 1-run
    const numZeros = i - this.length;
    const numOnes = length;
    this.run(numZeros, numOnes);
  }

  rank1(i) {
    if (i < 0) return 0;
    if (i >= this.length) return this.length - this.Z[this.Z.length - 1];

    // Number of complete zero-runs up to virtual index i
    const j = sparseRank1(this.ZO, i);

    // Number of zeros preceding the current 1-run
    let numPrecedingZeros = j === 0 ? 0 : sparseSelect1(this.Z, j);

    // Ensure that `i` does not point beyond the last 1 in the current 1-run
    if (j < this.Z.length) {
      // Number of zeros in the next (0, 1)-run
      const numZerosNextRun = sparseSelect1(this.Z, j + 1) - numPrecedingZeros;
      // Virtual index of the start of the next 1-run
      const nextOneRunIndex = sparseSelect1(this.ZO, j + 1);
      // Virtual index of the end of the current 1-run
      const lastOneRunIndex = nextOneRunIndex - numZerosNextRun - 1;
      if (i > lastOneRunIndex) i = lastOneRunIndex;
    }
    return i - numPrecedingZeros + 1;
  }

  rank0(i) {
    if (i < 0) return 0;
    if (i >= this.length) return this.Z[this.Z.length - 1];
    return i - this.rank1(i) + 1;
  }

  select0(i) {
    if (i < 1 || i > this.Z[this.Z.length - 1]) return -1;

    // The i-th zero is in the j-th (0, 1)-block.
    const j = sparseRank1(this.Z, i - 1);

    // If we're in the first (0, 1)-block, the i-th zero is at index i - 1.
    if (j === 0) return i - 1;

    // Number of zeros preceding the jth (0, 1)-block
    const numPrecedingZeros = sparseSelect1(this.Z, j);

    // Number of zeros up to and including the jth (0, 1)-block
    const cumulativeZeros = sparseSelect1(this.Z, j + 1);

    // Number of zeros in the j-th (0, 1) block
    const numZeros = cumulativeZeros - numPrecedingZeros;

    // Index of the first 1 in the j-th (0, 1)-block
    const firstOneIndex = sparseSelect1(this.ZO, j + 1);

    // Index of the first zero of the j-th (0, 1)-block
    const firstZeroIndex = firstOneIndex - numZeros;

    // Return the index of the (i - numPrecedingZeros)th zero in the j-th (0, 1)-block.
    return firstZeroIndex + (i - numPrecedingZeros) - 1;
  }

  select1(i) {
    if (i < 1 || i > this.length - this.Z[this.Z.length - 1]) return -1;

    // The i-th one is in the j-th (0, 1)-block.
    const j =
      binarySearchAfter((k) => sparseSelect1(this.ZO, k + 1) - sparseSelect1(this.Z, k + 1), i - 1, 0, this.Z.length) -
      1;

    // Index of the first 1 in the j-th (0, 1)-block
    const firstOneIndex = sparseSelect1(this.ZO, j + 1);

    // Number of zeros up to and including the jth (0, 1)-block
    const cumulativeZeros = sparseSelect1(this.Z, j + 1);

    // Number of ones in the blocks preceding the j-th block (ZO[j+1] - Z[j+1])
    const numPrecedingOnes = firstOneIndex - cumulativeZeros;

    // Return the index of the (i - numPrecedingOnes)th one in the j-th (0, 1)-block.
    return firstOneIndex + (i - numPrecedingOnes) - 1;
  }

  function finish() {
    this.Z = new Uint32Array(this.Z)
    this.ZO = new Uint32Array(this.ZO)
  }
}

// Returns the rightmost insertion index for T in A in order to maintain A's sorted order.
// Searches the index range [L, R).
function binarySearchAfter(access, T, L, R) {
  while (L < R) {
    // Note: This midpoint calculation will return incorrect results for arrays with length > 2^30
    // Correct (but slow) alternative: Math.floor(L + (R - L) / 2)
    const m = (L + R) >>> 1;
    if (access(m) > T) R = m;
    else L = m + 1;
  }
  return R;
}

// Rank1 on a sparse bitvector represented by an array of 1-bit positions
function sparseRank1(bv, i) {
  return binarySearchAfter((i) => bv[i], i, 0, bv.length);
}

// Select on a sparse bitvector represented by an array of 1-bit positions
function sparseSelect1(bv, i) {
  return bv[i - 1];
}
