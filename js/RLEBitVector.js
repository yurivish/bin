export class RLEBitVector {
  constructor() {
    this.Z = [];
    this.ZO = [];
    this.length = 0;
    this.numZeros = 0;
    this.numOnes = 0;
  }

  // Encodes a run of `numZeros` zeros followed by `numOnes` ones
  run(numZeros, numOnes) {
    if (numZeros === 0 && numOnes === 0) return;
    const len = this.Z.length;
    this.numZeros += numZeros;
    this.numOnes += numOnes;
    this.length += numZeros + numOnes;
    if (numZeros === 0 && len > 0) {
      // This run consists of only ones; coalesce it with the
      // previous run (since all runs contain ones at their end).
      this.ZO[len - 1] += numOnes;
    } else if (numOnes === 0 && this.lastBlockContainsOnlyZeros()) {
      // This run consists of only zeros; coalesce it with the
      // previous run (since it turns out to consist of only zeros).
      this.Z[len - 1] += numZeros;
      this.ZO[len - 1] += numZeros;
    } else {
      // No coalescing is possible; create a new block of runs.
      // Append the cumulative number of zeros to the Z array
      this.Z.push(this.numZeros);
      // Append the index of the last one of this run to the ZOE array
      this.ZO.push(this.length);
    }
  }

  lastBlockContainsOnlyZeros() {
    const len = this.Z.length;
    if (len === 0) return false;
    if (len === 1) return this.Z[0] === this.ZO[0];
    const lastBlockLength = this.ZO[len - 1] - this.ZO[len - 2];
    const lastBlockNumZeros = this.Z[len - 1] - this.Z[len - 2];
    return lastBlockLength === lastBlockNumZeros;
  }

  finish() {
    // todo: only convert to typed arrays if their value is < 2^32
    this.Z = new Uint32Array(this.Z);
    this.ZO = new Uint32Array(this.ZO);
  }

  rank0(i) {
    if (i < 0) return 0;
    if (i >= this.length) return this.numZeros;
    return i - this.rank1(i) + 1;
  }

  rank1(i) {
    if (i < 0) return 0;
    if (i >= this.length) return this.numOnes;

    // Number of complete 01-runs up to and including virtual index i
    const j = sparseRank1(this.ZO, i);

    // Number of zeros including the j-th block
    const numCumulativeZeros = sparseSelect1(this.Z, j + 1);

    // Number of zeros preceding the j-th block
    const numPrecedingZeros = j === 0 ? 0 : sparseSelect1(this.Z, j);

    // Number of zeros in the j-th block
    const numZeros = numCumulativeZeros - numPrecedingZeros;

    // Start index of the j-th block
    const blockStart = j === 0 ? 0 : sparseSelect1(this.ZO, j);

    // Number of ones preceding the j-th block
    const numPrecedingOnes = blockStart - numPrecedingZeros;

    // Start index of ones in the j-th block
    const onesStart = blockStart + numZeros;

    return numPrecedingOnes + Math.max(0, i - onesStart + 1);
  }

  access(i) {
    // Quick hack. Can do better.
    return this.rank1(i) - this.rank1(i - 1)
  }

  alignedRank0(i) {
    if (i < 0) return 0;
    if (i >= this.length) return this.numZeros;

    // Number of complete 01-runs up to virtual index i
    const j = sparseRank1(this.ZO, i);

    // Number of zeros preceding the (aligned) index i
    return sparseSelect1(this.Z, j + 1);
  }

  alignedRank1(i) {
    if (i < 0) return 0;
    if (i >= this.length) return this.numOnes;
    return i - this.alignedRank0(i) + 1;
  }

  select0(i) {
    if (i < 1 || i > this.numZeros) return -1;

    // The i-th zero is in the j-th 01-block.
    const j = sparseRank1(this.Z, i - 1);

    // If we're in the first 01-block, the i-th zero is at index i - 1.
    if (j === 0) return i - 1;

    // Start index of the j-th block
    const blockStart = sparseSelect1(this.ZO, j);

    // Number of zeros preceding the jth block
    const numPrecedingZeros = sparseSelect1(this.Z, j);

    // Return the index of the (i - numPrecedingZeros)th zero in the j-th 01-block.
    return blockStart + (i - numPrecedingZeros) - 1;
  }

  select1(i) {
    if (i < 1 || i > this.numOnes) return -1;

    // The i-th one is in the j-th 01-block.
    const j = binarySearchAfter(
      (k) => sparseSelect1(this.ZO, k + 1) - sparseSelect1(this.Z, k + 1),
      i - 1,
      0,
      this.Z.length,
    );

    // Number of zeros up to and including the jth block
    const numCumulativeZeros = sparseSelect1(this.Z, j + 1);

    return numCumulativeZeros + i - 1;
  }

  approxSizeInBits() {
    // ignores fixed-size fields
    const zBits = 8 * this.Z.length * this.Z.BYTES_PER_ELEMENT;
    const zoBits = 8 * this.ZO.length * this.ZO.BYTES_PER_ELEMENT;
    return zBits + zoBits;
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

// Returns the rightmost insertion index for T in A in order to maintain A's sorted order.
// Searches the index range [L, R).
function binarySearchAfterArray(A, T, L, R) {
  while (L < R) {
    // Note: This midpoint calculation will return incorrect results for arrays with length > 2^30
    // Correct (but slow) alternative: Math.floor(L + (R - L) / 2)
    const m = (L + R) >>> 1;
    if (A[m] > T) R = m;
    else L = m + 1;
  }
  return R;
}


// Rank1 on a sparse bitvector represented by an array of 1-bit positions
function sparseRank1(bv, i) {
  return binarySearchAfter((i) => bv[i], i, 0, bv.length);
}

// Select on a sparse bitvector represented by an array of 1-bit positions
function sparseSelect1(bv, n) {
  return n < 1 || n > bv.length ? -1 : bv[n - 1];
}