import { binarySearchAfter, binarySearchAfterAccess } from './util.js';
import { SparseOneBitVector } from './SparseBitVector.js';

export class RLEBitVector {
  constructor() {
    this.Z = [];
    this.ZO = [];
    this.length = 0;
    this.numZeros = 0;
    this.numOnes = 0;
    this.z = this.zo = null;
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
      // Append the cumulative number of ones and zeros to the ZO array
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
    const { Z, ZO, length, numOnes } = this;
    this.z = new SparseOneBitVector(length, this.Z).finish();
    this.zo = new SparseOneBitVector(length, this.ZO).finish();
    this.Z = this.ZO = null; // drop references to the original buffers
    return this;
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
    const j = this.zo.rank1(i);

    // Number of zeros including the j-th block
    const numCumulativeZeros = this.z.select1(j + 1);

    // Number of zeros preceding the j-th block
    const numPrecedingZeros = j === 0 ? 0 : this.z.select1(j);

    // Number of zeros in the j-th block
    const numZeros = numCumulativeZeros - numPrecedingZeros;

    // Start index of the j-th block
    const blockStart = j === 0 ? 0 : this.zo.select1(j);

    // Number of ones preceding the j-th block
    const numPrecedingOnes = blockStart - numPrecedingZeros;

    // Start index of ones in the j-th block
    const onesStart = blockStart + numZeros;

    return numPrecedingOnes + Math.max(0, i - onesStart + 1);
  }

  // Aligned variants are guaranteed to be correct *only* at the end of
  // each 0-run and 1-run. Results will be returned for other indices,
  // but may be wildly incorrect – even negative. So be very careful.
  alignedRank0(i) {
    if (i < 0) return 0;
    if (i >= this.length) return this.numZeros;

    // Number of complete 01-runs up to virtual index i
    const j = this.zo.rank1(i);

    // Number of zeros preceding the (aligned) index i
    return this.z.select1(j + 1);
  }

  alignedRank1(i) {
    if (i < 0) return 0;
    if (i >= this.length) return this.numOnes;
    return i - this.alignedRank0(i) + 1;
  }

  access(i) {
    // Quick hack. Can do better. Will need to account for flip.
    return this.rank1(i) - this.rank1(i - 1);
  }

  select0(n) {
    if (n < 1 || n > this.numZeros) return -1;

    // The i-th zero is in the j-th 01-block.
    const j = this.z.rank1(n - 1);

    // If we're in the first 01-block, the i-th zero is at index i - 1.
    if (j === 0) return n - 1;

    // Start index of the j-th block
    const blockStart = this.zo.select1(j);

    // Number of zeros preceding the jth block
    const numPrecedingZeros = this.z.select1(j);

    // Return the index of the (i - numPrecedingZeros)th zero in the j-th 01-block.
    return blockStart + (n - numPrecedingZeros) - 1;
  }

  select1(n) {
    if (n < 1 || n > this.numOnes) return -1;

    // The n-th one is in the j-th 01-block.
    const j = binarySearchAfterAccess((i) => this.zo.select1(i + 1) - this.z.select1(i + 1), n - 1, 0, this.z.numOnes);

    // Number of zeros up to and including the jth blocka
    const numCumulativeZeros = this.z.select1(j + 1);

    return numCumulativeZeros + n - 1;
  }

  approxSizeInBits() {
    // ignores fixed-size fields
    if (this.z) return this.z.approxSizeInBits() + this.zo.approxSizeInBits();
    return this.sparse.approxSizeInBits();
  }
}
