import { binarySearchAfter, binarySearchAfterAccess } from './util.js';
import { SparseBitVector } from './SparseBitVector.js';

export class RLEBitVector {
  constructor() {
    this.Z = [];
    this.ZO = [];
    this.length = 0;
    this.numZeros = 0;
    this.numOnes = 0;
    this.z = this.zo = this.sparseOnes = this.sparseZeros = null;
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

  // Returns a nonempty info if all 1-runs are length 1 (or length 0, in the case of the final 01-run).
  // Otherwise returns null. As a special case, we return null if this.length is zero.
  isSparse() {
    if (this.length === 0) return null;
    let pz = 0;
    let pzo = 0;
    let numOnes = 0;
    let numZeros = 0;
    let oneSparse = true;
    let zeroSparse = true;
    for (let i = 0; i < this.Z.length; i++) {
      let z = this.Z[i];
      let zo = this.ZO[i];
      let numZeros = z - pz;
      let numZerosAndOnes = zo - pzo;
      numOnes = numZerosAndOnes - numZeros;
      if (numOnes > 1) oneSparse = false;
      if (numZeros > 1) zeroSparse = false;
      if (!oneSparse && !zeroSparse) return null;
      pz = z;
      pzo = zo;
    }
    return {
      oneSparse,
      zeroSparse,
      hasTrailingZeros: this.length > 0 && numOnes === 0
    };
  }

  finish() {
    const { Z, ZO, length, numOnes } = this;

    // "Sparse" is not quite the right word; we mean that all 0- or 1-runs are length one.
    const { oneSparse, zeroSparse, hasTrailingZeros } = this.isSparse() ?? {};
    if (oneSparse) {
      // All 1-runs are length one, so we can use a single sparse bitvector to represent them.
      for (let i = 0; i < ZO.length; i++) ZO[i] -= 1;
      // construct after the loop since the type depends on the maximum element.
      // If there are trailing zeros, exclude the final element from the array
      // type calculation since it won't be used.
      const lastIndex = ZO.length - 1 - hasTrailingZeros;
      let data = new (typedArrayType(ZO[lastIndex]))(this.ZO);
      if (hasTrailingZeros) data = data.subarray(0, -1);
      this.sparseOnes = new SparseBitVector(data, length);
    } else if (zeroSparse) {
      // All 0-runs are length one, so we can use a single sparse bitvector to represent them.
      const leadingZero = length > 0 && Z[0] > 0;
      const lastIndex = ZO.length - 1;
      // If there is a leading zero, we want store ZO into data offset by 1
      // so that data[0] === 0 and data[data.length-1]  === ZO[ZO.length-2].
      let data = new (typedArrayType(ZO[lastIndex]))(ZO);
      if (leadingZero) {
        data.set(data.subarray(0, -1), 1);
        data[0] = 0;
      } else {
        data = data.subarray(0, -1);
      }
      this.sparseZeros = new SparseBitVector(data, length);
    } else {
      // Turn each of Z and ZO into a u32 or f64 array depending on the maximum value
      // and construct sparse bitvectors for both.
      const lastIndex = Z.length - 1;
      this.z = new SparseBitVector(
        new (typedArrayType(Z[lastIndex]))(this.Z),
        this.Z.length
      );
      this.zo = new SparseBitVector(
        new (typedArrayType(ZO[lastIndex]))(this.ZO),
        this.ZO.length
      );
    }
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
    if (this.sparseOnes) return this.sparseOnes.rank1(i);
    if (this.sparseZeros) return this.sparseZeros.rank0(i);
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
    if (this.sparseOnes) return this.sparseOnes.rank0(i);
    if (this.sparseZeros) return this.sparseZeros.rank1(i);

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
    // Quick hack. Can do better.
    return this.rank1(i) - this.rank1(i - 1);
  }

  select0(n) {
    if (n < 1 || n > this.numZeros) return -1;
    if (this.sparseOnes) return this.sparseOnes.select0(n);
    if (this.sparseZeros) return this.sparseZeros.select1(n);

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
    if (this.sparseOnes) return this.sparseOnes.select1(n);
    if (this.sparseZeros) return this.sparseZeros.select0(n);

    // The n-th one is in the j-th 01-block.
    const j = binarySearchAfterAccess(
      (i) => this.zo.select1(i + 1) - this.z.select1(i + 1),
      n - 1,
      0,
      this.z.length
    );

    // Number of zeros up to and including the jth blocka
    const numCumulativeZeros = this.z.select1(j + 1);

    return numCumulativeZeros + n - 1;
  }

  approxSizeInBits() {
    // ignores fixed-size fields
    if (this.z) return this.z.approxSizeInBits() + this.zo.approxSizeInBits();
    if (this.sparseOnes) return this.sparseOnes.approxSizeInBits();
    return this.sparseZeros.approxSizeInBits();
  }
}

function typedArrayType(maxVal) {
  if (maxVal < 2 ** 32) return Uint32Array;
  return Float64Array;
}