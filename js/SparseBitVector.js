import { binarySearchAfter, binarySearchAfterAccess } from './util.js';

// Stores a collection of one bits by explicitly representing them in an array.
export class SparseBitVector {
  constructor(ones, length) {
    if (length === undefined)
      throw new Error("length must be specified explicitly");
    // todo: error if ones is not sorted

    if (length < ones[ones.length - 1])
      throw new Error(
        "length must be >= the largest one index: " +
          length +
          "," +
          ones[ones.length - 1]
      );

    let prev = ones[0];
    for (let i = 1; i < ones.length; i++) {
      const cur = ones[i];
      if (prev >= cur) {
        throw new Error("ones must be supplied in strictly ascending order");
      }
    }

    // Each array element contains the index of a one bit, sorted in increasing order.
    this.ones = ones;
    this.length = length;
    this.numZeros = this.length - this.ones.length;
    this.numOnes = this.length - this.numZeros;
  }
  rank1(i) {
    return binarySearchAfter(this.ones, i, 0, this.ones.length);
  }
  rank0(i) {
    if (i < 0) return 0;
    if (i >= this.length) return this.numZeros;
    return i - this.rank1(i) + 1;
  }
  select1(n) {
    return n < 1 || n > this.numOnes ? -1 : this.ones[n - 1];
  }
  select0(n) {
    if (n < 1 || n > this.numZeros) return -1;
    if (this.ones.length === 0) return n - 1;
    const i = binarySearchAfterAccess(
      (i) => this.ones[i] - i,
      n - 1,
      0,
      this.ones.length
    );
    // i is the index of the block containing the i-th zero (it's in the run preceding the 1).
    // There are i ones before it.
    return n + i - 1;
  }
  access(i) {
    // Quick hack. Can do better.
    return this.rank1(i) - this.rank1(i - 1);
  }
  approxSizeInBits() {
    return 8 * this.ones.length * this.ones.BYTES_PER_ELEMENT;
  }
}