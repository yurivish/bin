import { binarySearchAfter, binarySearchAfterAccess } from './util.js';

export class SparseBitVector {
  constructor(data, length) {
    if (length === undefined) throw new Error("length must be specified explicitly");
    // Each array element contains the index of a one bit, sorted in increasing order.
    this.data = data;
    this.length = length;
    this.numZeros = this.length - this.data.length;
    this.numOnes = this.length - this.numZeros;
  }
  rank1(i) {
    return binarySearchAfter(this.data, i, 0, this.data.length);
  }
  rank0(i) {
    if (i < 0) return 0;
    if (i >= this.length) return this.numZeros;
    return i - this.rank1(i) + 1;
  }
  select1(n) {
    return n < 1 || n > this.numOnes ? -1 : this.data[n - 1];
  }
  select0(n) {
    if (n < 1 || n > this.numZeros) return -1;
    if (this.data.length === 0) return n - 1;
    const i = binarySearchAfterAccess(
      (i) => this.data[i] - i,
      n - 1,
      0,
      this.data.length
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
    return 8 * this.data.length * this.data.BYTES_PER_ELEMENT;
  }
}
