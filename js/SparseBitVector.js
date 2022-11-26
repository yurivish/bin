import { binarySearchAfter, binarySearchAfterAccess } from './util.js';

// Stores a collection of one bits by explicitly representing them in an array.
export class SparseOneBitVector {
  constructor(length, ones) {
    if (length === undefined) throw new Error('length must be specified explicitly');
    this.length = length
    this.ones = ones ?? [];
  }

  one(i) {
    this.ones.push(i);
  }

  finish() {
    const ones = this.ones;
    const len = this.ones.length;

    if (len > 0) {
      // check for errors
      const max = ones[len - 1] 
      if (max > this.length) {
        throw new Error('each one index must be < length');
      }
      let prev = this.ones[0];
      for (let i = 0; i < len; i++) {
        const cur = ones[i];
        if (prev > cur) {
          throw new Error('ones must be supplied in strictly ascending order');
        }
        prev = cur;
      }
      // turn the ones array into a typed array if needed
      // duck type check for whether ones is already a typed array
      if (ones.subarray === undefined) {
        const T = max < 2 ** 32 ? Uint32Array : Float64Array;
        this.ones = new T(ones);
      }
    }

    this.numZeros = this.length - this.ones.length;
    this.numOnes = this.length - this.numZeros;
    return this;
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
    const i = binarySearchAfterAccess((i) => this.ones[i] - i, n - 1, 0, this.ones.length);
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

// TODO: Can we remove the indirection associated with inheritance?
// Maybe make a standard old-style constructor function with a `flip` argument.
export class SparseZeroBitVector extends SparseOneBitVector {
  zero(i) {
    return super.one(i);
  }

  rank1(i) {
    return super.rank0(i);
  }

  rank0(i) {
    return super.rank1(i);
  }

  select1(n) {
    return super.select0(n);
  }

  select0(n) {
    return super.select1(n);
  }

  access(i) {
    // Quick hack. Can do better.
    return super.rank0(i) - super.rank0(i - 1);
  }
}
