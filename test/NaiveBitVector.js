export class NaiveBitVector {
  constructor(ones, length) {
    // ones: Sorted array containing the index of every 1-bit
    // length: Length of the bit vector in bits
    this.ones = ones;
    this.set = new Set(ones);
    this.sum = new Uint32Array(length);
    this.sum[0] = this.set.has(0);
    for (let i = 1; i < length; i++) {
      this.sum[i] = this.sum[i - 1] + this.set.has(i);
    }
  }
  rank1(i) {
    // Return the number of 1-bits up to and including the bit at index i
    return this.sum[i];
  }
  select1(i) {
    // Return the index of the i-th 1-bit
    return this.ones[i - 1];
  }
  access(i) {
    // Return the value of the bit at index i
    return +this.set.has(i);
  }
}