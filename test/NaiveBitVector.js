export class NaiveBitVector {
  constructor(length) {
    this.length = length;
    this.ones = [];
  }
  one(i) {
    this.ones.push(i);
  }
  finish() {
    this.onesSet = new Set(this.ones);
    // filter duplicates & sort
    this.ones = [];
    this.zeros = [];
    for (let i = 0; i < this.length; i++) {
      if (this.onesSet.has(i)) this.ones.push(i);
      else this.zeros.push(i);
    }
    this.onesSum = new Uint32Array(this.length);
    this.onesSum[0] = this.onesSet.has(0);
    for (let i = 1; i < this.length; i++) {
      this.onesSum[i] = this.onesSum[i - 1] + this.onesSet.has(i);
    }
    this.numOnes = this.onesSum[this.onesSum.length - 1];
  }
  rank1(i) {
    // Return the number of 1-bits up to and including the bit at index i
    if (i < 0) return 0;
    if (i >= this.length) return this.numOnes;
    return this.onesSum[i];
  }
  rank0(i) {
    // Return the number of 0-bits up to and including the bit at index i
    if (i < 0) return 0;
    if (i >= this.length) return this.length - this.numOnes;
    return i + 1 - this.rank1(i);
  }
  select1(i) {
    // Return the index of the i-th 1-bit
    if (i < 1 || i > this.ones.length) return -1;
    return this.ones[i - 1];
  }
  select0(i) {
    // Return the index of the i-th 0-bit
    if (i < 1 || i > this.zeros.length) return -1;
    return this.zeros[i - 1];
  }
  access(i) {
    if (i < 0 || i >= this.length) throw new Error('access: out of bounds')
    // Return the value of the bit at index i
    return +this.onesSet.has(i);
  }
}
