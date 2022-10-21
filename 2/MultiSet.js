import { RankBitVector } from './RankBitVector';
import { ZeroCompressedBitVector } from './ZeroCompressedBitVector';

// represents a sorted list of integer elements (must be constructed from a sorted list)
// todo: describe how it works (it keeps one bitmap of unique elements, and another of repetitions)
export class RankMultiSet {
  constructor(data) {
    // stores a 1 for every unique value in the universe [0, maximum(data)]
    this.occupancy = new RankBitVector(data.length === 0 ? 0 : data[data.length - 1] + 1);
    // stores a trailing 1 for every unique value in the universe [0, length(data) - 1]
    // e.g. 5555 88 999 becomes 0001 01 001
    this.multiplicity = new ZeroCompressedBitVector(data.length, { select: true });
    this.maxValue = 0;
    if (data.length > 0) {
      let cur = data[0];
      this.occupancy.one(cur);
      for (let index = 1, prev = cur; index < data.length; index++) {
        prev = cur;
        cur = data[index];
        if (prev !== cur) {
          this.occupancy.one(cur);
          this.multiplicity.one(index - 1);
        }
      }
      this.multiplicity.one(data.length - 1);
      this.maxValue = cur;
    }
    this.occupancy.finish();
    this.multiplicity.finish();
    this.length = data.length;
  }
  rank1(i) {
    if (!this.occupancy.rank1) throw new Error('occupancy must support rank1');
    if (!this.multiplicity.select1) throw new Error('multiplicity must support select1');
    if (i > this.maxValue) return this.length;
    const n = this.occupancy.rank1(i);
    if (n === 0) return 0;
    return this.multiplicity.select1(n) + 1;
  }
  approxSizeInBits() {
    return this.occupancy.approxSizeInBits() + this.multiplicity.approxSizeInBits();
  }
}

export class AccessMultiSet {
  constructor(data) {
    // stores a 1 for every unique value in the universe [0, maximum(data)]
    this.occupancy = new ZeroCompressedBitVector(data.length === 0 ? 0 : data[data.length - 1] + 1, { select: true });
    // stores a leading 1 for every unique value in the universe [0, length(data) - 1]
    // e.g. 5555 88 999 becomes 1000 10 100
    this.multiplicity = new RankBitVector(data.length);
    this.maxValue = 0;
    if (data.length > 0) {
      let cur = data[0];
      this.occupancy.one(cur);
      this.multiplicity.one(0);
      for (let index = 1, prev = cur; index < data.length; index++) {
        prev = cur;
        cur = data[index];
        if (prev !== cur) {
          this.occupancy.one(cur);
          this.multiplicity.one(index);
        }
      }
      this.maxValue = cur;
    }
    this.occupancy.finish();
    this.multiplicity.finish();
    this.length = data.length;
  }
  access(i) {
    if (!this.occupancy.select1) throw new Error('occupancy must support select1');
    if (!this.multiplicity.rank1) throw new Error('multiplicity must support rank1');
    if (i < 0) throw new Error('i must be >= 0');
    if (i >= this.length) throw new Error('i must be < length');
    // how many unique elements there are up to element i
    const n = this.multiplicity.rank1(i);
    // tells us what the nth unique element is
    return this.occupancy.select1(n);
  }
  approxSizeInBits() {
    return this.occupancy.approxSizeInBits() + this.multiplicity.approxSizeInBits();
  }
}
