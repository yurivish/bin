import { RankBitVector } from './RankBitVector';
import { SelectBitVector } from './SelectBitVector';

// represents a sorted list of integer elements (must be constructed from a sorted list)
// todo: describe how it works (it keeps one bitmap of unique elements, and another of repetitions)
export class MultiSet {
  // for rank1, set occupancyType = RankBitVector and multiplicityType = SelectBitVector
  // for access, set occupancyType = SelectBitVector and multiplicityType = RankBitVector
  constructor(data, { occupancyType = RankBitVector, multiplicityType = SelectBitVector } = {}) {
    this.occupancy = new occupancyType(data.length === 0 ? 0 : data[data.length - 1] + 1);
    this.multiplicity = new multiplicityType(data.length);
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
    this.length = data.length;
    this.occupancy.finish();
    this.multiplicity.finish();
  }
  rank1(i) {
    if (i > this.maxValue) return this.length;
    const n = this.occupancy.rank1(i);
    if (n === 0) return 0;
    return this.multiplicity.select1(n) + 1;
  }
  access(i) {
    if (i < 0) throw new Error('i must be >= 0')
    if (i >= this.length) throw new Error('i must be < length')
    // how many unique elements there are up to element i
    const n = multiplicity.rank1(i) + 1
    // tells us what the nth unique element is
    return occupancy.select1(j) 
  }
}
