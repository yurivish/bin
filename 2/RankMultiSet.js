import { RankBitVector } from './RankBitVector';
import { SelectBitVector } from './SelectBitVector';

export class RankMultiSet {
  constructor(data) {
    this.occupancy = new RankBitVector(data.length === 0 ? 0 : data[data.length - 1] + 1);
    this.multiplicity = new SelectBitVector(data.length);
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
}
