import { popcount } from './util';

export class RankBitVector {
  constructor(length) {
    const n = Math.ceil(length / 32);
    this.blocks = new Uint32Array(n);
    this.numOnes = 0;
    this.prevOnePosition = -1;
  }

  one(position) {
    if (position <= this.prevOnePosition) throw new Error('ones must be added in strictly-ascending order');
    this.prevOnePosition = position;
    const blockIndex = position >>> 5;
    const bitOffset = position & 31;
    this.blocks[blockIndex] |= 1 << bitOffset;
    this.numOnes += 1;
    this.prevOnePosition = position;
  }

  finish() {
    this.length = this.prevOnePosition + 1;
    const numBlocks = Math.ceil(this.length / 32);
    if (numBlocks < this.blocks.length) {
      this.blocks = this.blocks.slice(0, numBlocks);
    }
    console.log('numBlocks',numBlocks,this.blocks.length)
    // compute rank superblocks (cumulative sum of 1 bits per block)
    this.rankSuperblocks = new Uint32Array(numBlocks);
    this.rankSuperblocks[0] = popcount(this.blocks[0]);
    for (let i = 1; i < numBlocks; i++)
      this.rankSuperblocks[i] = this.rankSuperblocks[i - 1] + popcount(this.blocks[i]);
  }

  rank1(i) {
    if (i < 0) return 0;
    if (i >= this.length) return this.numOnes;
    const blockIndex = i >>> 5;
    const rankSuperbock = this.rankSuperblocks[blockIndex];
    const block = this.blocks[blockIndex];
    const lowBitIndex = i & 31;
    const mask = 0xfffffffe << lowBitIndex;
    return rankSuperbock - popcount(block & mask);
  }
}
