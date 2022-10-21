import { popcount } from './util';

// todo: consider the space required during construction; can we reduce it?
// todo: use the same { rank: true, select: true } constructor to support both ops;
// though here, select implies rank. and we may want larger rank blocks for select.
export class RankBitVector {
  constructor(length) {
    // todo: interleave rank and bits blocks for improved rank performance (access slows down)
    // todo: make a RankSelectBitVector for noninterleaved noncompressed rank/select
    const n = Math.ceil(length / 32);
    this.blocks = new Uint32Array(n);
    this.numOnes = 0;
    this.maxOnePosition = -1;
    this.length = length;
  }

  one(position) {
    if (position >= this.length) throw new Error('position must be < length')
    const blockIndex = position >>> 5;
    const bitOffset = position & 31;
    this.blocks[blockIndex] |= 1 << bitOffset;
    this.numOnes += 1;
    if (position > this.maxOnePosition) this.maxOnePosition = position;
  }

  finish() {
    this.storedLength = this.maxOnePosition + 1;
    const numBlocks = Math.ceil(this.storedLength / 32);
    if (numBlocks < this.blocks.length) {
      this.blocks = this.blocks.slice(0, numBlocks);
    }
    // compute rank superblocks (cumulative sum of 1 bits per block)
    this.rankSuperblocks = new Uint32Array(numBlocks);
    this.rankSuperblocks[0] = popcount(this.blocks[0]);
    for (let i = 1; i < numBlocks; i++)
      this.rankSuperblocks[i] = this.rankSuperblocks[i - 1] + popcount(this.blocks[i]);
  }

  rank1(i) {
    if (i < 0) return 0;
    if (i >= this.storedLength) return this.numOnes;
    const blockIndex = i >>> 5;
    const rankSuperbock = this.rankSuperblocks[blockIndex];
    const block = this.blocks[blockIndex];
    const lowBitIndex = i & 31;
    const mask = 0xfffffffe << lowBitIndex;
    return rankSuperbock - popcount(block & mask);
  }

  rank0(i) {
    if (i < 0) return 0;
    // We check against length here, rather than storedLength,
    // since we need to count the implicitly-represented zeros.
    if (i >= this.length) return this.length - this.numOnes;
    // note: the final block is padded with zeros so rank0 will return
    // incorrect results if called with an out-of-bounds index that is
    // within the final block. So we do the bounds checks here too.
    // Can optimize via copy-pasting the rank1 impl in here.
    return i - this.rank1(i) + 1;
  }

  access(i) {
    if (i < 0 || i > this.length) throw new Error('access: out of bounds');
    const blockIndex = i >>> 5;
    const bitOffset = i & 31;
    const block = this.blocks[blockIndex];
    const targetMask = 1 << bitOffset; // mask out the target bit
    return +Boolean(block & targetMask);
  }

  approxSizeInBits() {
    // ignores fixed-size fields
    const blockBits = 8 * this.blocks.length * this.blocks.BYTES_PER_ELEMENT
    const rankSuperblockBits = 8 * this.rankSuperblocks.length * this.rankSuperblocks.BYTES_PER_ELEMENT
    return blockBits + rankSuperblockBits
  }
}
