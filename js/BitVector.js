import { popcount } from './util.js';

// todo: consider the space required during construction; can we reduce it?
// todo: use the same { rank: true, select: true } constructor to support both ops;
// though here, select implies rank. and we may want larger rank blocks for select.
export class BitVector {
  // todo: allow initializing from blocks/rankSuperblocks; can use select
  // to determine the true maxOneIndex (pre-init it to length since select uses it)
  constructor(length) {
    const n = Math.ceil(length / 32);
    this.blocks = new Uint32Array(n);
    this.numOnes = 0;
    this.maxOneIndex = -1;
    this.length = length;
  }

  one(i) {
    if (i < 0 || i >= this.length) throw new Error('i must be in [0, length)');
    const blockIndex = i >>> 5;
    const bitOffset = i & 31;
    this.blocks[blockIndex] |= 1 << bitOffset;
    this.numOnes += 1;
    if (i > this.maxOneIndex) this.maxOneIndex = i;
  }

  finish() {
    const numBlocks = Math.ceil((this.maxOneIndex + 1) / 32);
    // If the number of used blocks is less than x% of the
    // allocated blocks, resize downwards
    const resizeProportion = 0.9;
    if (numBlocks < resizeProportion * this.blocks.length) {
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
    if (i > this.maxOneIndex) return this.numOnes;
    const blockIndex = i >>> 5;
    const rankSuperblock = this.rankSuperblocks[blockIndex];
    const block = this.blocks[blockIndex];
    const lowBitIndex = i & 31;
    const mask = 0xfffffffe << lowBitIndex;
    return rankSuperblock - popcount(block & mask);
  }

  rank0(i) {
    if (i < 0) return 0;
    if (i >= this.length) return this.length - this.numOnes;
    return i - this.rank1(i) + 1;
  }

  // These select1 and select0 implementations use binary search over the array
  // without a select-based acceleration index, and are thus O(log(length)).
  // Both support hinted search: the search range can be specified through
  // input arguments for those cases where the sought-after bit is known to be
  // confined to a particular index range.
  // Sampled select blocks could similarly cut down the search range but are not
  // implemented here (and would cost additional space, which would be configurable
  // through sample rate tuning).
  // Other optimization opportunities: Binary search over the 32 blocks when performing
  // zero-compressed select. Perform exponential rather than binary search (determine
  // when and why this might be more efficient).
  select1(i, L = 0, R = this.length) {
    if (i < 1 || i > this.numOnes) return -1;
    while (L < R) {
      const m = (L + R) >>> 1;
      if (this.rank1(m) < i) L = m + 1;
      else R = m;
    }
    return L;
  }

  select0(i, L = 0, R = this.length) {
    const numZeros = this.length - this.numOnes;
    if (i < 1 || i > numZeros) return -1;
    while (L < R) {
      const m = (L + R) >>> 1;
      if (this.rank0(m) < i) L = m + 1;
      else R = m;
    }
    return L;
  }

  access(i) {
    if (i < 0 || i >= this.length)
      throw new Error('access: out of bounds');
    if (i > this.maxOneIndex) return 0;
    const blockIndex = i >>> 5;
    const bitOffset = i & 31;
    const block = this.blocks[blockIndex];
    const targetMask = 1 << bitOffset; // mask out the target bit
    return +Boolean(block & targetMask);
  }

  approxSizeInBits() {
    // ignores fixed-size fields
    const blockBits = 8 * this.blocks.length * this.blocks.BYTES_PER_ELEMENT;
    const rankSuperblockBits = 8 * this.rankSuperblocks.length * this.rankSuperblocks.BYTES_PER_ELEMENT;
    return blockBits + rankSuperblockBits;
  }
}
