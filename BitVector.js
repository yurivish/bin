import { popcount } from './util.js';

// todo: consider the space required during construction; can we reduce it?
// todo: use the same { rank: true, select: true } constructor to support both ops;
// though here, select implies rank. and we may want larger rank blocks for select.
// todo: rename position to i/index
export class BitVector {
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
    return i + 1 - this.rank1(i);
  }

  // These select1 and select0 implementations use binary search over the array
  // without a select-based acceleration index, and are thus O(log(length)).
  // Both support hinted binary search: the search range can be specified through
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
    // if (i < 1) throw new Error('out of bounds: i < 1');
    // if (i > this.numOnes) {
    //   throw new Error(`out of bounds: i (${i}) > numOnes (${this.numOnes})`);
    // }
    // Search based on the structure of binarySearchBefore
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
    // Search based on the structure of binarySearchBefore
    while (L < R) {
      const m = (L + R) >>> 1;
      if (this.rank0(m) < i) L = m + 1;
      else R = m;
    }
    return L;
  }

  access(i) {
    if (i < 0 || i > this.length) throw new Error('access: out of bounds at index ' + i +' with length '+this.length);
    if (i >= this.storedLength) return 0;
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
