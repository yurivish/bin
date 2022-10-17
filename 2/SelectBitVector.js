import { popcount, trailing0 } from './util';

export class SelectBitVector {
  constructor(length) {
    const n = Math.ceil(length / 32);
    this.blocks = new Uint32Array(n);
    this.precedingZeroBlocks = new Uint32Array(n);
    this.selectSuperblocks = new Uint32Array(n);
    this.numOnes = 0;
    this.numZeroBlocks = 0;
    this.prevUncompressedBlockIndex = -1;
    this.prevOnePosition = -1;
    this.length = length;
  }
  
  // maybe we should do the interleave select superblock creation with block creation, and lift the restriction on strictly ascending ones required 
  // (like the rank bitvector)
  one(position) {
    if (position <= this.prevOnePosition) throw new Error('ones must be added in strictly-ascending order');
    this.prevOnePosition = position;

    let blockIndex;
    let uncompressedBlockIndex = position >>> 5;
    if (uncompressedBlockIndex > this.prevUncompressedBlockIndex) {
      const prevBlockIndex = this.prevUncompressedBlockIndex - this.numZeroBlocks;
      const numNewZeroBlocks = uncompressedBlockIndex - this.prevUncompressedBlockIndex - 1;
      this.numZeroBlocks += numNewZeroBlocks;
      this.prevUncompressedBlockIndex = uncompressedBlockIndex;
      blockIndex = prevBlockIndex + 1;
      this.precedingZeroBlocks[blockIndex] = this.numZeroBlocks;
    } else {
      blockIndex = this.prevUncompressedBlockIndex - this.numZeroBlocks;
    }

    const bitOffset = position & 31;
    if (this.numOnes % 32 === 0) {
      const selectSampleIndex = this.numOnes >>> 5;
      // number of 1 bits in this block before the (32k+1)th bit
      const adjustment = popcount(this.blocks[blockIndex]);
      this.selectSuperblocks[selectSampleIndex] = (blockIndex << 5) | adjustment;
    }
    this.blocks[blockIndex] |= 1 << bitOffset;
    this.prevBlockIndex = blockIndex;
    this.numOnes += 1;
  }

  finish() {
    this.storedLength = this.maxOnePosition + 1;
    const numBlocks = Math.ceil(this.storedLength / 32);    
    if (numBlocks < this.blocks.length) {
      this.blocks = this.blocks.slice(0, numBlocks);
      this.precedingZeroBlocks = this.precedingZeroBlocks.slice(0, numBlocks);
    }

    const numSelectSuperblocks = Math.ceil(this.numOnes / 32);
    if (numSelectSuperblocks < this.selectSuperblocks.length) {
      this.selectSuperblocks = this.selectSuperblocks.slice(0, numSelectSuperblocks);
    }

    // compute rank superblocks (cumulative sum of 1 bits per block)
    this.rankSuperblocks = new Uint32Array(numBlocks);
    this.rankSuperblocks[0] = popcount(this.blocks[0]);
    for (let i = 1; i < numBlocks; i++) {
      this.rankSuperblocks[i] = this.rankSuperblocks[i - 1] + popcount(this.blocks[i]);
    }
  }

  select1(i) {
    if (i < 1) throw new Error('out of bounds: i < 1');
    if (i > this.numOnes) {
      throw new Error('out of bounds: i > numOnes: ' + i + ' vs. ' + this.numOnes);
    }
    const selectSuperblockIndex = (i - 1) >>> 5;
    const selectSuperblock = this.selectSuperblocks[selectSuperblockIndex];
    let blockIndex = selectSuperblock >>> 5;
    let blockRank = this.rankSuperblocks[blockIndex]; // cumulative rank up to and including blockIndex
    const adjustment = selectSuperblock & 31; // number of 1 bits in this block before the (32k+1)th bit
    let prevBlockRank = (selectSuperblockIndex << 5) - adjustment; // cumulative rank of previous block
    while (blockRank < i) {
      prevBlockRank = blockRank;
      blockIndex = blockIndex + 1;
      blockRank = this.rankSuperblocks[blockIndex];
    }
    let block = this.blocks[blockIndex];
    if (block === undefined) throw new Error('undef block');
    for (let r = prevBlockRank + 1; r < i; r++) block &= block - 1;
    return ((this.precedingZeroBlocks[blockIndex] + blockIndex) << 5) + trailing0(block);
  }
}
