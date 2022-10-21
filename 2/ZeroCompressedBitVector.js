import { RankBitVector } from './RankBitVector';
import { popcount, trailing0 } from './util';

// todo: do the zero compression after all ones exist in a subsequent filter pass
// that's also the time to do select block indexing.
// this would lift the restriction on "all ones must be added in order", which 
// currently prevents the use of this bitvector with the wavelet matrix.
// this would let us add all indexing structures in the same pass, as we walk
// along and compress.
// todo: is there a way to factor out the common functionality to
// avoid the conditional execution for rank or select indices?
// todo: consider the space required during construction; can we reduce it?
// a: yes, we can avoid storing `precedingZeroBlocks` if `isZeroBlock` has select1
// (will be slower; but will replace `this.precedingZeroBlocks[blockIndex]` with
// `this.isZeroBlock.select1(blockIndex)`)
//
// The rank index `isZeroBlock` stores an additional bitvector of `numBlocks` bits,
// for one additional bit per block in the original-length bitvector.
// The rank index `rankSuperblocks` stores an additional block per nonzero block.
// The select index `precedingZeroBlocks` stores an additional block per nonzero block.
// The select index `selectSuperblocks` stores an additional bit for every original 1 bit.
// Selecting uses rank superblocks, so we build those in both cases.
// This bit vector is configurable based on whether you want to support one  rank/select.
export class ZeroCompressedBitVector {
  constructor(length, { rank = false, select = false } = {}) {
    const initialNumBlocks = Math.ceil(length / 32);
    // stores nonzero blocks
    this.blocks = new Uint32Array(initialNumBlocks);
    if (select) {
      // for each nonzero block, how many zero blocks came before it
      this.precedingZeroBlocks = new Uint32Array(initialNumBlocks);
      // superblocks indexing every 32 1 bits; for details, see
      // the paper "Fast, Small, Simple Rank/Select on Bitmaps".
      this.selectSuperblocks = new Uint32Array(initialNumBlocks);
    }
    if (rank) {
      // for each original block, whether it is a zero block
      this.isZeroBlock = new RankBitVector(initialNumBlocks);
    }
    this.numOnes = 0;
    this.numZeroBlocks = 0;
    this.prevUncompressedBlockIndex = -1;
    this.prevOnePosition = -1;
    this.length = length;
    this.rank = rank;
    this.select = select;
  }

  // maybe we should do the interleave select superblock creation with block creation, and lift the restriction on strictly ascending ones required
  // (like the rank bitvector)
  one(position) {
    if (position < 0) throw new Error('position must be >= 0');
    if (position >= this.length) throw new Error('position must be < length');
    if (position <= this.prevOnePosition) throw new Error('ones must be added in strictly-ascending order');
    this.prevOnePosition = position;

    let blockIndex;
    let uncompressedBlockIndex = position >>> 5;
    if (uncompressedBlockIndex > this.prevUncompressedBlockIndex) {
      const prevBlockIndex = this.prevUncompressedBlockIndex - this.numZeroBlocks;
      const numNewZeroBlocks = uncompressedBlockIndex - this.prevUncompressedBlockIndex - 1;
      if (this.rank) {
          for (let i = 0; i < numNewZeroBlocks; i++) {
          this.isZeroBlock.one(this.prevUncompressedBlockIndex + i + 1);
        }
      }
      this.numZeroBlocks += numNewZeroBlocks;
      this.prevUncompressedBlockIndex = uncompressedBlockIndex;
      blockIndex = prevBlockIndex + 1;
      if (this.select) this.precedingZeroBlocks[blockIndex] = this.numZeroBlocks;
    } else {
      blockIndex = this.prevUncompressedBlockIndex - this.numZeroBlocks;
    }

    const bitOffset = position & 31;
    if (this.select) {
      if (this.numOnes % 32 === 0) {
        const selectSampleIndex = this.numOnes >>> 5;
        // number of 1 bits in this block before the (32k+1)th bit
        const adjustment = popcount(this.blocks[blockIndex]);
        this.selectSuperblocks[selectSampleIndex] = (blockIndex << 5) | adjustment;
      }
    }
    this.blocks[blockIndex] |= 1 << bitOffset;
    this.prevBlockIndex = blockIndex;
    this.numOnes += 1;
  }

  finish() {
    if (this.rank) this.isZeroBlock.finish();
    this.storedLength = this.prevOnePosition + 1;
    const numBlocks = Math.ceil(this.storedLength / 32) - this.numZeroBlocks;
    if (numBlocks < this.blocks.length) {
      this.blocks = this.blocks.slice(0, numBlocks);
      if (this.select) this.precedingZeroBlocks = this.precedingZeroBlocks.slice(0, numBlocks);
    }

    if (this.select) {
      const numSelectSuperblocks = Math.ceil(this.numOnes / 32);
      if (numSelectSuperblocks < this.selectSuperblocks.length) {
        this.selectSuperblocks = this.selectSuperblocks.slice(0, numSelectSuperblocks);
      }
    }

    if (this.rank || this.select) {
      // compute rank superblocks (cumulative sum of 1 bits per block)
      this.rankSuperblocks = new Uint32Array(numBlocks);
      this.rankSuperblocks[0] = popcount(this.blocks[0]);
      for (let i = 1; i < numBlocks; i++) {
        this.rankSuperblocks[i] = this.rankSuperblocks[i - 1] + popcount(this.blocks[i]);
      }
    }
  }

  rank1(i) {
    if (!this.rank) throw new Error('no indexing structures for rank (see constructor)');
    if (i < 0) return 0;
    const uncompressedBlockIndex = i >>> 5;
    const numPrecedingZeroBlocks = this.isZeroBlock.rank1(uncompressedBlockIndex - 1);
    // todo: find a way to avoid this extra memory access in most cases
    // (it's often in the same block as the call to rank1 above)
    const isZeroBlock = this.isZeroBlock.access(uncompressedBlockIndex); // returns 0 or 1
    i -= numPrecedingZeroBlocks << 5;
    if (i >= this.storedLength) return this.numOnes;
    const blockIndex = i >>> 5;
    // whether i points inside a zero block; if yes, we will need 
    // to adjust lowBitIndex so that we point to the last bit in 
    // the closest nonzero block that precedes it.
    const lowBitIndex = isZeroBlock ? 0 : i & 31; // could be (1 - isZeroBlock) * (i & 31)
    const rankSuperbock = this.rankSuperblocks[blockIndex];
    const block = this.blocks[blockIndex];
    const mask = 0xfffffffe << lowBitIndex;
    return rankSuperbock - popcount(block & mask);
  }

  access(i) {
    const uncompressedBlockIndex = i >>> 5;
    const numPrecedingZeroBlocks = this.isZeroBlock.rank1(uncompressedBlockIndex - 1);
    const isZeroBlock = this.isZeroBlock.access(uncompressedBlockIndex);
    i -= numPrecedingZeroBlocks << 5;
    if (i < 0 || i > this.length) throw new Error('access: out of bounds');
    const blockIndex = i >>> 5;
    const bitOffset = isZeroBlock ? 0 : i & 31;
    const block = this.blocks[blockIndex];
    const targetMask = 1 << bitOffset; // mask out the target bit
    return Boolean(block & targetMask);
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

  select1(i) {
    if (!this.select) throw new Error('no indexing structures for select (see constructor)');
    if (i < 1) throw new Error('out of bounds: i < 1');
    if (i > this.numOnes) {
      throw new Error(`out of bounds: i (${i}) > numOnes (${this.numOnes})`);
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
    const numPrecedingZeroBlocks = this.precedingZeroBlocks[blockIndex];
    return ((numPrecedingZeroBlocks + blockIndex) << 5) + trailing0(block);
  }

  approxSizeInBits() {
    // ignores fixed-size fields
    const blockBits = 8 * this.blocks.length * this.blocks.BYTES_PER_ELEMENT;
    let numBits = blockBits;
    if (this.rankSuperblockBits) {
      numBits += 8 * this.rankSuperblocks.length * this.rankSuperblocks.BYTES_PER_ELEMENT;
    }
    if (this.selectSuperblockBits) {
      numBits += 8 * this.selectSuperblocks.length * this.selectSuperblocks.BYTES_PER_ELEMENT;
    }
    if (this.precedingZeroBlocks) {
      numBits += 8 * this.precedingZeroBlocks.length * this.precedingZeroBlocks.BYTES_PER_ELEMENT;
    }
    if (this.isZeroBlock) {
      numBits += this.isZeroBlock.approxSizeInBits();
    }
    return numBits;
  }
}
