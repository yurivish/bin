class NaiveBitVector {
  constructor(ones, length) {
    // ones: Sorted array containing the index of every 1-bit
    // length: Length of the bit vector in bits
    this.ones = ones;
    this.set = new Set(ones);
    this.sum = new Uint32Array(length);
    this.sum[0] = this.set.has(0);
    for (let i = 1; i < length; i++) {
      this.sum[i] = this.sum[i - 1] + this.set.has(i);
    }
  }
  rank1(i) {
    // Return the number of 1-bits up  to and including the bit at index i
    return this.sum[i];
  }
  select1(i) {
    // Return the index of the i-th 1-bit
    return this.ones[i - 1];
  }
  access(i) {
    // Return the value of the bit at index i
    return +this.set.has(i);
  }
}

// fast rank
class SimpleBitVector {
  constructor(ones, length) {
    // `ones` should be a sorted array of the unique locations of 1 bits
    const n = Math.ceil(length / 32 /* does this work for large numbers? */); // can also use u64
    const blocks = new Uint32Array(n); // packed bit vector blocks
    for (let i = 0; i < ones.length; i++) {
      const pos = ones[i]; // bit position to set
      const blockIndex = pos >>> 5; // = floor(value / 32)
      const bitOffset = pos & 31; // bit offset within the block
      blocks[blockIndex] |= 1 << bitOffset;
    }
    const rankSuperblocks = new Uint32Array(n); // cumulative sum of 1 bits per block
    rankSuperblocks[0] = popcount(blocks[0]);
    for (let i = 1; i < n; i++)
      rankSuperblocks[i] = rankSuperblocks[i - 1] + popcount(blocks[i]);
    this.blocks = blocks;
    this.rankSuperblocks = rankSuperblocks;
    this.numOnes = ones.length;
    this.length = length;
  }
  access(i) {
    const blockIndex = i >>> 5;
    const bitOffset = i & 31;
    const block = this.blocks[blockIndex];
    const target = 1 << bitOffset; // mask out the target bit
    return (block & target) === target; // return true if the masked bit is set
  }
  rank1(i) {
    if (i < 0) return 0;
    // if (i < 0) throw new Error("i < 0");
    if (i >= this.length) return this.numOnes;

    const blockIndex = i >>> 5;
    const rankSuperbock = this.rankSuperblocks[blockIndex];
    const block = this.blocks[blockIndex];
    const lowBitIndex = i & 31;
    const mask = 0xfffffffe << lowBitIndex;
    const v = rankSuperbock - popcount(block & mask);
    // todo: check bounds
    // if (Number.isNaN(v)) `throw new Error("nan result from simplevector.rank1");
    return v;
  }
  rank0(i) {
    // note: the final block is padded with zeros so rank0 will return
    // incorrect results if called with an out-of-bounds index that is
    // within the final block.
    return i - this.rank1(i) + 1;
  }
  // rank1_early_out(i) {
  //   const blockIndex = i >>> 5;
  //   const rankSuperbock = this.rankSuperblocks[blockIndex];
  //   const block = this.blocks[blockIndex]; // declaring here improves random access perf
  //   const lowBitIndex = i & 31;
  //   if (lowBitIndex === 31) return rankSuperbock;
  //   const mask = 0xfffffffe << lowBitIndex;
  //   return rankSuperbock - popcount(block & mask);
  // }
}

// fast rank and select (in progress)
class SimpleSelectBitVector {
  select1(i) {
    // arg is not an index; rename from i?
    // if (i < 1) return 0; // NOTE: unconventional!
    if (i < 1) throw new Error("out of bounds: i < 1");
    if (i > this.numOnes)
      throw new Error(
        "out of bounds: i > numOnes: " + i + " vs. " + this.numOnes
      );
    const selectSuperblockIndex = (i - 1) >>> 5; // >>> SsBitsPow2
    const selectSuperblock = this.selectSuperblocks[selectSuperblockIndex];
    const adjustment = selectSuperblock & 31; // number of 1 bits in this block before the (32k+1)th bit
    let prevBlockRank = (selectSuperblockIndex << 5) - adjustment; // cumulative rank of previous block
    let blockIndex = selectSuperblock >>> 5; // >>> SsBitsPow2; index of next block to scan
    let blockRank = this.rankSuperblocks[blockIndex]; // cumulative rank up to and including blockIndex
    while (blockRank < i) {
      prevBlockRank = blockRank;
      blockIndex = blockIndex + 1;
      blockRank = this.rankSuperblocks[blockIndex];
    }
    let block = this.blocks[blockIndex];
    if (block === undefined) throw new Error("undef block");
    for (let r = prevBlockRank + 1; r < i; r++) block &= block - 1;
    return (blockIndex << 5) /* << bitsPerBlockPow2 */ + trailing0(block);
    // hint can be {blockIndex, prevBlockRank} encoded like a select block, or similar; to allow us to avoid the select superblock lookup and having to iterate forwards from it...
  }

  constructor(ones, length) {
    // console.clear();
    // `ones` should be a sorted array of the locations of 1 bits
    const SrBits = 32; // rank sample rate (also the size of the basic block)
    const SsBits = 32; // select sample rate (number of 1 bits per select superblock)
    const selectSuperblocks = new Uint32Array(Math.ceil(ones.length / SsBits));
    let selectSampleIndex = 0;
    const n = Math.ceil(length / SrBits); // can also use u64
    const blocks = new Uint32Array(n); // packed bit vector blocks
    for (let i = 0; i < ones.length; i++) {
      const pos = ones[i]; // bit position to set
      const blockIndex = pos >>> 5; // = floor(value / 32) (32 === bitsPerBlock)
      const bitOffset = pos & 31; // bit offset within the block (31 == bitsPerBlock - 1)
      if (i % SsBits === 0) {
        // number of 1 bits in this block before the (32k+1)th bit
        const adjustment = popcount(blocks[blockIndex]);
        selectSuperblocks[selectSampleIndex] = (blockIndex << 5) | adjustment;
        selectSampleIndex++;
      }
      blocks[blockIndex] |= 1 << bitOffset;
    }
    const rankSuperblocks = new Uint32Array(n); // cumulative sum of 1 bits per block
    let cumulativeRank = 0;
    for (let i = 0; i < n; i++) {
      cumulativeRank += popcount(blocks[i]);
      rankSuperblocks[i] = cumulativeRank;
    }
    this.blocks = blocks;
    this.rankSuperblocks = rankSuperblocks;
    this.selectSuperblocks = selectSuperblocks;
    this.numOnes = ones.length;
  }

  access(i) {
    const blockIndex = i >>> 5;
    const bitOffset = i & 31;
    const block = this.blocks[blockIndex];
    const target = 1 << bitOffset; // mask out the target bit
    return (block & target) === target; // return true if the masked bit is set
  }
  rank1_early_out(i) {
    const blockIndex = i >>> 5;
    const rankSuperbock = this.rankSuperblocks[blockIndex];
    const block = this.blocks[blockIndex]; // declaring here improves random access perf
    const lowBitIndex = i & 31;
    if (lowBitIndex === 31) return rankSuperbock;
    const mask = 0xfffffffe << lowBitIndex;
    return rankSuperbock - popcount(block & mask);
  }
  rank1(i) {
    const blockIndex = i >>> 5;
    const rankSuperbock = this.rankSuperblocks[blockIndex];
    const block = this.blocks[blockIndex];
    const lowBitIndex = i & 31;
    const mask = 0xfffffffe << lowBitIndex;
    return rankSuperbock - popcount(block & mask);
  }
  rank0(i) {
    // note: the final block is padded with zeros so rank0 will return
    // incorrect results if called with an out-of-bounds index that is
    // within the final block.
    return i - this.rank1(i) + 1;
  }
}

// returns block indices of all zero 32-blocks
const zeroBlocks = (a, universeSize) => {
  const blocks = new Uint8Array(Math.ceil(universeSize / 32));
  for (let i = 0; i < a.length; i++) blocks[a[i] >>> 5] = 1;
  // return d3.sum(blocks);
  const ones = new Uint32Array(d3.sum(blocks, (d) => !d));
  let n = 0;
  for (let i = 0; i < blocks.length; i++) if (blocks[i] === 0) ones[n++] = i;
  return ones;
}

const uniq = (a) => {
  let unique = new Uint32Array(a); // filter in place
  let prev = unique[0];
  let n = 1;
  for (let i = 1; i < unique.length; i++) {
    const cur = unique[i];
    if (prev !== cur) {
      unique[n++] = cur;
      prev = cur;
    }
  }
  return new Uint32Array(unique.slice(0, n));
}

class PlainMultiSet {
  constructor(data, universeSize) {
    // data is sorted, may contain repetitions; todo: special-case for speed when no repetitions?
    let unique = new Uint32Array(data); // filter in place
    let prev = unique[0];
    let n = 1;
    for (let i = 1; i < unique.length; i++) {
      const cur = unique[i];
      if (prev !== cur) {
        unique[n++] = cur;
        prev = cur;
      }
    }
    unique = new Uint32Array(unique.subarray(0, n));
    this.occupancy = new SimpleBitVector(unique, universeSize);
    this.multiplicity = new SimpleSelectBitVector(
      multiplicityOnes(data),
      data.length
    );
    this.max = data[data.length - 1]; // todo: does not handle zero data;
    this.length = data.length;
    this.universeSize = universeSize;
  }
  rank1(i) {
    if (i < 0) return 0; //throw new Error("i < 0");
    if (i > this.max) return this.length;
    const offset = this.occupancy.rank1(i);
    if (offset === 0) return 0;
    return this.multiplicity.select1(offset); // searchRight
  }
}

class MultiSet {
  constructor(data, universeSize) {
    // data is sorted, may contain repetitions; todo: special-case for speed when no repetitions?
    let unique = uniq(data);
    const zeroPositions = zeroBlocks(unique, Math.ceil(unique.length / 32));
    this.zeros = new SimpleBitVector(
      zeroPositions,
      Math.ceil(unique.length / 32)
    );
    let nz = 0;
    for (let i = 0; i < unique.length; i++) {
      const d = unique[i];
      const block = d >>> 5;
      // how many zero blocks were there before this block?
      while (zeroPositions[nz] < block && nz < zeroPositions.length) nz++;
      unique[i] -= nz << 5;
    }
    // return unique;

    this.occupancy = new SimpleBitVector(unique, universeSize);
    this.multiplicity = new SimpleSelectBitVector(
      multiplicityOnes(data),
      data.length
    );
    this.max = data[data.length - 1]; // todo: does not handle zero data;
    this.length = data.length;
    this.universeSize = universeSize;
  }
  rank1(i) {
    // figure out how to support the optimization where we don't rerun the select if the answer is the same as the previous rank answer
    // if (i < 0) return 0; // throw new Error("i < 0");
    if (i > this.max) return this.length;
    i -= this.zeros.rank1(i >>> 5) << 5;
    // if (i < 0) return 0;
    const offset = this.occupancy.rank1(i);
    if (offset === 0) return 0;
    
    return this.multiplicity.select1(offset); // searchRight
  }
}

// fast rank and select w/ configurable sampling rates (slowest)
class BitVector {
  select1(i) {
    if (i < 1) throw new Error("out of bounds: i < 1");
    if (i > this.numOneBits) throw new Error("out of bounds: i > numOneBits");
    const selectSuperblockIndex = (i - 1) >>> this.SsBitsPow2;
    const selectSuperblock = this.sampleSelect[selectSuperblockIndex];
    // blockIndex points to the next block to scan;
    let blockIndex = selectSuperblock >>> this.bitsPerBlockPow2;
    const correction = selectSuperblock & this.bitsPerBlockMask;
    // rank represented by the superblock
    let superblockRank = selectSuperblockIndex * this.SsBits;
    // rank up to the preceding block boundary, ie. up to but not including blocks[blockIndex]
    let r = superblockRank - correction;
    const blocks = this.blocks;
    if (r < i) {
      // Accelerate searches by stepping through rank superblocks
      const sampleRank = this.sampleRank;
      let rankSuperblockIndex = blockIndex >>> this.SrPow2;
      while (rankSuperblockIndex < this.sampleRank.length) {
        const nextRank = this.sampleRank[rankSuperblockIndex]; // note: may be less than r!
        const nextBlockIndex = (rankSuperblockIndex << this.SrPow2) + 1;
        if (nextRank >= i) break;
        if (nextBlockIndex > blockIndex) {
          // note: can only be lessthan the first time
          r = nextRank;
          blockIndex = nextBlockIndex;
        }
        rankSuperblockIndex += 1;
      }
    }

    // todo: is there a more efficient way to phrase this (and similar) loops?
    // the trouble is that we can't index into blocks until we do the index check,
    // so the while loop has two place where it can conditionally break out.
    while (blockIndex < blocks.length) {
      const next = r + popcount(blocks[blockIndex]);
      if (next >= i) break;
      r = next;
      blockIndex += 1;
    }

    let block = blocks[blockIndex];
    let target = (i - r) >>> 0;
    if (popcount(block) < target) throw new Error("wat"); // implies implementation error
    for (r++; r < i; r++) block &= block - 1;
    return (blockIndex << this.bitsPerBlockPow2) + trailing0(block);
  }

  rank1(i) {
    if (i < 0) throw new Error("out of bounds: i < 0");
    if (i >= this.length) throw new Error("out of bounds: i > length");
    const blocks = this.blocks;
    const lastBlockIndex = i >>> this.bitsPerBlockPow2;
    const rankSuperblockIndex = lastBlockIndex >>> this.SrPow2;
    let r = this.sampleRank[rankSuperblockIndex];
    // our rank blocks cover the last block if they can; thus,
    // blockIndex is the index of the next block that we haven't scanned yet
    let blockIndex = (rankSuperblockIndex << this.SrPow2) + 1;
    const blockBitOffset = i & this.bitsPerBlockMask;
    const target = 1 << blockBitOffset;
    const loBits = target | (target - 1);
    const hiBits = ~loBits; // todo: can we initialize r to lastBlock - extraBlockBits?
    const lastBlock = blocks[lastBlockIndex];
    const extraBlockBits = popcount(lastBlock & hiBits);

    if (blockIndex < lastBlockIndex) {
      const sampleSelect = this.sampleSelect;
      // index of the next superblock to scan; pick the one that represents rank r,
      // since that select block will point to the block with a 1 greater than r.
      // note that this index may be beyond sampleSelect.length, since we store the
      // jth block only when there are Ss*j+1 ones.
      // For example, if there are zero ones, we store no select blocks at all, but
      // selectSuperblockIndex will be 0. Similarly, if there are exactly Ss
      // ones, selectSuperblockIndex will be 1 but there will be no block at that index
      // since there is no block with the next 1 for it to point to.
      let selectSuperblockIndex = r >>> this.SsBitsPow2;
      while (selectSuperblockIndex < sampleSelect.length) {
        const next = sampleSelect[selectSuperblockIndex];
        const nextBlockIndex = next >>> this.bitsPerBlockPow2; // note: may be less than blockIndex!
        if (nextBlockIndex > lastBlockIndex) break;
        if (nextBlockIndex > blockIndex) {
          blockIndex = nextBlockIndex;
          const correction = next & this.bitsPerBlockMask;
          r = (selectSuperblockIndex << this.SsBitsPow2) - correction;
        }
        selectSuperblockIndex += 1;
      }
    }

    for (let i = blockIndex; i <= lastBlockIndex; i++) {
      r += popcount(blocks[i]);
    }

    return r - extraBlockBits;
  }

  constructor(
    ones,
    length,
    {
      blockType = Uint32Array,
      SrPow2 = 0, // Rank superblock sampling rate, specified as a power of two
      SsPow2 = 0 // Select superblock sampling rate, specified as a power of two
    } = {}
  ) {
    // note: ones must be sorted ascending and unique; otherwise it will currently silently break an invariant.
    // can run through to check

    this.blockType = blockType;
    this.SrPow2 = SrPow2;
    this.SsPow2 = SsPow2;
    if (length === undefined) throw new Error("missing length");
    if (ones.length > 0 && length < ones[ones.length - 1])
      throw new Error("length must be > the highest one bit");
    this.bitsPerBlockPow2 = Math.log2(8 * this.blockType.BYTES_PER_ELEMENT);
    this.bitsPerBlock = 2 ** this.bitsPerBlockPow2;
    this.bitsPerBlockMask = this.bitsPerBlock - 1;
    this.Sr = 2 ** this.SrPow2;
    this.Ss = 2 ** this.SsPow2;
    this.SrBits = this.Sr * this.bitsPerBlock;
    this.SsBits = this.Ss * this.bitsPerBlock;
    const blocks = new this.blockType(Math.ceil(length / this.bitsPerBlock));
    const sampleRank = new Uint32Array(Math.ceil(length / this.SrBits));
    // Each select block represents k * `SsBits` 1 bits, and
    // points to the block containing the next 1 bit after that block.
    // Each select block also stores correction information that allows
    // one to recover the cumulative rank up to the preceding block.
    // We add a select block only when there is a 1 bit after it, so if
    // numOneBits % SsBits == 0, they will not be represented
    // by a select block.
    const sampleSelect = new Uint32Array(Math.ceil(ones.length / this.SsBits));
    let selectSampleIndex = 0;
    let prevBitIndex = -1;
    for (let i = 0; i < ones.length; i++) {
      const bitIndex = ones[i];
      if (bitIndex === prevBitIndex)
        throw new Error("duplicate 1-bit positions not allowed"); // just skipping would throw off the math
      prevBitIndex = bitIndex;
      const blockIndex = bitIndex >>> this.bitsPerBlockPow2; // block containing the i-th bit
      const blockBitOffset = bitIndex & this.bitsPerBlockMask; // bitIndex % this.bitsPerBlock; // bit offset within the block
      if (i % this.SsBits === 0) {
        const blockRankBefore = popcount(blocks[blockIndex]);
        sampleSelect[selectSampleIndex] =
          (blockIndex << this.bitsPerBlockPow2) | blockRankBefore;
        selectSampleIndex += 1;
      }
      blocks[blockIndex] |= 1 << blockBitOffset; // set the bit
    }

    // sampleRank[j] = rank1(B, j · Sr)
    // - sampleRank[0] = up to block 0, inclusive
    // - sampleRank[1] = up to block Sr, inclusive
    // - ...
    // - as a consequence, will store zero blocks if the length is zero
    // sampleSelect[j] = select1(B, j · Ss + 1)
    // - sampleSelect[0] = represents 0 1-bits, records the block with the 1st one-bit
    // - sampleSelect[1] = represents Ss 1-bits, records the block with the (Ss+1)th one-bit
    // - ...
    // - added lazily, when we see the (Ss+1)th one-bit, so # blocks = ceil((#ones)/Ss)
    let cumulativeRank = 0;
    let blockIndex = 0;
    for (let j = 0; j < sampleRank.length; j++) {
      const lastBlockIndex = j * this.Sr;
      while (blockIndex <= lastBlockIndex) {
        cumulativeRank += popcount(blocks[blockIndex]);
        blockIndex += 1;
      }
      sampleRank[j] = cumulativeRank;
    }

    this.blocks = blocks;
    this.sampleRank = sampleRank;
    this.sampleSelect = sampleSelect;
    this.SrBits = this.Sr * this.bitsPerBlock;
    this.SsBits = this.Ss * this.bitsPerBlock;
    this.SrBitsPow2 = this.SrPow2 + this.bitsPerBlockPow2;
    this.SsBitsPow2 = this.SsPow2 + this.bitsPerBlockPow2;
    this.length = length;
  }
}

// popcount = bit.popcount
// Based on an implementation by @ashaffer: https://github.com/micro-js/popcount
function popcount(x) {
  x -= (x >>> 1) & 0x55555555;
  x = (x & 0x33333333) + ((x >>> 2) & 0x33333333);
  x = (x + (x >>> 4)) & 0x0f0f0f0f;
  x += x >>> 8;
  x += x >>> 16;
  return x & 0x7f;
}

// trailing0 = bit.trailing0
// Based on an implementation by @mikolalysenko: https://github.com/mikolalysenko/count-trailing-zeros
function trailing0(v) {
  var c = 32;
  v &= -v;
  if (v) c--;
  if (v & 0x0000ffff) c -= 16;
  if (v & 0x00ff00ff) c -= 8;
  if (v & 0x0f0f0f0f) c -= 4;
  if (v & 0x33333333) c -= 2;
  if (v & 0x55555555) c -= 1;
  return c;
}