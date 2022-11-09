// util.js
function popcount(x) {
  x -= x >>> 1 & 1431655765;
  x = (x & 858993459) + (x >>> 2 & 858993459);
  x = x + (x >>> 4) & 252645135;
  x += x >>> 8;
  x += x >>> 16;
  return x & 127;
}
function trailing0(v) {
  var c2 = 32;
  v &= -v;
  if (v)
    c2--;
  if (v & 65535)
    c2 -= 16;
  if (v & 16711935)
    c2 -= 8;
  if (v & 252645135)
    c2 -= 4;
  if (v & 858993459)
    c2 -= 2;
  if (v & 1431655765)
    c2 -= 1;
  return c2;
}
function reverseBits(v, numBits) {
  return reverseBits32(v) >>> 32 - numBits;
}
function reverseBits32(v) {
  v = v >>> 1 & 1431655765 | (v & 1431655765) << 1;
  v = v >>> 2 & 858993459 | (v & 858993459) << 2;
  v = v >>> 4 & 252645135 | (v & 252645135) << 4;
  v = v >>> 8 & 16711935 | (v & 16711935) << 8;
  v = v >>> 16 | v << 16;
  return v;
}
function enumerateOnes(block, callback) {
  block = block >>> 0;
  while (block != 0) {
    callback(trailing0(block));
    const t = block & -block;
    block = block ^ t;
  }
}
function onesArray(block, ret = []) {
  ret.length = 0;
  enumerateOnes(block, ret.push.bind(ret));
  return ret;
}

// BitVector.js
var BitVector = class {
  constructor(length) {
    const n = Math.ceil(length / 32);
    this.blocks = new Uint32Array(n);
    this.numOnes = 0;
    this.maxOnePosition = -1;
    this.length = length;
  }
  one(position) {
    if (position >= this.length)
      throw new Error("position must be < length");
    const blockIndex = position >>> 5;
    const bitOffset = position & 31;
    this.blocks[blockIndex] |= 1 << bitOffset;
    this.numOnes += 1;
    if (position > this.maxOnePosition)
      this.maxOnePosition = position;
  }
  finish() {
    this.storedLength = this.maxOnePosition + 1;
    const numBlocks = Math.ceil(this.storedLength / 32);
    if (numBlocks < this.blocks.length) {
      this.blocks = this.blocks.slice(0, numBlocks);
    }
    this.rankSuperblocks = new Uint32Array(numBlocks);
    this.rankSuperblocks[0] = popcount(this.blocks[0]);
    for (let i = 1; i < numBlocks; i++)
      this.rankSuperblocks[i] = this.rankSuperblocks[i - 1] + popcount(this.blocks[i]);
  }
  rank1(i) {
    if (i < 0)
      return 0;
    if (i >= this.storedLength)
      return this.numOnes;
    const blockIndex = i >>> 5;
    const rankSuperbock = this.rankSuperblocks[blockIndex];
    const block = this.blocks[blockIndex];
    const lowBitIndex = i & 31;
    const mask = 4294967294 << lowBitIndex;
    return rankSuperbock - popcount(block & mask);
  }
  rank0(i) {
    if (i < 0)
      return 0;
    if (i >= this.length)
      return this.length - this.numOnes;
    return i - this.rank1(i) + 1;
  }
  select1(i, L = 0, R = this.length) {
    if (i < 1)
      throw new Error("out of bounds: i < 1");
    if (i > this.numOnes) {
      throw new Error(`out of bounds: i (${i}) > numOnes (${this.numOnes})`);
    }
    while (L < R) {
      const m = L + R >>> 1;
      if (this.rank1(m) < i)
        L = m + 1;
      else
        R = m;
    }
    return L;
  }
  select0(i, L = 0, R = this.length) {
    if (i < 1)
      throw new Error("out of bounds: i < 1");
    const numZeros = this.length - this.numOnes;
    if (i > numZeros) {
      throw new Error(`out of bounds: i (${i}) > numZeros (${numZeros})`);
    }
    while (L < R) {
      const m = L + R >>> 1;
      if (this.rank0(m) < i)
        L = m + 1;
      else
        R = m;
    }
    return L;
  }
  access(i) {
    if (i < 0 || i > this.length)
      throw new Error("access: out of bounds at index " + i + " with length " + this.length);
    const blockIndex = i >>> 5;
    const bitOffset = i & 31;
    const block = this.blocks[blockIndex];
    const targetMask = 1 << bitOffset;
    return +Boolean(block & targetMask);
  }
  approxSizeInBits() {
    const blockBits = 8 * this.blocks.length * this.blocks.BYTES_PER_ELEMENT;
    const rankSuperblockBits = 8 * this.rankSuperblocks.length * this.rankSuperblocks.BYTES_PER_ELEMENT;
    return blockBits + rankSuperblockBits;
  }
};

// CBitVector.js
var c = null;
var CBitVector = class {
  constructor(length) {
    this.v = c.bitvector_init(length);
  }
  destroy() {
    c.bitvector_free(this.v);
  }
  one(i) {
    return c.bitvector_one(this.v, i | 0);
  }
  finish() {
    return c.bitvector_finish(this.v);
  }
  rank1(i) {
    return c.bitvector_rank1(this.v, i | 0);
  }
  rank0(i) {
    return c.bitvector_rank0(this.v, i | 0);
  }
  select1(i) {
    return c.bitvector_select1(this.v, i | 0);
  }
  select0(i) {
    return c.bitvector_select0(this.v, i | 0);
  }
  access(i) {
    return c.bitvector_access(this.v, i | 0);
  }
  approxSizeInBits() {
    return -1;
  }
};

// ZeroCompressedBitVector.js
var ZeroCompressedBitVector = class {
  constructor(length, { rank = false, select = false } = {}) {
    const initialNumBlocks = Math.ceil(length / 32);
    this.blocks = new Uint32Array(initialNumBlocks);
    if (select) {
      this.precedingZeroBlocks = new Uint32Array(initialNumBlocks);
      this.selectSuperblocks = new Uint32Array(initialNumBlocks);
    }
    if (rank) {
      this.isZeroBlock = new BitVector(initialNumBlocks);
    }
    this.numOnes = 0;
    this.numZeroBlocks = 0;
    this.prevUncompressedBlockIndex = -1;
    this.prevOnePosition = -1;
    this.length = length;
    this.rank = rank;
    this.select = select;
  }
  one(position) {
    if (position < 0)
      throw new Error("position must be >= 0");
    if (position >= this.length)
      throw new Error("position must be < length");
    if (position <= this.prevOnePosition)
      throw new Error("ones must be added in strictly-ascending order");
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
      if (this.select)
        this.precedingZeroBlocks[blockIndex] = this.numZeroBlocks;
    } else {
      blockIndex = this.prevUncompressedBlockIndex - this.numZeroBlocks;
    }
    const bitOffset = position & 31;
    if (this.select) {
      if (this.numOnes % 32 === 0) {
        const selectSampleIndex = this.numOnes >>> 5;
        const adjustment = popcount(this.blocks[blockIndex]);
        this.selectSuperblocks[selectSampleIndex] = blockIndex << 5 | adjustment;
      }
    }
    this.blocks[blockIndex] |= 1 << bitOffset;
    this.prevBlockIndex = blockIndex;
    this.numOnes += 1;
  }
  finish() {
    if (this.rank)
      this.isZeroBlock.finish();
    this.storedLength = this.prevOnePosition + 1;
    const numBlocks = Math.ceil(this.storedLength / 32) - this.numZeroBlocks;
    if (numBlocks < this.blocks.length) {
      this.blocks = this.blocks.slice(0, numBlocks);
      if (this.select)
        this.precedingZeroBlocks = this.precedingZeroBlocks.slice(0, numBlocks);
    }
    if (this.select) {
      const numSelectSuperblocks = Math.ceil(this.numOnes / 32);
      if (numSelectSuperblocks < this.selectSuperblocks.length) {
        this.selectSuperblocks = this.selectSuperblocks.slice(0, numSelectSuperblocks);
      }
    }
    if (this.rank || this.select) {
      this.rankSuperblocks = new Uint32Array(numBlocks);
      this.rankSuperblocks[0] = popcount(this.blocks[0]);
      for (let i = 1; i < numBlocks; i++) {
        this.rankSuperblocks[i] = this.rankSuperblocks[i - 1] + popcount(this.blocks[i]);
      }
    }
  }
  rank1(i) {
    if (!this.rank)
      throw new Error("no indexing structures for rank (see constructor)");
    if (i < 0)
      return 0;
    const uncompressedBlockIndex = i >>> 5;
    const numPrecedingZeroBlocks = this.isZeroBlock.rank1(uncompressedBlockIndex - 1);
    const isZeroBlock = this.isZeroBlock.access(uncompressedBlockIndex);
    i -= numPrecedingZeroBlocks << 5;
    if (i >= this.storedLength)
      return this.numOnes;
    const blockIndex = i >>> 5;
    const lowBitIndex = isZeroBlock ? 0 : i & 31;
    const rankSuperbock = this.rankSuperblocks[blockIndex];
    const block = this.blocks[blockIndex];
    const mask = 4294967294 << lowBitIndex;
    return rankSuperbock - popcount(block & mask);
  }
  access(i) {
    const uncompressedBlockIndex = i >>> 5;
    const numPrecedingZeroBlocks = this.isZeroBlock.rank1(uncompressedBlockIndex - 1);
    const isZeroBlock = this.isZeroBlock.access(uncompressedBlockIndex);
    i -= numPrecedingZeroBlocks << 5;
    if (i < 0 || i > this.length)
      throw new Error("access: out of bounds");
    const blockIndex = i >>> 5;
    const bitOffset = isZeroBlock ? 0 : i & 31;
    const block = this.blocks[blockIndex];
    const targetMask = 1 << bitOffset;
    return Boolean(block & targetMask);
  }
  rank0(i) {
    if (i < 0)
      return 0;
    if (i >= this.length)
      return this.length - this.numOnes;
    return i - this.rank1(i) + 1;
  }
  select1(i) {
    if (!this.select)
      throw new Error("no indexing structures for select (see constructor)");
    if (i < 1)
      throw new Error("out of bounds: i < 1");
    if (i > this.numOnes) {
      throw new Error(`out of bounds: i (${i}) > numOnes (${this.numOnes})`);
    }
    const selectSuperblockIndex = i - 1 >>> 5;
    const selectSuperblock = this.selectSuperblocks[selectSuperblockIndex];
    let blockIndex = selectSuperblock >>> 5;
    let blockRank = this.rankSuperblocks[blockIndex];
    const adjustment = selectSuperblock & 31;
    let prevBlockRank = (selectSuperblockIndex << 5) - adjustment;
    while (blockRank < i) {
      prevBlockRank = blockRank;
      blockIndex = blockIndex + 1;
      blockRank = this.rankSuperblocks[blockIndex];
    }
    let block = this.blocks[blockIndex];
    if (block === void 0)
      throw new Error("undef block");
    for (let r = prevBlockRank + 1; r < i; r++)
      block &= block - 1;
    const numPrecedingZeroBlocks = this.precedingZeroBlocks[blockIndex];
    return (numPrecedingZeroBlocks + blockIndex << 5) + trailing0(block);
  }
  approxSizeInBits() {
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
};

// MultiSet.js
var RankMultiSet = class {
  constructor(sortedData) {
    const data = sortedData;
    this.occupancy = new BitVector(data.length === 0 ? 0 : data[data.length - 1] + 1);
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
    if (!this.occupancy.rank1)
      throw new Error("occupancy must support rank1");
    if (!this.multiplicity.select1)
      throw new Error("multiplicity must support select1");
    if (i > this.maxValue)
      return this.length;
    const n = this.occupancy.rank1(i);
    if (n === 0)
      return 0;
    return this.multiplicity.select1(n) + 1;
  }
  approxSizeInBits() {
    return this.occupancy.approxSizeInBits() + this.multiplicity.approxSizeInBits();
  }
};
var AccessMultiSet = class {
  constructor(sortedData) {
    const data = sortedData;
    this.occupancy = new ZeroCompressedBitVector(data.length === 0 ? 0 : data[data.length - 1] + 1, { select: true });
    this.multiplicity = new BitVector(data.length);
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
    if (!this.occupancy.select1)
      throw new Error("occupancy must support select1");
    if (!this.multiplicity.rank1)
      throw new Error("multiplicity must support rank1");
    if (i < 0)
      throw new Error("i must be >= 0");
    if (i >= this.length)
      throw new Error("i must be < length");
    const n = this.multiplicity.rank1(i);
    return this.occupancy.select1(n);
  }
  approxSizeInBits() {
    return this.occupancy.approxSizeInBits() + this.multiplicity.approxSizeInBits();
  }
};

// NaiveBitVector.js
var NaiveBitVector = class {
  constructor(ones, length) {
    this.ones = ones;
    this.set = new Set(ones);
    this.sum = new Uint32Array(length);
    this.sum[0] = this.set.has(0);
    for (let i = 1; i < length; i++) {
      this.sum[i] = this.sum[i - 1] + this.set.has(i);
    }
  }
  rank1(i) {
    return this.sum[i];
  }
  select1(i) {
    return this.ones[i - 1];
  }
  access(i) {
    return +this.set.has(i);
  }
};

// WaveletMatrix.js
var WaveletMatrix = class {
  constructor(data, alphabetSize, opts = {}) {
    const { largeAlphabet = alphabetSize > data.length } = opts;
    if (largeAlphabet)
      return this.constructLargeAlphabet(data, alphabetSize, opts);
    const numLevels = Math.ceil(Math.log2(alphabetSize));
    const maxLevel = numLevels - 1;
    const hist = new Uint32Array(2 ** numLevels);
    const borders = new Uint32Array(2 ** numLevels);
    const levels = new Array(numLevels);
    for (let i = 0; i < numLevels; i++) {
      levels[i] = new CBitVector(data.length);
    }
    const level = levels[0];
    const levelBitMask = 1 << maxLevel;
    for (let i = 0; i < data.length; i++) {
      const d = data[i];
      hist[d] += 1;
      if (d & levelBitMask)
        level.one(i);
    }
    for (let l = maxLevel; l > 0; l--) {
      const m = 2 ** l;
      for (let i = 0; i < m; i++) {
        hist[i] = hist[2 * i] + hist[2 * i + 1];
      }
      borders[0] = 0;
      for (let i = 1; i < m; i++) {
        const prevIndex = reverseBits(i - 1, l);
        borders[reverseBits(i, l)] = borders[prevIndex] + hist[prevIndex];
      }
      const level2 = levels[l];
      const levelBit = maxLevel - l;
      const levelBitMask2 = 1 << levelBit;
      const bitPrefixMask = 4294967294 << levelBit;
      const bitPrefixShift = levelBit + 1;
      for (let i = 0; i < data.length; i++) {
        const d = data[i];
        const nodeIndex = (d & bitPrefixMask) >>> bitPrefixShift;
        const p = borders[nodeIndex];
        borders[nodeIndex] += 1;
        if (d & levelBitMask2)
          level2.one(p);
      }
    }
    const numZeros = new Uint32Array(numLevels);
    for (let i = 0; i < numLevels; i++) {
      levels[i].finish();
      numZeros[i] = levels[i].rank0(levels[i].length);
    }
    this.levels = levels;
    this.alphabetSize = alphabetSize;
    this.numZeros = numZeros;
    this.numLevels = numLevels;
    this.maxLevel = maxLevel;
    this.length = data.length;
    const sz = Math.min(2 ** this.numLevels, this.alphabetSize);
    this.F = new Uint32Array(sz);
    this.L = new Uint32Array(sz);
    this.S = new Uint32Array(sz);
    this.C = new Uint32Array(sz);
    this.C2 = new Uint32Array(sz);
  }
  constructLargeAlphabet(data, alphabetSize, opts = {}) {
    let next = new Uint32Array(data.length);
    data = new Uint32Array(data);
    const numLevels = Math.ceil(Math.log2(alphabetSize));
    const maxLevel = numLevels - 1;
    const levels = new Array(numLevels);
    for (let i = 0; i < numLevels; i++) {
      levels[i] = new CBitVector(data.length);
    }
    const numZeros = new Uint32Array(numLevels);
    const walk = new ArrayWalker(0, data.length);
    for (let l = 0; l < maxLevel; l++) {
      const level2 = levels[l];
      const levelBit = maxLevel - l;
      const levelBitMask2 = 1 << levelBit;
      for (let i = 0; i < data.length; i++) {
        const d = data[i];
        if (d & levelBitMask2) {
          next[walk.nextBackIndex()] = d;
          level2.one(i);
        } else {
          next[walk.nextFrontIndex()] = d;
        }
        numZeros[l] = walk.frontIndex;
      }
      walk.reset(next);
      const tmp = data;
      data = next;
      next = tmp;
    }
    const level = levels[maxLevel];
    const levelBitMask = 1 << 0;
    for (let i = 0; i < data.length; i++) {
      if (data[i] & levelBitMask)
        level.one(i);
    }
    numZeros[maxLevel] = level.rank0(level.length);
    for (let l = 0; l < numLevels; l++)
      levels[l].finish();
    this.levels = levels;
    this.alphabetSize = alphabetSize;
    this.numZeros = numZeros;
    this.numLevels = numLevels;
    this.maxLevel = maxLevel;
    this.length = data.length;
    const sz = Math.min(2 ** this.numLevels, this.alphabetSize, this.length);
    this.F = new Uint32Array(sz);
    this.L = new Uint32Array(sz);
    this.S = new Uint32Array(sz);
    this.C = new Uint32Array(sz);
    this.C2 = new Uint32Array(sz);
  }
  symbol(index) {
    if (index < 0 || index > this.length)
      throw new Error("symbol: out of bounds");
    let symbol = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const index1 = level.rank1(index - 1);
      if (level.access(index) === 0) {
        index = index - index1;
      } else {
        const nz = this.numZeros[l];
        index = nz + index1;
        const levelBitMask = 1 << this.maxLevel - l;
        symbol |= levelBitMask;
      }
    }
    return symbol;
  }
  countSymbol(first, last, symbol, groupBits) {
    const indices = this.symbolIndices(first, last, symbol, groupBits);
    return indices.last - indices.first;
  }
  find(first, last, symbol, n) {
    const indices = this.symbolIndices(first, last, symbol, 0);
    if (indices.first === indices.last)
      throw new Error("symbol does not appear in index range");
    let index = indices.first + n - 1;
    for (let l = this.numLevels; l > 0; ) {
      l -= 1;
      const level = this.levels[l];
      const nz = this.numZeros[l];
      if (index < nz) {
        const n2 = index + 1;
        index = level.select0(n2);
      } else {
        const n2 = index - nz + 1;
        index = level.select1(n2);
      }
    }
    return index;
  }
  symbolIndices(first, last, symbol, groupBits = 0) {
    const symbolGroupSize = 1 << groupBits;
    if (symbol % symbolGroupSize !== 0) {
      throw new Error("symbol must evenly divide the block size implied by groupBits");
    }
    if (symbol >= this.alphabetSize)
      throw new Error("symbol must be < alphabetSize");
    if (first > last)
      throw new Error("last must be <= first");
    if (first === last)
      return 0;
    if (first > this.length)
      throw new Error("first must be < wavelet matrix length");
    if (this.numLevels === 0)
      return 0;
    const numLevels = this.numLevels - groupBits;
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const first1 = level.rank1(first - 1);
      const last1 = level.rank1(last - 1);
      const levelBitMask = 1 << this.maxLevel - l;
      if ((symbol & levelBitMask) === 0) {
        first = first - first1;
        last = last - last1;
      } else {
        const nz = this.numZeros[l];
        first = nz + first1;
        last = nz + last1;
      }
    }
    return { first, last };
  }
  countLessThan(first, last, symbol) {
    if (first < 0)
      throw new Error("first must be >= 0");
    if (first > last)
      throw new Error("first must be <= last");
    if (last > this.length)
      throw new Error("last must be < wavelet matrix length");
    if (symbol <= 0)
      return 0;
    if (symbol >= this.alphabetSize)
      return last - first;
    let count = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const first1 = level.rank1(first - 1);
      const last1 = level.rank1(last - 1);
      const levelBitMask = 1 << this.maxLevel - l;
      if ((symbol & levelBitMask) === 0) {
        first = first - first1;
        last = last - last1;
      } else {
        count += last - last1 - (first - first1);
        const nz = this.numZeros[l];
        first = nz + first1;
        last = nz + last1;
      }
    }
    return count;
  }
  count(first, last, lower, upper) {
    return this.countLessThan(first, last, upper) - this.countLessThan(first, last, lower);
  }
  countSymbolBatch(first, last, sortedSymbols, groupBits = 0) {
    const symbolGroupSize = 1 << groupBits;
    for (const symbol of sortedSymbols) {
      if (symbol % symbolGroupSize !== 0)
        throw new Error("symbol must evenly divide the block size implied by groupBits");
    }
    if (first > last)
      throw new Error("first must be <= last");
    if (last > this.length)
      throw new Error("last must be < wavelet matrix length");
    const { F, L, C, S } = this;
    const walk = new ArrayWalker(sortedSymbols.length === 0 ? 0 : 1, F.length);
    const nextIndex = walk.nextFrontIndex();
    F[nextIndex] = first;
    L[nextIndex] = last;
    C[nextIndex] = sortedSymbols.length;
    S[nextIndex] = 0;
    walk.reset(F, L, C, S);
    let nRankCalls = 0;
    const numLevels = this.numLevels - groupBits;
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBitMask = 1 << this.maxLevel - l;
      const symbolsPerNode = (levelBitMask << 1) - 1;
      const symbolBitMask = 4294967295 << this.maxLevel - l;
      let k = 0;
      for (let i = 0; i < walk.len; i++) {
        const first2 = F[i];
        const first1 = level.rank1(first2 - 1);
        const first0 = first2 - first1;
        const last2 = L[i];
        const last1 = level.rank1(last2 - 1);
        const last0 = last2 - last1;
        nRankCalls += 2;
        const symbolCount = C[i];
        const symbol = S[i];
        const a = symbol & symbolBitMask;
        const b = a + symbolsPerNode;
        const m = a + b >>> 1;
        const lo = k;
        const hi = k + symbolCount;
        const splitIndex = binarySearchAfter(sortedSymbols, m, lo, hi);
        const numGoLeft = splitIndex - lo;
        const numGoRight = symbolCount - numGoLeft;
        k += symbolCount;
        if (numGoLeft > 0) {
          const nextIndex2 = walk.nextBackIndex();
          F[nextIndex2] = first0;
          L[nextIndex2] = last0;
          C[nextIndex2] = numGoLeft;
          S[nextIndex2] = symbol;
        }
        if (numGoRight > 0) {
          const nextIndex2 = walk.nextBackIndex();
          F[nextIndex2] = nz + first1;
          L[nextIndex2] = nz + last1;
          C[nextIndex2] = numGoRight;
          S[nextIndex2] = symbol | levelBitMask;
        }
      }
      walk.reset(F, L, C, S);
    }
    for (let i = 0; i < walk.len; i++)
      L[i] -= F[i];
    const counts = L.subarray(0, walk.len).slice();
    const symbols = S.subarray(0, walk.len).slice();
    return { symbols, counts, nRankCalls };
  }
  subcodeIndicator(subcodeSizesInBits) {
    let indicator = 0;
    let offset = 0;
    for (const sz of subcodeSizesInBits) {
      if (sz === 0)
        throw "cannot have zero-sized field";
      indicator |= 1 << sz - 1 + offset;
      offset += sz;
    }
    return indicator >>> 0;
  }
  encodeSubcodes(indicator, values) {
    if (indicator === 0) {
      if (values.length !== 1) {
        throw new Error("number of values must be one if the indicator is zero");
      }
      return values[0];
    }
    if (popcount(indicator) !== values.length) {
      throw new Error("number of values must be equal to the number of 1 bits in the indicator");
    }
    let code = 0;
    let offset = 0;
    let i = 0;
    while (indicator > 0) {
      const subcodeSize = trailing0(indicator) + 1;
      code |= values[i] << offset;
      i += 1;
      offset += subcodeSize;
      indicator >>>= subcodeSize;
    }
    return code >>> 0;
  }
  counts(first, last, lower, upper, { groupBits = 0, subcodeIndicator = 0, sort = true } = {}) {
    const symbolGroupSize = 1 << groupBits;
    if (lower % symbolGroupSize !== 0)
      throw new Error("lower must evenly divide the symbol block size implied by groupBits");
    if ((upper + 1) % symbolGroupSize !== 0)
      throw new Error("(upper + 1) must evenly divide the symbol block size implied by groupBits");
    const numLevels = this.numLevels - groupBits;
    const { F, L, S } = this;
    const walk = new ArrayWalker(1, F.length);
    const nextIndex = walk.nextFrontIndex();
    F[nextIndex] = first;
    L[nextIndex] = last;
    S[nextIndex] = 0;
    walk.reset(F, L, S);
    let nRankCalls = 0;
    let subcodeMask = 0;
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBit = this.maxLevel - l;
      const levelBitMask = 1 << levelBit;
      if ((subcodeIndicator & levelBitMask) === 0)
        subcodeMask |= levelBitMask;
      else
        subcodeMask = levelBitMask;
      const subcodeLower = lower & subcodeMask;
      const subcodeUpper = upper & subcodeMask;
      for (let i = 0; i < walk.len; i++) {
        const first2 = F[i];
        const first1 = level.rank1(first2 - 1);
        const first0 = first2 - first1;
        const last2 = L[i];
        const last1 = level.rank1(last2 - 1);
        const last0 = last2 - last1;
        const symbol = S[i];
        nRankCalls += 2;
        const leftCount = last0 - first0;
        if (leftCount > 0) {
          const a = symbol & subcodeMask;
          const b = (a | levelBitMask - 1) & subcodeMask;
          if (intervalsOverlapInclusive(a, b, subcodeLower, subcodeUpper)) {
            const nextIndex2 = sort ? walk.nextBackIndex() : walk.nextFrontIndex();
            F[nextIndex2] = first0;
            L[nextIndex2] = last0;
            S[nextIndex2] = symbol;
          }
        }
        const rightCount = last1 - first1;
        if (rightCount > 0) {
          const a = (symbol | levelBitMask) & subcodeMask;
          const b = (a | levelBitMask - 1) & subcodeMask;
          if (intervalsOverlapInclusive(a, b, subcodeLower, subcodeUpper)) {
            const nextIndex2 = walk.nextBackIndex();
            F[nextIndex2] = nz + first1;
            L[nextIndex2] = nz + last1;
            S[nextIndex2] = symbol | levelBitMask;
          }
        }
      }
      walk.reset(F, L, S);
    }
    for (let i = 0; i < walk.len; i++)
      L[i] -= F[i];
    const counts = L.subarray(0, walk.len).slice();
    const symbols = S.subarray(0, walk.len).slice();
    return { symbols, counts, nRankCalls };
  }
  quantile(first, last, index) {
    if (first > last)
      throw new Error("first must be <= last");
    if (last > this.length)
      throw new Error("last must be < wavelet matrix length");
    if (index < 0 || index >= last - first)
      throw new Error("index cannot be less than zero or exceed length of range [first, last)");
    let symbol = 0;
    let nRankCalls = 0;
    for (let l = 0; l < this.numLevels; l++) {
      const level = this.levels[l];
      const first0 = level.rank0(first - 1);
      const last0 = level.rank0(last - 1);
      const count = last0 - first0;
      nRankCalls += 2;
      if (index < count) {
        first = first0;
        last = last0;
      } else {
        const nz = this.numZeros[l];
        first = nz + (first - first0);
        last = nz + (last - last0);
        const levelBitMask = 1 << this.maxLevel - l;
        symbol |= levelBitMask;
        index -= count;
      }
    }
    return { symbol, count: last - first, nRankCalls };
  }
  quantileBatch(first, last, sortedIndices, groupBits = 0) {
    const symbolGroupSize = 1 << groupBits;
    if (first > last)
      throw new Error("first must be <= last");
    if (last > this.length)
      throw new Error("last must be < wavelet matrix length");
    for (let i = 1; i < sortedIndices.length; i++) {
      if (!(sortedIndices[i - 1] <= sortedIndices[i]))
        throw new Error("sorted indices must be sorted");
    }
    if (sortedIndices[0] < 0 || sortedIndices[sortedIndices.length - 1] >= last - first)
      throw new Error("sortedIndex cannot be less than zero or exceed length of range [first, last)");
    const { F, L, S, C } = this;
    const walk = new ArrayWalker(sortedIndices.length === 0 ? 0 : 1, F.length);
    const nextIndex = walk.nextFrontIndex();
    F[nextIndex] = first;
    L[nextIndex] = last;
    S[nextIndex] = 0;
    C[nextIndex] = sortedIndices.length;
    walk.reset(F, L, S, C);
    const I = this.C2.subarray(0, sortedIndices.length);
    I.set(sortedIndices);
    let nRankCalls = 0;
    const numLevels = this.numLevels - groupBits;
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBitMask = 1 << this.maxLevel - l;
      let k = 0;
      for (let i = 0; i < walk.len; i++) {
        const first2 = F[i];
        const first1 = level.rank1(first2 - 1);
        const first0 = first2 - first1;
        const last2 = L[i];
        const last1 = level.rank1(last2 - 1);
        const last0 = last2 - last1;
        const leftChildCount = last0 - first0;
        const symbol = S[i];
        const sortedIndexCount = C[i];
        nRankCalls += 2;
        const lo = k;
        const hi = k + sortedIndexCount;
        const splitIndex = binarySearchBefore(I, leftChildCount, lo, hi);
        const numGoLeft = splitIndex - lo;
        const numGoRight = sortedIndexCount - numGoLeft;
        k += sortedIndexCount;
        for (let n = splitIndex; n < hi; n++)
          I[n] -= leftChildCount;
        if (numGoLeft > 0) {
          const nextIndex2 = walk.nextBackIndex();
          F[nextIndex2] = first0;
          L[nextIndex2] = last0;
          S[nextIndex2] = symbol;
          C[nextIndex2] = numGoLeft;
        }
        if (numGoRight > 0) {
          const nextIndex2 = walk.nextBackIndex();
          F[nextIndex2] = nz + (first2 - first0);
          L[nextIndex2] = nz + (last2 - last0);
          S[nextIndex2] = symbol | levelBitMask;
          C[nextIndex2] = numGoRight;
        }
      }
      walk.reset(F, L, S, C);
    }
    for (let i = 0; i < walk.len; i++)
      L[i] -= F[i];
    const counts = L.subarray(0, walk.len).slice();
    const symbols = S.subarray(0, walk.len).slice();
    const numSortedIndices = C.subarray(0, walk.len).slice();
    return { symbols, counts, numSortedIndices, nRankCalls };
  }
  quantiles(first, last, firstIndex, lastIndex, groupBits = 0) {
    const symbolGroupSize = 1 << groupBits;
    if (first > last)
      throw new Error("first must be <= last");
    if (last > this.length)
      throw new Error("last must be < wavelet matrix length");
    if (firstIndex > lastIndex)
      throw new Error("firstIndex must be <= lastIndex");
    if (firstIndex < 0 || lastIndex > last - first)
      throw new Error("sortedIndex cannot be less than zero or exceed length of range [first, last)");
    const { F, L, S, C, C2 } = this;
    const walk = new ArrayWalker(firstIndex === lastIndex ? 0 : 1, F.length);
    const nextIndex = walk.nextFrontIndex();
    F[nextIndex] = first;
    L[nextIndex] = last;
    S[nextIndex] = 0;
    C[nextIndex] = firstIndex;
    C2[nextIndex] = lastIndex;
    walk.reset(F, L, S, C, C2);
    let nRankCalls = 0;
    const numLevels = this.numLevels - groupBits;
    for (let l = 0; l < numLevels; l++) {
      const level = this.levels[l];
      const nz = this.numZeros[l];
      const levelBitMask = 1 << this.maxLevel - l;
      for (let i = 0; i < walk.len; i++) {
        const first2 = F[i];
        const first1 = level.rank1(first2 - 1);
        const first0 = first2 - first1;
        const last2 = L[i];
        const last1 = level.rank1(last2 - 1);
        const last0 = last2 - last1;
        const symbol = S[i];
        const c2 = C[i];
        const c22 = C2[i];
        nRankCalls += 2;
        const leftChildCount = last0 - first0;
        if (c2 < leftChildCount) {
          const nextIndex2 = walk.nextBackIndex();
          F[nextIndex2] = first0;
          L[nextIndex2] = last0;
          S[nextIndex2] = symbol;
          C[nextIndex2] = c2;
          C2[nextIndex2] = leftChildCount < c22 ? leftChildCount : c22;
        }
        if (c22 > leftChildCount) {
          const nextIndex2 = walk.nextBackIndex();
          F[nextIndex2] = nz + (first2 - first0);
          L[nextIndex2] = nz + (last2 - last0);
          S[nextIndex2] = symbol | levelBitMask;
          C[nextIndex2] = leftChildCount < c2 ? c2 - leftChildCount : 0;
          C2[nextIndex2] = c22 - leftChildCount;
        }
      }
      walk.reset(F, L, S, C, C2);
    }
    for (let i = 0; i < walk.len; i++)
      L[i] -= F[i];
    const counts = L.subarray(0, walk.len).slice();
    const symbols = S.subarray(0, walk.len).slice();
    return { symbols, counts, nRankCalls };
  }
  majority(first, last) {
    const index = last + first >>> 1;
    const q = this.quantile(first, last, index);
    const count = last - first;
    const half = count >>> 1;
    if (q.count > half)
      return q.symbol;
    return null;
  }
  approxSizeInBits() {
    let size = 0;
    for (const level of this.levels) {
      size += level.approxSizeInBits();
    }
    return size;
  }
  symbolPath(i) {
    if (i < 0 || i > this.length)
      throw new Error("symbol: out of bounds");
    let l = 0;
    let a = 0;
    let b = (1 << this.numLevels) - 1;
    const path = [];
    while (a !== b) {
      const level = this.levels[l];
      const m = a + b >>> 1;
      path.push({ index: i, bit: level.access(i) });
      if (level.access(i) === 0) {
        i = level.rank0(i - 1);
        b = m;
      } else {
        const nz = this.numZeros[l];
        i = nz + level.rank1(i - 1);
        a = m + 1;
      }
      l += 1;
    }
    return { symbol: a, path, virtualLeafIndex: i };
  }
};
function binarySearchAfter(A, T, L, R) {
  while (L < R) {
    const m = L + R >>> 1;
    if (A[m] > T)
      R = m;
    else
      L = m + 1;
  }
  return R;
}
function binarySearchBefore(A, T, L, R) {
  while (L < R) {
    const m = L + R >>> 1;
    if (A[m] < T)
      L = m + 1;
    else
      R = m;
  }
  return L;
}
var ArrayWalker = class {
  constructor(len, cap) {
    this.len = len;
    this.cap = cap;
    this.frontIndex = 0;
    this.backIndex = cap;
  }
  nextFrontIndex() {
    const index = this.frontIndex;
    this.frontIndex += 1;
    return index;
  }
  nextBackIndex() {
    return this.backIndex -= 1;
  }
  reset(...arrays) {
    for (let i = 0; i < arrays.length; i++) {
      const arr = arrays[i];
      arr.set(arr.subarray(this.backIndex, this.cap).reverse(), this.frontIndex);
    }
    this.len = this.frontIndex + (this.cap - this.backIndex);
    this.frontIndex = 0;
    this.backIndex = this.cap;
  }
};
function intervalsOverlapInclusive(aLo, aHi, bLo, bHi) {
  return aLo <= bHi && bLo <= aHi;
}
export {
  AccessMultiSet,
  BitVector,
  CBitVector,
  NaiveBitVector,
  RankMultiSet,
  WaveletMatrix,
  ZeroCompressedBitVector,
  enumerateOnes,
  onesArray,
  popcount,
  reverseBits,
  reverseBits32,
  trailing0
};
