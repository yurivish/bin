import { WaveletMatrix } from './../js/WaveletMatrix.js';
import { BitVector } from './../js/BitVector.js';
import { CBitVector } from './../js/CBitVector.js';
import { ZeroCompressedBitVector } from './../js/ZeroCompressedBitVector.js';
import { NaiveBitVector } from './NaiveBitVector.js';
import assert from 'node:assert/strict';
import { readFileSync } from 'fs';

const doTestBitVector = true; // ⚠️
const doTestWaveletMatrix = true; // ⚠️

// const data = [0, 1, 3, 7, 1, 5, 4, 2, 6, 3];
// const wm = new WaveletMatrix(data, 8);

// Wavelet matrix tests
// [x] test with sort true/false (check equality of arr.sort())
// [x] test with construction w/ largeAlphabet / not
// [x] test alphabetsize >> number of present symbols
// [x] test alphabetsize >> max symbol
// [x] counts sort
// [x] access
// [x] quantile
// [x] majority
// [x] find
// [x] counts w/ subcodes
// [x] counts w/ groupBits
// [x] countSymbolBatch with groupBits
// [x] subcodeIndicator
// [x] encodeSubcodes
// todo
// - annotate tests with comments stating what they're testing
// - test sparse alphabet construction (just that it completes successfully)
// - test that errors are thrown for OOB input indices
// - test multiplicity construction
// - test multiplicity functions
// - test error behavior for invalid or OOB first/last/lower/upper

function testWaveletMatrix(alphabetSizePadding, wmOpts) {
  describe('wavelet matrix', function () {
    const data = [0, 1, 3, 7, 1, 5, 4, 2, 6, 3];
    const wm = new WaveletMatrix(data, 7 + alphabetSizePadding, wmOpts);
    it('counts correctly', function () {
      assert.equal(wm.countSymbol(0, wm.length, 3), 2);
      assert.equal(wm.countSymbol(0, wm.length - 1, 3), 1);
      assert.equal(wm.countSymbol(3, wm.length - 1, 3), 0);
      assert.equal(wm.countSymbol(3, wm.length - 1, 5), 1);

      // test groupBits:
      assert.equal(wm.countSymbol(0, wm.length, 0, { groupBits: 3 }), data.length);
      assert.equal(wm.countSymbol(0, wm.length, 0, { groupBits: 2 }), 6);
      assert.equal(wm.countSymbol(0, wm.length, 4, { groupBits: 2 }), 4);

      assert.equal(wm.countLessThan(3, wm.length - 1, 5), 3);
      assert.equal(wm.countLessThan(0, wm.length, 8), 10);
      assert.equal(wm.countLessThan(0, wm.length, 7), 9);
      assert.equal(wm.count(0, wm.length, 1, 7), 8);
      assert.equal(wm.count(0, wm.length, 7, 7), 0);

      {
        const res = wm.countSymbolBatch(0, wm.length, [1, 2, 3]);
        assert.deepStrictEqual(res.symbols, new Uint32Array([1, 2, 3]));
        assert.deepStrictEqual(res.counts, new Uint32Array([2, 1, 2]));
      }

      {
        const res = wm.countSymbolBatch(5, wm.length, [1, 2, 3]);
        assert.deepStrictEqual(res.symbols, new Uint32Array([1, 2, 3]));
        assert.deepStrictEqual(res.counts, new Uint32Array([0, 1, 1]));
      }

      {
        const res = wm.countSymbolBatch(0, wm.length, [0], { groupBits: 3 });
        assert.deepStrictEqual(res.symbols, new Uint32Array([0]));
        assert.deepStrictEqual(res.counts, new Uint32Array([10]));
      }

      {
        const res = wm.countSymbolBatch(0, wm.length, [0], { groupBits: 2 });
        assert.deepStrictEqual(res.symbols, new Uint32Array([0]));
        assert.deepStrictEqual(res.counts, new Uint32Array([6]));
      }

      {
        // allows duplicate inputs
        const res = wm.countSymbolBatch(0, wm.length, [1, 1, 1, 2, 3]);
        assert.deepStrictEqual(res.symbols, new Uint32Array([1, 2, 3]));
        assert.deepStrictEqual(res.counts, new Uint32Array([2, 1, 2]));
      }

      {
        // test subcode queries
        const wm = new WaveletMatrix([0, 1, 2, 3, 4, 5, 6, 7], 7 + alphabetSizePadding, wmOpts);
        // 0b00001111
        // 0b00110011 <-- we're querying for elements with a 1 bit at this level
        // 0b01010101
        const res = wm.counts(1, wm.length, 0b010, 0b111, { subcodeIndicator: 0b111 });
        assert.deepStrictEqual(res.symbols, new Uint32Array([2, 3, 6, 7]));
        assert.deepStrictEqual(res.counts, new Uint32Array([1, 1, 1, 1]));
      }

      {
        // subcodes
        assert.throws(() => wm.subcodeIndicator([2, 0, 2]));
        assert.equal(wm.subcodeIndicator([2, 1, 2]), 0b10110);
        assert.equal(wm.subcodeIndicator([3, 1, 3]), 0b1001100);
        assert.equal(wm.subcodeIndicator([3, 1, 4]), 0b10001100);
      }

      for (const sort of [true, false]) {
        {
          const res = wm.counts(5, wm.length, 1, 3, { sort });
          assert.deepStrictEqual(res.symbols, new Uint32Array([2, 3]));
          assert.deepStrictEqual(res.counts, new Uint32Array([1, 1]));
        }

        {
          // error: groupBits must evenly divide lower
          assert.throws(() => wm.counts(5, wm.length, 1, 3, { sort, groupBits: 1 }));

          const res = wm.counts(5, wm.length, 0, 3, { sort, groupBits: 1 });
          assert.deepStrictEqual(res.symbols, new Uint32Array([2]));
          assert.deepStrictEqual(res.counts, new Uint32Array([2]));
        }

        for (const upper of [7, 8, 100]) {
          const res = wm.counts(0, wm.length, 0, upper, { sort });
          if (sort) {
            assert.deepStrictEqual(res.symbols, new Uint32Array([0, 1, 2, 3, 4, 5, 6, 7]));
            assert.deepStrictEqual(res.counts, new Uint32Array([1, 2, 1, 2, 1, 1, 1, 1]));
          } else {
            assert.deepStrictEqual(res.symbols.sort(), new Uint32Array([0, 1, 2, 3, 4, 5, 6, 7]));
            assert.deepStrictEqual(res.counts.sort(), new Uint32Array([1, 1, 1, 1, 1, 1, 2, 2]));
          }
        }
      }
    });
    it('quantiles correctly', function () {
      {
        const ret = wm.quantile(0, wm.length, 0);
        assert.equal(ret.symbol, 0);
        assert.equal(ret.count, 1);
      }
      {
        const ret = wm.quantile(0, wm.length, 0);
        assert.equal(ret.symbol, 0);
        assert.equal(ret.count, 1);
      }
      {
        const ret = wm.quantile(0, wm.length, 0);
        assert.equal(ret.symbol, 0);
        assert.equal(ret.count, 1);
      }
      {
        const ret = wm.quantile(1, wm.length, 0);
        assert.equal(ret.symbol, 1);
        assert.equal(ret.count, 2);
      }
      {
        const ret = wm.quantile(0, wm.length, 1);
        assert.equal(ret.symbol, 1);
        assert.equal(ret.count, 2);
      }
      {
        const ret = wm.quantile(0, wm.length, 9);
        assert.equal(ret.symbol, 7);
        assert.equal(ret.count, 1);
      }
      {
        const ret = wm.quantileBatch(0, wm.length, [0, 1, 9]);
        assert.deepStrictEqual(ret.symbols, new Uint32Array([0, 1, 7]));
        assert.deepStrictEqual(ret.counts, new Uint32Array([1, 2, 1]));
      }
      {
        const ret = wm.quantiles(0, wm.length, 0, 1 + 1); // todo: should upper index be inclusive?
        assert.deepStrictEqual(ret.symbols, new Uint32Array([0, 1]));
        assert.deepStrictEqual(ret.counts, new Uint32Array([1, 2]));
      }
      {
        const ret = wm.quantiles(0, wm.length, 0, 0)
        assert.deepStrictEqual(ret.symbols, new Uint32Array([]));
        assert.deepStrictEqual(ret.counts, new Uint32Array([]));
      }

      //
    });
    it('accesses correctly', function () {
      for (let i = 0; i < data.length; i++) {
        assert.equal(wm.access(i), data[i]);
      }
      assert.throws(() => wm.access(-1), `access(-1)`);
      assert.throws(() => wm.access(data.length), `access(data.length)`);
    });

    it('selects correctly', function () {
      // data: [0, 1, 3, 7, 1, 5, 4, 2, 6, 3];
      // init with a larger alphabet size so we can search for non-appearing symbols:
      const wm = new WaveletMatrix(data, 10 + alphabetSizePadding, wmOpts);
      assert.equal(wm.select(0, wm.length, 0, -1), -1);
      assert.equal(wm.select(0, wm.length, 0, 0), -1);
      assert.equal(wm.select(0, wm.length, 0, 1), 0);
      assert.equal(wm.select(0, wm.length, 1, 1), 1);
      assert.equal(wm.select(0, wm.length, 1, 0), -1);
      assert.equal(wm.select(0, wm.length, 3, 1), 2);
      assert.equal(wm.select(0, wm.length, 3, 2), 9);
      assert.equal(wm.select(3, wm.length, 3, 1), 9);
      assert.equal(wm.select(0, wm.length, 3, 3), -1);
      assert.equal(wm.select(0, wm.length, 0, 2), -1);
      assert.equal(wm.select(0, wm.length, 10, 1), -1);
      assert.equal(wm.select(1, wm.length, 0, 1), -1);
    });
    it('majorities correctly', function () {
      let wm = new WaveletMatrix(data, 10 + alphabetSizePadding, wmOpts);
      assert.equal(wm.simpleMajority(0, wm.length), null);
      assert.equal(wm.simpleMajority(0, 1).symbol, 0);
      assert.equal(wm.simpleMajority(0, 1).count, 1);
      wm = new WaveletMatrix([1, 1, 2, 1, 2, 2, 2], 10 + alphabetSizePadding, wmOpts);
      assert.equal(wm.simpleMajority(0, wm.length).symbol, 2);
      assert.equal(wm.simpleMajority(0, wm.length).count, 4);
      assert.equal(wm.simpleMajority(0, 4).symbol, 1);
      assert.equal(wm.simpleMajority(0, 4).count, 3);
      assert.deepStrictEqual(wm.majority(0, 4, 2), { symbols: new Uint32Array([1]), counts: new Uint32Array([3]) });
      assert.deepStrictEqual(wm.majority(0, 4, 1), { symbols: new Uint32Array([]), counts: new Uint32Array([]) });
      assert.deepStrictEqual(wm.majority(0, 4, 5), {
        symbols: new Uint32Array([1, 2]),
        counts: new Uint32Array([3, 1]),
      });
      assert.deepStrictEqual(wm.majority(0, 3, 4), {
        symbols: new Uint32Array([1, 2]),
        counts: new Uint32Array([2, 1]),
      });
      assert.deepStrictEqual(wm.majority(0, wm.length, 2), {
        symbols: new Uint32Array([2]),
        counts: new Uint32Array([4]),
      });
    });
  });
}

if (doTestWaveletMatrix) {
  // Ensure we're testing both wavelet matrix construction algorithms,
  // since which one to use would otherwise be decided by a heuristic.
  testWaveletMatrix(0, { largeAlphabet: false });
  testWaveletMatrix(0, { largeAlphabet: true });
  testWaveletMatrix(100, { largeAlphabet: false });
  testWaveletMatrix(100, { largeAlphabet: true });
} else {
  console.log('skipping wavelet matrix tests')
}


// Bit vector tests

const C = await WebAssembly.instantiate(readFileSync('./dist/bitvector.wasm')).then((r) => r.instance.exports);

function testBitVector(T) {
  describe(T.prototype.constructor.name, function () {
    it('constructor should return a BitVector of the specified length', function () {
      for (const length of [0, 10, 1000]) {
        const v = new T(length, { rank: true, select: true, C });
        assert.equal(v.length, length);
        if (v.destroy) v.destroy();
      }
    });

    it('should support zero-length bitvectors', function () {
      const length = 0;
      const v = new T(length, { rank: true, select: true, C });
      v.finish();
      for (const j of [0, 100]) {
        assert.equal(v.rank1(j), 0, `rank1(${j})`);
        assert.equal(v.rank0(j), 0, `rank0(${j})`);
        if (v.select1) assert.equal(v.select1(j), -1, `select1(${j})`);
        if (v.select0) assert.equal(v.select0(j), -1, `select0(${j})`);
        assert.throws(() => v.access(j), `access(${j})`);
      }
      if (v.destroy) v.destroy();
    });

    it('should work if the bitvector is all 0', function () {
      const length = 1e4;
      const v = new T(length, { rank: true, select: true, C });
      v.finish();
      for (let j = 0; j < length; j++) {
        assert.equal(v.rank1(j), 0, `rank1(${j})`);
        assert.equal(v.rank0(j), j + 1, `rank0(${j})`);
        assert.equal(v.select1(j), -1, `select1(${j})`);
        if (v.select0) assert.equal(v.select0(j + 1), j, `select0(${j})`);
        assert.equal(v.access(j), 0, `access(${j})`);
      }
      if (v.destroy) v.destroy();
    });

    it(`should work if the bitvector is all 1`, function () {
      const length = 1e4;
      const v = new T(length, { rank: true, select: true, C });
      for (let i = 0; i < length; i++) v.one(i);
      v.finish();
      for (let j = 0; j < length; j++) {
        assert.equal(v.rank1(j), j + 1, `rank1(${j})`);
        assert.equal(v.rank0(j), 0, `rank0(${j})`);
        assert.equal(v.select1(j + 1), j, `select1(${j})`);
        if (v.select0) assert.equal(v.select0(j), -1, `select0(${j})`);
        assert.equal(v.access(j), 1, `access(${j})`);
      }
      if (v.destroy) v.destroy();
    });

    it('should work when a 1-bit is on a block boundary', function () {
      const v = new ZeroCompressedBitVector(760, { rank: true, select: true, C });
      v.one(32);
      v.finish();
      assert.equal(v.rank1(0), 0, 'rank1(0)');
      assert.equal(v.access(0), 0, 'access(0)');
    });

    describe('rank, select, and access', function () {
      // enumerate all permutations of nBits bits,
      // construct bitvectors for each, and
      // check all valid ranks, selects, and accesess
      // against NaiveBitVector. Also check specific
      // out-of-bounds values that are just beyond/far
      // beyond the valid values.
      const nBits = 6;
      const limit = 2 ** nBits - 1;
      for (let i = 0; i < limit; i++) {
        for (let offset = 0; offset < 2; offset++) {
          // test ones at a starting offst
          for (let spacing = 1; spacing < 101; spacing += 25) {
            // test ones densely packed and spaced out
            const length = offset + nBits * spacing;
            it(`should work on vectors with offset ${offset}, spacing ${spacing}, length ${length}, and bits ${i} = b${i.toString(
              2,
            )}`, function () {
              // construct
              const v = new T(length, { rank: true, select: true, C });
              const w = new NaiveBitVector(length);
              for (let j = 0; j < nBits; j++) {
                // for each 1 bit in i, add the respective one to the bitvector
                if (i & (1 << j)) {
                  v.one(j * spacing);
                  w.one(j * spacing);
                }
              }
              v.finish();
              w.finish();
              for (let j = 0; j < length; j++) {
                // test in-bounds values
                assert.equal(v.rank1(j), w.rank1(j), `rank1(${j})`);
                assert.equal(v.rank0(j), w.rank0(j), `rank0(${j})`);
                assert.equal(v.select1(j), w.select1(j), `select1(${j})`);
                if (v.select0) assert.equal(v.select0(j), w.select0(j), `select0(${j})`);
                assert.equal(v.access(j), w.access(j)), `access(${j})`;
              }
              for (const j of [-1000, -1, length + 1, length + 1000]) {
                // test out-of-bounds results
                assert.equal(v.rank1(j), w.rank1(j), `rank1(${j})`);
                assert.equal(v.rank0(j), w.rank0(j), `rank0(${j})`);
                assert.equal(v.select1(j), w.select1(j), `select1(${j})`);
                if (v.select0) assert.equal(v.select0(j), w.select0(j), `select0(${j})`);
                assert.throws(() => v.access(j), `access(${j})`);
                assert.throws(() => w.access(j), `access(${j})`);
              }
              if (v.destroy) v.destroy();
            });
          }
        }
      }
    });

    // todo: deterministic random
    it('should return correct results in randomized tests', function () {
      const length = 1e5;
      for (let i = 0; i < 5; i++) {
        const v = new T(length, { rank: true, select: true, C });
        const nv = new NaiveBitVector(length, { rank: true, select: true, C });
        const n = Math.random();
        for (let j = 0; j < n; j++) {
          const i = Math.floor(length * Math.random());
          v.one(i);
          nv.one(i);
        }
        v.finish();
        nv.finish();
        for (let j = 0; j < length; j++) {
          assert.equal(v.rank1(j), nv.rank1(j), `rank1(${j})`);
          assert.equal(v.rank0(j), nv.rank0(j), `rank0(${j})`);
          assert.equal(v.select1(j), nv.select1(j), `select1(${j})`);
          assert.equal(v.access(j), nv.access(j), `access(${j})`);
        }
        if (v.destroy) v.destroy();
      }
    });
  });
}

if (doTestBitVector) {
  testBitVector(BitVector);
  testBitVector(CBitVector);
  testBitVector(ZeroCompressedBitVector);
} else {
  console.log('skipping bitvector tests')
}
