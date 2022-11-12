import { BitVector } from './../BitVector.js';
import { CBitVector } from './../CBitVector.js';
import { ZeroCompressedBitVector } from './../ZeroCompressedBitVector.js';
import { NaiveBitVector } from './NaiveBitVector.js';
import assert from 'node:assert/strict';
import { readFileSync } from 'fs';

const C = await WebAssembly.instantiate(readFileSync('./dist/bitvector.wasm')).then((r) => r.instance.exports);

function testBitVector(T) {
  describe('bit vector: ' + T.prototype.constructor.name, function () {
    describe('constructor()', function () {
      it('should return a BitVector of the specified length', function () {
        for (const length of [0, 10, 1000]) {
          const v = new T(length, { rank: true, select: true, C });
          assert.equal(v.length, length);
          if (v.destroy) v.destroy();
        }
      });
    });

    describe('zero length', function () {
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

    describe('all full, all empty', function () {
      const length = 1e4;
      it(`should work with all zeros`, function () {
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
      it(`should work with all ones`, function () {
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
    });

    describe('rank, select, and access', function () {
      it('works when a 1-bit is on a block boundary', function() {
        const v = new ZeroCompressedBitVector(760, { rank: true, select: true, C });
        v.one(32);
        v.finish();
        assert.equal(v.rank1(0), 0, 'rank1(0)')
        assert.equal(v.access(0), 0, 'access(0)')
      });

      // enumerate all permutations of nBits bits,
      // construct bitvectors for each, and
      // check all valid ranks, selects, and accesess
      // against NaiveBitVector. Also check specific
      // out-of-bounds values that are just beyond/far
      // beyond the valid values.
      const nBits = 7; 
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
  });

  // todo: deterministic random
  describe('randomized', function () {
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
}

testBitVector(BitVector);
testBitVector(CBitVector);
testBitVector(ZeroCompressedBitVector);