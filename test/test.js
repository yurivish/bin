import { BitVector } from './../BitVector.js';
import { CBitVector } from './../CBitVector.js';
import { ZeroCompressedBitVector } from './../ZeroCompressedBitVector.js';
import { NaiveBitVector } from './NaiveBitVector.js';
import assert from 'node:assert/strict';
import { readFileSync } from 'fs';

function testBitVector(name, construct, methods, destroy) {
  describe('bit vector: ' + name, function () {
    describe('constructor()', function () {
      it('should return a BitVector of the specified length', function () {
        for (const length of [0, 10, 1000]) {
          const v = construct(length, { rank: true, select: true });
          assert.equal(v.length, length);
          if (destroy) destroy(v);
        }
      });
    });

    describe('zero length', function () {
      const length = 0;
      const v = construct(length, { rank: true, select: true });
      v.finish();
      for (const j of [0, 100]) {
        if (methods.rank1) assert.equal(v.rank1(j), 0, `rank1(${j})`);
        if (methods.rank0) assert.equal(v.rank0(j), 0, `rank0(${j})`);
        if (methods.select1) assert.equal(v.select1(j), -1, `select1(${j})`);
        if (methods.select0) assert.equal(v.select0(j), -1, `select0(${j})`);
        assert.throws(() => v.access(j), `access(${j})`);
      }
      if (destroy) destroy(v);
    });

    describe('rank, select, and access', function () {
      // enumerate all permutations of nBits bits,
      // construct bitvectors for each, and
      // check all valid ranks, selects, and accesess
      // against NaiveBitVector. Also check specific
      // out-of-bounds values that are just beyond/far
      // beyond the valid values.
      const nBits = 8;
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
              const v = construct(length, { rank: true, select: true });
              const w = new NaiveBitVector(length);
              // for each 1 bit in i, add the respective one to the bitvector
              for (let j = 0; j < nBits; j++) {
                if (i & (1 << j)) {
                  v.one(j * spacing);
                  w.one(j * spacing);
                }
              }
              v.finish();
              w.finish();
              for (let j = 0; j < length; j++) {
                if (methods.rank1) assert.equal(v.rank1(j), w.rank1(j), `rank1(${j})`);
                if (methods.rank0) assert.equal(v.rank0(j), w.rank0(j), `rank0(${j})`);
                if (methods.select1) assert.equal(v.select1(j), w.select1(j), `select1(${j})`);
                if (methods.select0) assert.equal(v.select0(j), w.select0(j), `select0(${j})`);
                assert.equal(v.access(j), w.access(j)), `access(${j})`;
              }
              // test out-of-bounds results
              for (const j of [-1000, -1, length + 1, length + 1000]) {
                if (methods.rank1) assert.equal(v.rank1(j), w.rank1(j), `rank1(${j})`);
                if (methods.rank0) assert.equal(v.rank0(j), w.rank0(j), `rank0(${j})`);
                if (methods.select1) assert.equal(v.select1(j), w.select1(j), `select1(${j})`);
                if (methods.select0) assert.equal(v.select0(j), w.select0(j), `select0(${j})`);
                assert.throws(() => v.access(j), `access(${j})`);
                assert.throws(() => w.access(j), `access(${j})`);
              }
              if (destroy) destroy(v);
            });
          }
        }
      }
    });
  });
}

testBitVector('BitVector', (length, opts) => new BitVector(length, opts), {
  rank1: true,
  rank0: true,
  select1: true,
  select0: true,
});

testBitVector('ZeroCompressedBitVector', (length, opts) => new ZeroCompressedBitVector(length, opts), {
  rank1: true,
  rank0: true,
  select1: true,
  select0: false, // not implemented
});

const c_bitvector_wasm = readFileSync('./dist/bitvector.wasm');
const C = await WebAssembly.instantiate(c_bitvector_wasm).then((r) => r.instance.exports);

testBitVector(
  'CBitVector',
  (length, opts) => new CBitVector(length, opts, C),
  {
    rank1: true,
    rank0: true,
    select1: true,
    select0: true,
  },
  (v) => v.destroy(),
);
