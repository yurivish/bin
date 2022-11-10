import { BitVector } from './../BitVector.js';
import { CBitVector } from './../CBitVector.js';
import { ZeroCompressedBitVector } from './../ZeroCompressedBitVector.js';
import { NaiveBitVector } from './NaiveBitVector.js';
import assert from 'node:assert/strict';

function testBitVector(BitVector, methods) {
  describe('BitVector', function () {
    describe('constructor()', function () {
      it('should return a BitVector of the specified length', function () {
        for (const length of [0, 10, 1000]) {
          const v = new BitVector(length, { rank: true, select: true });
          assert.equal(v.length, length);
        }
      });
    });
    describe('rank and select', function () {
      // enumerate all permutations of nBits bits,
      // construct bitvectors for each, and
      // check all ranks
      const nBits = 8;
      const limit = 2 ** nBits - 1;
      for (let i = 0; i < limit; i++) {
        if (i > 1) break;
        // todo: test offset too
        for (let spacing = 1; spacing < 100; spacing += 10) {
          it(`should work on vectors with spacing ${spacing} and bits ${i} (b${i.toString(2)})`, function () {
            // space out the ones
            const v = new BitVector(nBits * spacing, { rank: true, select: true });
            const w = new NaiveBitVector(nBits * spacing);
            for (let j = 0; j < nBits; j++) {
              // for each 1 bit in i, add the respective one to the bitvector
              if (i & j) {
                v.one(j * spacing);
                w.one(j * spacing);
              }
            }
            v.finish();
            w.finish();
            for (let j = 0; j < nBits * spacing; j++) {
              if (methods.rank1) assert.equal(v.rank1(j), w.rank1(j), `rank1(${j})`);
              if (methods.rank0) assert.equal(v.rank0(j), w.rank0(j), `rank0(${j})`);
              if (methods.select1) assert.equal(v.select1(j), w.select1(j), `select1(${j})`);
              if (methods.select0) assert.equal(v.select0(j), w.select0(j), `select0(${j})`);
              assert.equal(v.access(j), w.access(j)), `access(${j})`;
            }

            // // test out-of-bounds results
            // for (const j of [-1, 1e6]) {
            //   if (methods.rank1) assert.equal(v.rank1(j), w.rank1(j), `rank1(${j})`);
            //   if (methods.rank0) assert.equal(v.rank0(j), w.rank0(j), `rank0(${j})`);
            //   if (methods.select1) assert.equal(v.select1(j), w.select1(j), `select1(${j})`);
            //   if (v.select0) assert.equal(v.select0(j), w.select0(j), `select0(${j})`);
            //   assert.throws(() => v.access(j), `access(${j})`);
            //   assert.throws(() => w.access(j), `access(${j})`);
            // }
          });
        }
      }
    });
  });
}

testBitVector(BitVector, {
  rank1: true,
  rank0: true,
  select1: true,
  select0: true,
});

// todo: fix! has 🐛
// both rank1 and access need the fix, i think (rank1 calls access)
// has to do with the zero blocks and miscalculation of block index,
// possibly the errors are only in the case when blockIndex is one
// beyond the end of the array.
// testBitVector(ZeroCompressedBitVector, {
//   rank1: true,
//   rank0: true,
//   select1: true,
//   select0: false /* not implemented */,
// });

// todo: figure out how to load the wasm
// testBitVector(CBitVector)
