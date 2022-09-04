// Bin: Axis-aligned rectangular region of space
// Tile: Power-of-2-sized and power-of-2-aligned bin in Morton space

export function tileIndexRanges(out, lo, hi, level, searchBefore, searchAfter) {
  // return an array of arrays of index ranges in z-order. The indices for individual tiles can be found by looking at differences between successive range elements
  rangeCheck(lo, hi, level);
  const shift = 2 * level; // bitshift factor to downscale [lo, hi] so that each tile is represented by a single z-order index
  const tileSize = 2 ** shift; // a tile at `level` is 2^level along each side, so contains 2^(2*level) elements
  let n = 0;
  let nRanges = 0;
  splitIntoContiguousRanges(lo >>> shift, hi >>> shift, (range) => {
    const loCode = (range.min << shift) >>> 0; // first code of the first tile in this range
    const hiCode = ((range.max << shift) >>> 0) + (tileSize - 1); // last code of the last tile
    out[n] = range.min; // encode range min
    out[n + 1] = range.max; // encode range max
    n += 2;
    // encode (min - max + 1) range elements (left tile edges)
    for (let code = loCode; code <= hiCode; code += tileSize) {
      out[n] = searchBefore(code);
      n += 1;
    }
    // encode final right tile edge
    out[n] = searchAfter(hiCode);
    n += 1;
    nRanges += 1;
  });
  return nRanges; // used to not need to reset the array after the fact, ie. can have junk but we know how many subranges
}

// todo: remove dependency on xMaxTiles, and simply treat fill the prefix of the `out` array
// as if (lo, hi) indicate the entire range
export function tileCounts(out, lo, hi, level, xMaxTiles, searchBefore, searchAfter) {
  rangeCheck(lo, hi, level);
  const shift = 2 * level;
  const tileSize = 2 ** shift; // a tile at `level` is 2^level along each side
  const loCodeScaled = lo >>> shift;
  const hiCodeScaled = hi >>> shift;
  const xTileLo = decode2x(loCodeScaled);
  const xTileHi = decode2y(loCodeScaled);
  splitIntoContiguousRanges(loCodeScaled, hiCodeScaled, (range) => {
    const loCode = (range.min << shift) >>> 0; // first code of the first tile in this range
    const hiCode = ((range.max << shift) >>> 0) + (tileSize - 1); // last code of the last tile
    let tileCode = range.min;
    let loIndex = searchBefore(loCode);
    for (
      let code = loCode + tileSize;
      code <= hiCode;
      code += tileSize, tileCode += 1
    ) {
      const hiIndex = searchBefore(code);
      const x = decode2x(tileCode) - xTileLo;
      const y = decode2y(tileCode) - xTileHi;
      const tileIndex = xMaxTiles * y + x;
      out[tileIndex] = hiIndex - loIndex;
      loIndex = hiIndex;
    }
    const hiIndex = searchAfter(hiCode);
    const x = decode2x(tileCode) - xTileLo;
    const y = decode2y(tileCode) - xTileHi;
    const tileIndex = xMaxTiles * y + x;
    out[tileIndex] = hiIndex - loIndex;
  });
  return out;
}

// If these requirements don't hold, there is no way to completely cover
// the given [lo, hi] interval with the tiles of the chosen power of 2.
function rangeCheck (lo, hi, level) {
  const step = 2 ** level;
  const [xLo, yLo] = decode2(lo);
  const [xHi, yHi] = decode2(hi + 1);
  if (!(xLo % step === 0 && yLo % step === 0))
    throw new Error("lo should evenly divide step");
  if (!(xHi % step === 0 && yHi % step === 0))
    throw new Error("hi+1 should evenly divide step");
}

// Find the leftmost insertion index in A for T
// in order to maintain A's sorted order.
// Searches the index range [L, R).
export function binarySearchBefore(A, T, L, R) {
  while (L < R) {
    // This midpoint calculation will return incorrect results for large arrays (>2^30)
    // By that point we should switch to Zig. Correct alternative: Math.floor(L + (R - L) / 2);
    const m = (L + R) >>> 1;
    if (A[m] < T) L = m + 1;
    else R = m;
  }
  return L;
}

// Find the rightmost insertion index in A for T
// in order to maintain A's sorted order.
// Searches the index range [L, R).
export function binarySearchAfter(A, T, L, R) {
  while (L < R) {
    // This midpoint calculation will return incorrect results for large arrays (>2^30)
    // By that point we should switch to Zig. Correct alternative: Math.floor(L + (R - L) / 2);
    const m = (L + R) >>> 1 
    if (A[m] > T) R = m;
    else L = m + 1;
  }
  return R;
}

export function encode2(x, y) {
  return  ((Part1By1(y) << 1) + Part1By1(x)) >>> 0
}

export function decode2x(code) {
  return  Compact1By1(code >> 0)
}

export function decode2y(code) {
  return  Compact1By1(code >> 1)
}

export function decode2(d) {
  return  [decode2x(d), decode2y(d)]
}

// Morton code subroutines based on Fabian Geisen's blog post:
// https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
// "Insert" a 0 bit after each of the 16 low bits of x
function Part1By1(x) {
  x &= 0x0000ffff; // x = ---- ---- ---- ---- fedc ba98 7654 3210
  x = (x ^ (x << 8)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
  x = (x ^ (x << 4)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
  x = (x ^ (x << 2)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
  x = (x ^ (x << 1)) & 0x55555555; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
  return x;
}

function Compact1By1(x) {
  x &= 0x55555555; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
  x = (x ^ (x >> 1)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
  x = (x ^ (x >> 2)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
  x = (x ^ (x >> 4)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
  x = (x ^ (x >> 8)) & 0x0000ffff; // x = ---- ---- ---- ---- fedc ba98 7654 3210
  return x;
}

// From https://twitter.com/jonahharris/status/1337087177591820290/photo/1
// Used with permission from Jonah, who can't remember where he got it but
// says he obtained it under the BSD license.
// The basic idea is to determine the MSB, then split the range below that,
// using the common prefix together with  a calculation for the new y/x
// positions indicating the split point.
// See also: https://snorrwe.onrender.com/posts/morton-table/#range-query-splitting
function litMaxBigMin(uMin, uMax) {
  const xor = uMin ^ uMax;
  const uMSBD = 1 << (31 - Math.clz32(xor));
  const xMask = 0x55555555;
  const yMask = ~xMask;
  const splitXAxis = uMSBD & xMask;
  const splitMask = splitXAxis ? xMask : yMask;
  const uMSMask = (uMSBD - 1) & splitMask;
  const uLSMask = (uMSBD - 1) & ~splitMask;
  const uBSCommon = uMin & ~(uMSBD + uMSBD - 1);
  const uLitMax = uBSCommon | uMSMask | (uLSMask & uMax);
  const uBigMin = uBSCommon | uMSBD | (uLSMask & uMin);
  return { litMax: uLitMax >>> 0, bigMin: uBigMin >>> 0 };
}

const stack = []; // preallocated stack for splitIntoContiguousRanges

// Computes the minimal set of contiguous ranges that cover the
// (possibly non-power-of-two) rectangle specified by the Morton
// codes `min` and `max`.
// Q: why does it always work to split up to three times but no more?
function splitIntoContiguousRanges(min, max, emit) {
  // assert that max > min
  stack.length = 0;
  stack.push({ min, max, iter: 0 }); // endpoints are inclusive
  let range = { min, max: min - 1 }; // start the first range at the passed-in min
  while (stack.length) {
    const { min, max, iter } = stack.pop();
    // if this range is still contiguous after splitting,
    // try splitting again up to two times; in all of my
    // testing this seems to guarantee that the resulting
    // range will lie within the original (min, max) bounds.
    // todo: find out if this is guaranteed to be the case!
    // (first time: iter = 0, so we split again;
    // second time: iter = 1, so we split again;
    // the third time, iter = 2 so we emit the range.)
    if (iter < 2 && max - min > 0) {
      const { litMax, bigMin } = litMaxBigMin(min, max);
      if (litMax + 1 == bigMin) {
        stack.push({ min: bigMin, max, iter: iter + 1 });
        stack.push({ min, max: litMax, iter: iter + 1 });
      } else {
        stack.push({ min: bigMin, max, iter: 0 });
        stack.push({ min, max: litMax, iter: 0 });
      }
    } else {
      // coalesce contiguous ranges as they are emitted
      if (range.max + 1 == min) {
        range.max = max;
      } else {
        emit(range);
        range = { min, max };
      }
    }
  }
  // emit the final range being built up in the loop
  if (range.max >= range.min) emit(range);
  // emit a final range up to the passed-in max
  if (range.max < max) emit({ min: range.max + 1, max });
};
