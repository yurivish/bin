// Based on an implementation by @ashaffer: https://github.com/micro-js/popcount
export function popcount(x) {
  x -= (x >>> 1) & 0x55555555;
  x = (x & 0x33333333) + ((x >>> 2) & 0x33333333);
  x = (x + (x >>> 4)) & 0x0f0f0f0f;
  x += x >>> 8;
  x += x >>> 16;
  return x & 0x7f;
}

// trailing0 = bit.trailing0
// Based on an implementation by @mikolalysenko: https://github.com/mikolalysenko/count-trailing-zeros
export function trailing0(v) {
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

// Reverse the lowest `numBits` bits of `v`.
// E.g. reverseBits(0b0000100100, 6) === 0b0000001001
//                       ^^^^^^               ^^^^^^
export function reverseBits(v, numBits) {
  return reverseBits32(v) >>> (32 - numBits);
}

// Adapted from https://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel
export function reverseBits32(v) {
  // unsigned int v; // 32-bit word to reverse bit order
  // swap odd and even bits
  v = ((v >>> 1) & 0x55555555) | ((v & 0x55555555) << 1);
  // swap consecutive pairs
  v = ((v >>> 2) & 0x33333333) | ((v & 0x33333333) << 2);
  // swap nibbles ...
  v = ((v >>> 4) & 0x0f0f0f0f) | ((v & 0x0f0f0f0f) << 4);
  // swap bytes
  v = ((v >>> 8) & 0x00ff00ff) | ((v & 0x00ff00ff) << 8);
  // swap 2-byte long pairs
  v = (v >>> 16) | (v << 16);
  return v;
}

export function clamp(x, lo, hi) {
  return x < lo ? lo : x > hi ? hi : x;
}

// unused in this library right now, but convenient
export function enumerateOnes(block, callback) {
  block = block >>> 0
  while (block != 0) {
    callback(trailing0(block));
    const t = block & -block;
    block = block ^ t;
  }
}

// unused in this library right now, but convenient
export function onesArray(block, ret = []) {
  ret.length = 0;
  enumerateOnes(block, ret.push.bind(ret));
  return ret;
}

// Utility taken from htl.
export function isObjectLiteral(value) {
  return value && value.toString === Object.prototype.toString;
}

// Returns the rightmost insertion index for T in A in order to maintain A's sorted order.
// Searches the index range [L, R).
// This variant takes an access function instead of an array.
function binarySearchAfterAccess(access, T, L, R) {
  while (L < R) {
    // Note: This midpoint calculation will return incorrect results for arrays with length > 2^30
    // Correct (but slow) alternative: Math.floor(L + (R - L) / 2)
    const m = (L + R) >>> 1;
    if (access(m) > T) R = m;
    else L = m + 1;
  }
  return R;
}


// Returns the rightmost insertion index for T in A in order to maintain A's sorted order.
// Searches the index range [L, R).
function binarySearchAfter(A, T, L, R) {
  while (L < R) {
    // Note: This midpoint calculation will return incorrect results for arrays with length > 2^30
    // Correct (but slow) alternative: Math.floor(L + (R - L) / 2)
    const m = (L + R) >>> 1;
    if (A[m] > T) R = m;
    else L = m + 1;
  }
  return R;
}

// Returns the leftmost insertion index for T in A in order to maintain A's sorted order.
// Searches the index range [L, R).
function binarySearchBefore(A, T, L, R) {
  while (L < R) {
    // Note: This midpoint calculation will return incorrect results for arrays with length > 2^30
    // Correct (but slow) alternative: Math.floor(L + (R - L) / 2)
    const m = (L + R) >>> 1;
    if (A[m] < T) L = m + 1;
    else R = m;
  }
  return L;
}