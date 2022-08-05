const std = @import("std");
const math = std.math;
const Vec = @import("vec.zig");

const inf64 = math.inf(f64);

fn sqr(x: f64) f64 {
    return x * x;
}

// Compute the one-dimensional Euclidean Distance Transform.
// Accepts and returns squared Euclidean distances.
// out: output slice; may alias input.
// in: input slice
// offset, stride, len: array iteration specification that allows encoding
// both row and column traversals of a slice representing a matrix.
// v: holds indices into f
// f: holds in array values for the current row/col (needed in order to operate in-place)
// z: holds envelope interval boundaries
// See "Distance Transforms of Sampled Functions" by Felzenszwalb & Huttenlocher
// for an explanation of the algorithms (variable names are theirs).
//
// Subtleties:
// We always compute `s`. It will be:
// * NaN when y1 == y2 == inf
// * inf when f[0] is inf
// * -inf when y1 == inf and y2 is not
// * inf when y2 == inf and y1 is not
// Since x1 is always less than x2, we know the denominator is always nonzero.
// Since f is always nonnegative, we know the y values are always nonnegative.
// The loop break condition has the k == 0 check so that we break when the value
// of `s` is -inf or NaN.
fn edt1d(out: []f64, in: []f64, offset: usize, stride: usize, len: usize, v: []usize, f: []f64, z: []f64) void {
    v[0] = 0;
    f[0] = in[offset];
    z[0] = -inf64;
    z[1] = inf64;
    var k: usize = 0;
    var s: f64 = 0;
    var q: usize = 1;
    while (q < len) : (q += 1) {
        f[q] = in[offset + q * stride];
        while (true) : (k -= 1) {
            const vk = v[k];
            const x1 = @intToFloat(f64, vk);
            const y1 = f[vk];
            const x2 = @intToFloat(f64, q);
            const y2 = f[q];
            s = (x1 * x1 - x2 * x2 + (y1 - y2)) / (2 * (x1 - x2));
            if (s > z[k] or k == 0) break;
        }
        k += 1;
        v[k] = q;
        z[k] = s;
        z[k + 1] = inf64;
    }
    k = 0;
    q = 0;
    while (q < len) : (q += 1) {
        const qf = @intToFloat(f64, q);
        while (z[k + 1] < @intToFloat(f64, q)) k += 1;
        const vk = v[k];
        out[offset + q * stride] = f[vk] + sqr(qf - @intToFloat(f64, vk));
    }
}

pub fn edt2d(out: []f64, in: []f64, width: usize, height: usize, v: []usize, f: []f64, z: []f64) void {
    // iterate over rows
    var y: usize = 0;
    while (y < in.len) : (y += width) {
        const offset = y;
        const stride = 1;
        const len = width;
        edt1d(out, in, offset, stride, len, v, f, z);
    }
    // iterate over columns
    var x: usize = 0;
    while (x < width) : (x += 1) {
        const offset = x;
        const stride = width;
        const len = height;
        edt1d(out, in, offset, stride, len, v, f, z);
    }
}

// One-dimensional binary Euclidean distance transform. Optionally in-place.
// Input: Array of indicator values with 0 denoting a filled
// sample and any other value denoting an unfilled sample.
// Note: Unlike the parabola-based method, this function
// treats all non-zero input values as positive infinity.
// Output: Array of sqruared euclidean distances from each
// array sample to its nearest filled neighbor.
pub fn edtBinary1d(out: []f64, in: []f64) void {
    // forwards sweep:
    // compute the minimum distance to each sample's nearest
    // left zero-valued neighbor
    var tmp = inf64;
    for (in) |x, i| {
        tmp = if (x == 0) 0 else tmp + 1;
        out[i] = sqr(tmp);
    }
    // backwards sweep:
    // take the min between the distance to the nearest
    // left zero-valued neighbor and the nearest right
    // zero-valued neighbor
    tmp = inf64;
    var i = out.len;
    while (i > 0) {
        i -= 1;
        const x = out[i];
        tmp = if (x == 0) 0 else tmp + 1;
        if (sqr(tmp) < x) out[i] = sqr(tmp);
    }
}

pub fn edtBinary2d(out: []f64, in: []f64, width: usize, height: usize, v: []usize, f: []f64, z: []f64) void {
    // iterate over rows, treating all nonzero values as infinity
    var y: usize = 0;
    while (y < in.len) : (y += width) {
        const outRow = out[y .. y + width];
        const inRow = in[y .. y + width];
        edtBinary1d(outRow, inRow);
    }
    // iterate over columns, treating all values as squared euclidean distances
    var x: usize = 0;
    while (x < width) : (x += 1) {
        const offset = x;
        const stride = width;
        const len = height;
        edt1d(out, in, offset, stride, len, v, f, z);
    }
}
