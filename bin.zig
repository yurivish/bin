const std = @import("std");
const math = std.math;
const Vec = @import("vec.zig");
const assert = std.debug.assert;

pub const Rgba = packed struct { r: u8, g: u8, b: u8, a: u8 };

pub fn colorize(colors: []Rgba, vs: anytype, ramp: []const Rgba) void {
    // Compute a color for each bin by looking up the appropriate ramp value
    const unknown = Rgba{ .r = 255, .g = 0, .b = 255, .a = 255 };
    for (vs.data) |v, i| {
        colors[i] = if (vs.inBounds(v)) ramp[vs.bin(v)] else unknown;
    }
}

pub fn count2d(bins: []u32, xs: anytype, ys: anytype) void {
    const xBins = xs.nBins;
    const len = xs.data.len;
    var i: usize = 0;
    while (i < len) : (i += 1) {
        const x = xs.binPcBounds(i);
        const y = ys.binPcBounds(i);
        if (x.inBounds and y.inBounds) {
            const i_bin = @floatToInt(usize, y.index * xBins + x.index);
            bins[i_bin] += 1;
        }
    }
}

fn f2u(index: f32) usize {
    return @floatToInt(usize, index);
}

pub fn count3d(bins: []u32, xs: anytype, ys: anytype, zs: anytype) void {
    // caution: untested
    const xBins = f2u(xs.nBins);
    const xyBins = f2u(xs.nBins * ys.nBins);
    const len = xs.data.len;
    var i: usize = 0;
    while (i < len) : (i += 1) {
        const x = xs.binPcBounds(i);
        const y = ys.binPcBounds(i);
        const z = zs.binPcBounds(i);

        // need to do all this if float32 since f32 cannot represent 4096^3

        if (x.inBounds and y.inBounds and z.inBounds) {
            const i_bin = xyBins * f2u(z.index) + xBins * f2u(y.index) + x.index;
            bins[i_bin] += 1;
        }
    }
}

// pub

// pub fn count2d(bins: []u32, xs: anytype, ys: anytype) void {
//     const xBins = xs.nBins;
//     const len = xs.data.len;
//     var i: usize = 0;
//     while (i < len) : (i += 1) {
//         var f_bin = ys.binPc(i) * xBins + xs.binPc(i);
//         if (!math.isNan(f_bin)) {
//             const i_bin = @floatToInt(usize, f_bin);
//             bins[i_bin] += 1;
//         }
//     }
// }

// pub fn bin1d(bins: []f64, assignments: []usize, xs: Vec) void {
//     for (xs.data) |x, i| {
//         if (xs.inBounds(x)) {
//             const index = xs.bin(x);
//             bins[index] += 1;
//             assignments[i] = index;
//         } else {
//             assignments[i] = math.maxInt(usize);
//         }
//     }
// }

// pub fn bin2d(bins: []f64, assignments: []usize, xs: Vec, ys: Vec) void {
//     const xBins = xs.nBins;
//     for (xs.data) |x, i| {
//         const y = ys.data[i];
//         if ((xs.inBounds(x) and ys.inBounds(y))) {
//             const index = xBins * ys.bin(y) + xs.bin(x);
//             bins[index] += 1;
//             assignments[i] = index;
//         } else {
//             assignments[i] = math.maxInt(usize);
//         }
//     }
// }

// pub fn count(comptime T: type, bins: []u32, vecs: anytype) void {
//     const len = vecs[0].data.len;
//     var i: usize = 0;
//     while (i < len) : (i += 1) {
//         var f_bin: T = vecs[vecs.len - 1].binPc(i);
//         comptime var i_vec = vecs.len - 2;
//         inline while (i_vec >= 0) : (i_vec -= 1) {
//             const vec = vecs[i_vec];
//             f_bin = f_bin * vec.nBins + vec.binPc(i);
//         }
//         if (!math.isNan(f_bin)) {
//             const i_bin = @floatToInt(usize, f_bin);
//             bins[i_bin] += 1;
//         }
//     }
// }

// bins[i_bin] += all_in_bounds
// assignments[i] = i_bin + !all_in_bounds * math.maxInt(usize);

// assignments[i] = blk: {
//     const z = zs.data[i];
//     if (!zs.inBounds(z)) break :blk none;
//     var i_bin = zs.bin(z);

//     const y = ys.data[i];
//     if (!ys.inBounds(y)) break :blk none;
//     i_bin = i_bin * yBins + ys.bin(y);

//     const x = xs.data[i];
//     if (!xs.inBounds(x)) break :blk none;
//     i_bin = i_bin * xBins + xs.bin(x);

//     bins[i_bin] += w[b * i];
//     break :blk i_bin;
// };

// const y = ys.data[i];
// const z = zs.data[i];
// if ((xs.inBounds(x) and ys.inBounds(y) and zs.inBounds(z))) {
//     const i_bin = xyBins * zs.bin(z) + xBins * ys.bin(y) + xs.bin(x);
//     assignments[i] = i_bin;
//     bins[i_bin] += if (b) w[i] else w[0];
// }

// pub fn bin3d(bins: []f64, assignments: []usize, maybe_weights: ?[]const f64, xs: Vec, ys: Vec, zs: Vec) void {
//     const weights = if (maybe_weights) |weights| weights else &[1]f64{1};
//     const i_weight: usize = if (maybe_weights) |_| 1 else 0;
//     const xBins = xs.nBins;
//     const yBins = ys.nBins;
//     const xyBins = xBins * yBins;
//     for (xs.data) |x, i| {
//         const y = ys.data[i];
//         const z = zs.data[i];
//         if ((xs.inBounds(x) and ys.inBounds(y) and zs.inBounds(z))) {
//             const index = xyBins * zs.bin(z) + xBins * ys.bin(y) + xs.bin(x);
//             bins[index] += weights[i_weight * index];
//             assignments[i] = index;
//         } else {
//             assignments[i] = math.maxInt(usize);
//         }
//     }
// }

// https://github.com/ziglang/zig/issues/7364

// original formulation in bin3d:
// const index = xyBins * zs.bin(z) + xBins * ys.bin(y) + xs.bin(x);
// = xBins * (yBins * (zIndex) + yIndex) + xIndex
// = xBins*yBins*zIndex + xBins*yIndex + xIndex
//
// our approach:
//   index = 0;
//   index = index * zBins + zIndex;
//   index = index * yBins + yIndex;
//   index = index * xBins + xIndex;
// which results in...
//         = ((0 * zBins + zIndex) * yBins + yIndex) * xBins + xIndex;
//         = zIndex * yBins * xBins + yIndex * xBins + xIndex;
// as desired.
//
// for example:
//   i_bin = 0;
//   i_bin = i_bin * zs.nBins + zs.bin(val);
//   i_bin = i_bin * ys.nBins + ys.bin(val);
//   i_bin = i_bin * xs.nBins + xs.bin(val);
// is equivalent to:
//   i_bin = xyBins * zs.bin(z) + xBins * ys.bin(y) + xs.bin(x);

// pub fn count3d(bins: []u32, xs: Vec, ys: Vec, zs: Vec) void {
//     const xBins = xs.nBinsF64;
//     const yBins = ys.nBinsF64;
//     const len = xs.data.len;
//     var i: usize = 0;
//     while (i < len) : (i += 1) {
//         var f_bin = zs.binPc(i);
//         f_bin = f_bin * yBins + ys.binPc(i);
//         f_bin = f_bin * xBins + xs.binPc(i);
//         if (!math.isNan(f_bin)) {
//             const i_bin = @floatToInt(usize, f_bin);
//             bins[i_bin] += 1;
//         }
//     }
// }

// pub fn bin(bins: []f64, assignments: []usize, maybe_weights: ?[]const f64, vecs: anytype) void {
//     assert(vecs.len > 0);
//     // if the weight vector has length 1, broadcast it.
//     const weights = if (maybe_weights) |weights| weights else &[1]f64{1};
//     const i_weight: usize = if (maybe_weights) |_| 1 else 0;
//     const n_data = vecs[0].data.len;
//     var i_data: usize = 0;
//     while (i_data < n_data) : (i_data += 1) {
//         // iterate over dimensions in reverse to compute i_bin.
//         var i_bin: usize = 0;
//         var all_in_bounds = true;
//         comptime var i_vec = vecs.len - 1;
//         inline while (i_vec >= 0) : (i_vec -= 1) {
//             const vec = vecs[i_vec];
//             const val = vec.data[i_data];
//             all_in_bounds = all_in_bounds and vec.inBounds(val);
//             if (all_in_bounds) i_bin = i_bin * vec.nBins + vec.bin(val);
//         }
//         if (all_in_bounds) {
//             bins[i_bin] += weights[i_weight * i_data];
//             assignments[i_data] = i_bin;
//         } else {
//             assignments[i_data] = math.maxInt(usize);
//         }
//     }
// }

// fn count2dSimd(bins: []u32, xs: anytype, ys: anytype) void {
//     // not faster than non-simd with n=4, but n=8 is a bit faster
//     const n = 8;
//     const xBins = @splat(n, xs.nBins);
//     const len = xs.data.len;
//     var i: usize = 0;
//     const x_lo = @splat(n, xs.trueMin);
//     const x_hi = @splat(n, xs.trueMax);
//     const x_min = @splat(n, xs.min);
//     const x_scale = @splat(n, xs.scale);

//     const y_lo = @splat(n, ys.trueMin);
//     const y_hi = @splat(n, ys.trueMax);
//     const y_min = @splat(n, ys.min);
//     const y_scale = @splat(n, ys.scale);

//     while (i < len) : (i += n) {
//         const x = @Vector(n, f32){
//             xs.data[i],
//             xs.data[i + 1],
//             xs.data[i + 2],
//             xs.data[i + 3],
//             xs.data[i + 4],
//             xs.data[i + 5],
//             xs.data[i + 6],
//             xs.data[i + 7],
//         };
//         const x_in_bounds_lo = x_lo <= x;
//         const x_in_bounds_hi = x <= x_hi;
//         const x_index = @floor((x - x_min) * x_scale);

//         const y = @Vector(n, f32){
//             ys.data[i],
//             ys.data[i + 1],
//             ys.data[i + 2],
//             ys.data[i + 3],
//             ys.data[i + 4],
//             ys.data[i + 5],
//             ys.data[i + 6],
//             ys.data[i + 7],
//         };
//         const y_in_bounds_lo = y_lo <= y;
//         const y_in_bounds_hi = y <= y_hi;
//         const y_index = @floor((y - y_min) * y_scale);

//         const f_index = y_index * xBins + x_index;

//         var j: usize = 0;
//         while (j < n) : (j += 1) {
//             if (x_in_bounds_lo[j] and x_in_bounds_hi[j] and y_in_bounds_lo[j] and y_in_bounds_hi[j]) {
//                 const i_bin = @floatToInt(usize, f_index[j]);
//                 bins[i_bin] += 1;
//             }
//         }
//     }
// }
