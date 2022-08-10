const std = @import("std");
const math = std.math;
const Vec = @import("vec.zig");
const assert = std.debug.assert;

pub const Rgba = packed struct { r: u8, g: u8, b: u8, a: u8 };

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

// pub fn bin3d(bins: []f64, assignments: []usize, xs: Vec, ys: Vec, zs: Vec) void {
//     const xBins = xs.nBins;
//     const yBins = ys.nBins;
//     const xyBins = xBins * yBins;
//     for (xs.data) |x, i| {
//         const y = ys.data[i];
//         const z = zs.data[i];
//         if ((xs.inBounds(x) and ys.inBounds(y) and zs.inBounds(z))) {
//             const index = xyBins * zs.bin(z) + xBins * ys.bin(y) + xs.bin(x);
//             bins[index] += 1;
//             assignments[i] = index;
//         } else {
//             assignments[i] = math.maxInt(usize);
//         }
//     }
// }

pub fn bin1d(bins: []f64, assignments: []usize, xs: Vec) void {
    bin(bins, assignments, .{xs});
}

pub fn bin2d(bins: []f64, assignments: []usize, xs: Vec, ys: Vec) void {
    bin(bins, assignments, .{ xs, ys });
}

pub fn bin3d(bins: []f64, assignments: []usize, xs: Vec, ys: Vec, zs: Vec) void {
    bin(bins, assignments, .{ xs, ys, zs });
}

pub fn bin(bins: []f64, assignments: []usize, vecs: anytype) void {
    assert(vecs.len > 0);
    const n_data = vecs[0].data.len;
    var i_data: usize = 0;
    while (i_data < n_data) : (i_data += 1) {
        var i_bin: usize = 0;
        var all_in_bounds = true;
        // to compute i_bin, iterate over dimensions in reverse
        comptime var i_vec = vecs.len - 1;
        inline while (i_vec >= 0) : (i_vec -= 1) {
            const vec = vecs[i_vec];
            const val = vec.data[i_data];
            all_in_bounds = all_in_bounds and vec.inBounds(val);
            if (all_in_bounds) i_bin = i_bin * vec.nBins + vec.bin(val);
        }
        if (all_in_bounds) {
            bins[i_bin] += 1;
            assignments[i_data] = i_bin;
        } else {
            assignments[i_data] = math.maxInt(usize);
        }
    }
}

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

pub fn colorize(colors: []Rgba, vs: Vec, ramp: []Rgba) void {
    // Compute a color for each bin by looking up the appropriate ramp value
    const unknown = Rgba{ .r = 255, .g = 0, .b = 255, .a = 255 };
    for (vs.data) |v, i| {
        colors[i] = if (vs.inBounds(v)) ramp[vs.bin(v)] else unknown;
    }
}
