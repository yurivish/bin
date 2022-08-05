const std = @import("std");
const math = std.math;
const Vec = @import("vec.zig");
const assert = std.debug.assert;

pub const Rgba = packed struct { r: u8, g: u8, b: u8, a: u8 };

pub fn bin1d(bins: []f64, assignments: []usize, xs: Vec) void {
    for (xs.data) |x, i| {
        if (xs.inBounds(x)) {
            const index = xs.bin(x);
            bins[index] += 1;
            assignments[i] = index;
        } else {
            assignments[i] = math.maxInt(usize);
        }
    }
}

pub fn bin2d(bins: []f64, assignments: []usize, xs: Vec, ys: Vec) void {
    const xBins = xs.nBins;
    for (xs.data) |x, i| {
        const y = ys.data[i];
        if ((xs.inBounds(x) and ys.inBounds(y))) {
            const index = xBins * ys.bin(y) + xs.bin(x);
            bins[index] += 1;
            assignments[i] = index;
        } else {
            assignments[i] = math.maxInt(usize);
        }
    }
}


pub fn bin3d(bins: []f64, assignments: []usize, xs: Vec, ys: Vec, zs: Vec) void {
    const xBins = xs.nBins;
    const yBins = ys.nBins;
    const xyBins = xBins * yBins;
    for (xs.data) |x, i| {
        const y = ys.data[i];
        const z = zs.data[i];
        if ((xs.inBounds(x) and ys.inBounds(y) and zs.inBounds(z))) {
            const index = xyBins * zs.bin(z) + xBins * ys.bin(y) + xs.bin(x);
            bins[index] += 1;
            assignments[i] = index;
        } else {
            assignments[i] = math.maxInt(usize);
        }
    }
}

pub fn colorize(colors: []Rgba, vs: Vec, ramp: []Rgba) void {
    // Compute a color for each bin by looking up the appropriate ramp value
    const unknown = Rgba{ .r = 255, .g = 0, .b = 255, .a = 255 };
    for (vs.data) |v, i| {
        colors[i] = if (vs.inBounds(v)) ramp[vs.bin(v)] else unknown;
    }
}
