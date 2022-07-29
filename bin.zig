const std = @import("std");
const math = std.math;
const Vec = @import("vec.zig");
const assert = std.debug.assert;

pub const Rgba = packed struct { r: u8, g: u8, b: u8, a: u8 };

test {
    var bins = [_]f64{ 0, 0, 0 };
    var a = [_]f64{ 1, 2, 3 }; // , math.inf(f64)}; // also test empty array; should crash.
    var vec = Vec.initAuto(bins.len, a[0..], math.nan(f64), math.nan(f64));
    bin1d(bins[0..], vec);
}

pub fn bin1d(bins: []f64, xs: Vec, binAssignments: []usize) void {
    for (xs.data) |x, i| {
        if (xs.inBounds(x)) {
            const index = xs.bin(x);
            bins[index] += 1;
            binAssignments[i] = index;
        } else {
            binAssignments[i] = math.maxInt(usize);
        }
    }
}

pub fn bin2d(bins: []f64, xs: Vec, ys: Vec, binAssignments: []usize) void {
    const xBins = xs.nBins;
    // The (0, 0) bin is bottom left
    for (xs.data) |x, i| {
        const y = ys.data[i];
        if ((xs.inBounds(x) and ys.inBounds(y))) {
            const index = xBins * ys.bin(y) + xs.bin(x);
            bins[index] += 1;
            binAssignments[i] = index;
        } else {
            binAssignments[i] = math.maxInt(usize);
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
