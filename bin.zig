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

test "1d facets" {

    // note: number of bins = xBins * numFacets
    var bins_a = [_]f64{0, 0, 0, 0, 0, 0};
    var bins: []f64 = bins_a[0..];

    var assignments_a = [_]usize{0,0,0,0,0,0};
    var assignments: []usize = assignments_a[0..];

    var xs_a = [_]f64{1, 2, 2, 3, 3, 3};
    var xs = try Vec.init(3, xs_a[0..]);

    var facets_a = [_]usize{0, 0, 1, 1, 1, 1};
    var facets: []usize = facets_a[0..];

    binFacets1d(bins, assignments, xs, facets);

   std.debug.print("Hello World!", .{});

}

pub fn binFacets1d(bins: []f64, assignments: []usize, xs: Vec, facets: []usize) void {
    const nBinsPerFacet = xs.nBins;
    for (xs.data) |x, i| {
        if (xs.inBounds(x)) {
            const facetOffset = facets[i] * nBinsPerFacet;
            const index = facetOffset + xs.bin(x);
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

// 6715 bytes
pub fn binFacets2d(bins: []f64, assignments: []usize, xs: Vec, ys: Vec, facets: []usize) void {
    const xBins = xs.nBins;
    const yBins = ys.nBins;
    const nBinsPerFacet = xBins * yBins;
    for (xs.data) |x, i| {
        const y = ys.data[i];
        if ((xs.inBounds(x) and ys.inBounds(y))) {
            const facetOffset = facets[i] * nBinsPerFacet;
            const index = facetOffset + xBins * ys.bin(y) + xs.bin(x);
            bins[index]  += 1;
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
