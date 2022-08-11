const std = @import("std");
const math = std.math;
const bin = @import("bin.zig");
const edt = @import("edt.zig");
const Vec = @import("vec.zig");
const Rgba = bin.Rgba;

test "debug" {}

fn possiblyWeightedBin(bins: []f64, assignments: []usize, vecs: anytype, weights: ?[]f64) void {
    const ws = if (weights) |ws| ws else &[1]f64{1};
    bin.bin(bins, assignments, ws, vecs);
}

export fn bin1d(
    // zig fmt: off
    bins: [*]f64, 
    assignments: [*]usize,
    xBins: usize, xs: [*]f64, xDomain: *[2]f64,
    weights: ?[*]f64, len: usize,
    // zig fmt: on
) f64 {
    if (len == 0) return 0;
    const xVec = Vec.initAuto(xBins, xs[0..len], xDomain[0], xDomain[1]) catch |e| return Vec.errString(e);
    xDomain.* = .{ xVec.min, xVec.max };
    const ws = if (weights) |ws| ws[0..len] else null;
    possiblyWeightedBin(bins[0..xBins], assignments[0..len], .{xVec}, ws);
    return 0;
}

export fn bin2d(
    // zig fmt: off
    bins: [*]f64,
    assignments: [*]usize,
    xBins: usize, xs: [*]f64, xDomain: *[2]f64,
    yBins: usize, ys: [*]f64, yDomain: *[2]f64,
    weights: ?[*]f64, len: usize,
    // zig fmt: on
) f64 {
    if (len == 0) return 0;
    const xVec = Vec.initAuto(xBins, xs[0..len], xDomain[0], xDomain[1]) catch |e| return Vec.errString(e);
    const yVec = Vec.initAuto(yBins, ys[0..len], yDomain[0], yDomain[1]) catch |e| return Vec.errString(e);
    xDomain.* = .{ xVec.min, xVec.max };
    yDomain.* = .{ yVec.min, yVec.max };
    const ws = if (weights) |ws| ws[0..len] else null;
    possiblyWeightedBin(bins[0 .. xBins * yBins], assignments[0..len], .{ xVec, yVec }, ws);
    return 0;
}

export fn bin3d(
    // zig fmt: off
    bins: [*]f64,
    assignments: [*]usize,
    xBins: usize, xs: [*]f64, xDomain: *[2]f64,
    yBins: usize, ys: [*]f64, yDomain: *[2]f64,
    zBins: usize, zs: [*]f64, zDomain: *[2]f64,
    weights: ?[*]f64, len: usize,
    // zig fmt: on
) f64 {
    if (len == 0) return 0;
    const xVec = Vec.initAuto(xBins, xs[0..len], xDomain[0], xDomain[1]) catch |e| return Vec.errString(e);
    const yVec = Vec.initAuto(yBins, ys[0..len], yDomain[0], yDomain[1]) catch |e| return Vec.errString(e);
    const zVec = Vec.initAuto(zBins, zs[0..len], zDomain[0], zDomain[1]) catch |e| return Vec.errString(e);
    xDomain.* = .{ xVec.min, xVec.max };
    yDomain.* = .{ yVec.min, yVec.max };
    zDomain.* = .{ zVec.min, zVec.max };
    const ws = if (weights) |ws| ws[0..len] else null;
    possiblyWeightedBin(bins[0 .. xBins * yBins * zBins], assignments[0..len], .{ xVec, yVec, zVec }, ws);
    return 0;
}

export fn colorize(
    // zig fmt: off
    colors: [*]Rgba, 
    xs: [*]f64, xDomain: *[2]f64, len: usize, 
    ramp: [*]Rgba, lenRamp: usize, 
    // zig fmt: on
) f64 {
    if (len == 0) return 0;
    const xVec = Vec.initAuto(lenRamp, xs[0..len], xDomain[0], xDomain[1]) catch |e| return Vec.errString(e);
    xDomain.* = .{ xVec.min, xVec.max };
    bin.colorize(colors[0..len], xVec, ramp[0..lenRamp]);
    return 0;
}

export fn edt2d(
    // zig fmt: off
    out: [*]f64, in: [*]f64,
    width: usize, height: usize, binary: bool,
    v: [*]usize, f: [*]f64, z: [*]f64,
    // zig fmt: on
) void {
    const len = width * height;
    if (len == 0) return;
    const func = if (binary) edt.edtBinary2d else edt.edt2d;
    func(
        out[0..len],
        in[0..len],
        width,
        height,
        v[0..math.max(width, height)],
        f[0..math.max(width, height)],
        z[0 .. math.max(width, height) + 1],
    );
}
