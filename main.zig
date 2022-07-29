const std = @import("std");
const math = std.math;
const bin = @import("bin.zig");
const edt = @import("edt.zig");
const Vec = @import("vec.zig");
const Rgba = bin.Rgba;

// zig fmt: off
export fn bin1d(
    bins: [*]f64, 
    xBins: usize, xs: [*]f64, xDomain: *[2]f64,
    len: usize,
    binAssignments: [*]usize,
) f64 {
    if (len == 0) return 0;
    const xVec = Vec.initAuto(xBins, xs[0..len], xDomain[0], xDomain[1]) catch |e| return Vec.errString(e);
    xDomain.* = .{xVec.min, xVec.max};
    bin.bin1d(bins[0..xBins], xVec, binAssignments[0..len]);
    return 0;
}

export fn bin2d(
    bins: [*]f64,
    xBins: usize, xs: [*]f64, xDomain: *[2]f64,
    yBins: usize, ys: [*]f64, yDomain: *[2]f64,
    len: usize,
    binAssignments: [*]usize,
) f64 {
    if (len == 0) return 0;
    const xVec = Vec.initAuto(xBins, xs[0..len], xDomain[0], xDomain[1]) catch |e| return Vec.errString(e);
    const yVec = Vec.initAuto(yBins, ys[0..len], yDomain[0], yDomain[1]) catch |e| return Vec.errString(e);
    xDomain.* = .{xVec.min, xVec.max};
    yDomain.* = .{yVec.min, yVec.max};
    bin.bin2d(bins[0 .. xBins * yBins], xVec, yVec, binAssignments[0..len]);
    return 0;
}

export fn colorize(
    colors: [*]Rgba, 
    xs: [*]f64, xDomain: *[2]f64, len: usize, 
    ramp: [*]Rgba, lenRamp: usize, 
) f64 {
    if (len == 0) return 0;
    const xVec = Vec.initAuto(lenRamp, xs[0..len], xDomain[0], xDomain[1]) catch |e| return Vec.errString(e);
    xDomain.* = .{xVec.min, xVec.max};
    bin.colorize(colors[0..len], xVec, ramp[0..lenRamp]);
    return 0;
}

export fn edt2d(
    out: [*]f64, in: [*]f64,
    width: usize, height: usize, binary: bool,
    v: [*]usize, f: [*]f64, z: [*]f64,
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
// zig fmt: on
