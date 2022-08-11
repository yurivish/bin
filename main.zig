const std = @import("std");
const math = std.math;
const bin = @import("bin.zig");
const edt = @import("edt.zig");
const Vec = @import("vec.zig");
const Rgba = bin.Rgba;

test "debug" {}

// todo: apply const to readonly arrays?

export fn bin1d(
    // zig fmt: off
    bins: [*]f64, 
    assignments: [*]usize,
    xBins: usize, xs: [*]const f64, xDomain: *[2]f64,
    weights: ?[*]const f64, len: usize,
    // zig fmt: on
) f64 {
    if (len == 0) return 0;
    const xVec = Vec.initAuto(xBins, xs[0..len], xDomain[0], xDomain[1]) catch |e| return Vec.errString(e);
    xDomain.* = .{ xVec.min, xVec.max };
    const ws = if (weights) |ws| ws[0..len] else null;
    bin.bin(bins[0..xBins], assignments[0..len], ws, .{xVec});
    return 0;
}

export fn bin2d(
    // zig fmt: off
    binsPtr: [*]f64,
    assignmentsPtr: [*]usize,
    xBins: usize, xs: [*]const f64, xDomain: *[2]f64,
    yBins: usize, ys: [*]const f64, yDomain: *[2]f64,
    weightsPtr: ?[*]const f64, len: usize,
    binsEdtPtr: [*]f64, v: [*]usize, f: [*]f64, z: [*]f64,
    outlierPtr: [*]u8, outlierCutoff: f64, outlierRadius: f64, 
    // zig fmt: on
) f64 {
    if (len == 0) return 0;
    const bins = binsPtr[0 .. xBins * yBins];
    const assignments = assignmentsPtr[0..len];
    const outlier = outlierPtr[0..len];
    const xVec = Vec.initAuto(xBins, xs[0..len], xDomain[0], xDomain[1]) catch |e| return Vec.errString(e);
    const yVec = Vec.initAuto(yBins, ys[0..len], yDomain[0], yDomain[1]) catch |e| return Vec.errString(e);
    xDomain.* = .{ xVec.min, xVec.max };
    yDomain.* = .{ yVec.min, yVec.max };
    const ws = if (weightsPtr) |ws| ws[0..len] else null;
    bin.bin(bins, assignments, ws, .{ xVec, yVec });

    if (outlierCutoff > 0) {
        const binsEdt = binsEdtPtr[0 .. xBins * yBins];
        for (bins) |b, i| binsEdt[i] = @intToFloat(f64, @boolToInt(b <= outlierCutoff));
        const max = math.max(xBins, yBins);
        edt.edtBinary2d(binsEdt, binsEdt, xBins, yBins, v[0..max], f[0..max], z[0 .. max + 1]);
        const boolEdt = binsEdt; // @ptrCast([*]u8, binsEdtPtr)[0 .. xBins * yBins];
        for (binsEdt) |b, i| {
            boolEdt[i] = if (b > outlierRadius) 1 else 0;
        }
        for (assignments) |a, i| {
            outlier[i] = if (a != math.maxInt(usize) and boolEdt[a] == 1) 1 else 0;
        }
    }
    return 0;
}

export fn bin3d(
    // zig fmt: off
    binsPtr: [*]f64,
    assignmentsPtr: [*]usize,
    xBins: usize, xs: [*]f64, xDomain: *[2]f64,
    yBins: usize, ys: [*]f64, yDomain: *[2]f64,
    zBins: usize, zs: [*]f64, zDomain: *[2]f64,
    weights: ?[*]const f64, len: usize,
    binsEdtPtr: [*]f64, vPtr: [*]usize, fPtr: [*]f64, zPtr: [*]f64,
    outlierPtr: [*]u8, outlierCutoff: f64, outlierRadius: f64, 
    // zig fmt: on
) f64 {
    if (len == 0) return 0;
    const nBins = xBins * yBins * zBins;
    const bins = binsPtr[0..nBins];
    const assignments = assignmentsPtr[0..len];
    const outlier = outlierPtr[0..len];
    const xVec = Vec.initAuto(xBins, xs[0..len], xDomain[0], xDomain[1]) catch |e| return Vec.errString(e);
    const yVec = Vec.initAuto(yBins, ys[0..len], yDomain[0], yDomain[1]) catch |e| return Vec.errString(e);
    const zVec = Vec.initAuto(zBins, zs[0..len], zDomain[0], zDomain[1]) catch |e| return Vec.errString(e);
    xDomain.* = .{ xVec.min, xVec.max };
    yDomain.* = .{ yVec.min, yVec.max };
    zDomain.* = .{ zVec.min, zVec.max };
    const ws = if (weights) |ws| ws[0..len] else null;
    bin.bin(bins, assignments[0..len], ws, .{ xVec, yVec, zVec });
    if (outlierCutoff > 0) {
        var lo: usize = 0;
        const xy = xBins * yBins;
        const max = math.max(xBins, yBins);
        const v = vPtr[0..max];
        const f = fPtr[0..max];
        const z = zPtr[0 .. max + 1];
        while (lo < bins.len) : (lo += xy) {
            const binsEdt = binsEdtPtr[lo .. lo + xy];
            const xyBins = bins[lo .. lo + xy];
            for (xyBins) |b, i| binsEdt[i] = @intToFloat(f64, @boolToInt(b <= outlierCutoff));
            edt.edtBinary2d(binsEdt, binsEdt, xBins, yBins, v, f, z);
        }
        // can do this across all 2d slices
        const binsEdt = binsEdtPtr[0..nBins];
        const boolEdt = binsEdt; // @ptrCast([*]u8, binsEdtPtr)[0..nBins];
        for (binsEdt) |b, i| {
            boolEdt[i] = if (b > outlierRadius) 1 else 0;
        }
        for (assignments) |a, i| {
            outlier[i] = if (a != math.maxInt(usize) and boolEdt[a] == 1) 1 else 0;
        }
    }
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
