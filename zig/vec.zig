const std = @import("std");
const math = std.math;

// Vec stores a single dimension of data together with the information required
// to efficiently compute bin indices for data values.
// The data must be non-empty, since otherwse min and max can't be computed.
const Vec = @This();

// Stores the data from which we compute any NaN elements in the domain
data: []const f64,

// Stores the endpoints of the bin domain. Note that min and max are not
// necessarily the min and max of the data, since it can be specified manually.
// In particular, the domain extent may be reversed, with min greater than max.
min: f64,
max: f64,

// Scale multiplier to rescale (value - min) to [0, nbins - eps].
scale: f64,
nBins: usize,
nBinsF64: f64,

// Store the true min and max for faster bounds checks
trueMin: f64,
trueMax: f64,

pub noinline fn initAuto(nBins: usize, data: []const f64, min: f64, max: f64) !Vec {
    const nanMin = math.isNan(min);
    const nanMax = math.isNan(max);
    if (nanMin and nanMax) return Vec.init(nBins, data);
    if (nanMin) return Vec.initMax(nBins, data, max);
    if (nanMax) return Vec.initMin(nBins, data, min);
    return Vec.initMinMax(nBins, data, min, max);
}

pub fn init(nBins: usize, data: []const f64) !Vec {
    var min = math.inf(f64);
    var max = -math.inf(f64);
    for (data[0..]) |d| {
        min = math.min(d, min);
        max = math.max(d, max);
    }
    return Vec.initMinMax(nBins, data, min, max);
}

pub fn initMin(nBins: usize, data: []const f64, min: f64) !Vec {
    var max = -math.inf(f64);
    for (data) |d| max = math.max(d, max);
    return Vec.initMinMax(nBins, data, min, max);
}

pub fn initMax(nBins: usize, data: []const f64, max: f64) !Vec {
    var min = math.inf(f64);
    for (data) |d| min = math.min(d, min);
    return Vec.initMinMax(nBins, data, min, max);
}

pub noinline fn initMinMax(nBins: usize, data: []const f64, min: f64, max: f64) !Vec {
    if (!(math.isFinite(min) and math.isFinite(max))) return error.NonFiniteMinOrMax;

    // This eps value can be subtracted from an f64 < 32768^2 to decrease it, and
    // can be added to an f64 < prevfloat(32768^2) to increase it. Used to slightly
    // perturb bin positions downwards so that binIndex assigns the correct bin
    // to the maximum value. Note: Needs testing; not fully confident this works.
    const eps64 = 1e-7;
    // can be subtracted from f32 (2^12 = 4096) to decrease it, so we can have up to 4k bins in f32 mode?
    // const eps32 = 0.000244140625
    // Scale factor from [min, max] to [0..bins.len - eps]. The epsilon ensures
    // that maximum point gets mapped into the largest bin rather than beyond it.
    // It does squish all of the projected values slightly, but I haven't found this to be an issue for my use cases.
    // If max is equal to min, we want to map all value to the zeroth bin.
    const denom = if (max - min == 0) math.inf(f64) else max - min;
    const scale = (@intToFloat(f64, nBins) - eps64) / denom;
    return Vec{
        .data = data,
        .min = min,
        .max = max,
        .scale = scale,
        .nBins = nBins,
        .nBinsF64 = @intToFloat(f64, nBins),
        .trueMin = math.min(min, max),
        .trueMax = math.max(min, max),
    };
}

pub fn bin(vec: Vec, value: f64) usize {
    // Compute the bin index for the given value.
    const pc = (value - vec.min) * vec.scale;
    return @floatToInt(usize, pc);
}

pub fn binPc(vec: Vec, index: usize) f64 {
    const value = vec.data[index];
    if (!vec.inBounds(value)) return math.nan(f64);
    return @floor((value - vec.min) * vec.scale);
}

pub fn inBounds(vec: Vec, value: f64) bool {
    return vec.trueMin <= value and value <= vec.trueMax;
}

const Error = error{NonFiniteMinOrMax};

pub fn errString(e: Vec.Error) f64 {
    const Strings = struct {
        const NonFiniteMinOrMax = "min or max is non-finite (ie., infinity or NaN)";
    };
    const s = switch (e) {
        Vec.Error.NonFiniteMinOrMax => Strings.NonFiniteMinOrMax,
    };
    // return pointer to static string.
    // note: assumes that we're on a 32-bit platform.
    return ptrLenF64(@intCast(u32, @ptrToInt(&s[0])), @intCast(u32, s.len));
}

fn ptrLenF64(ptr: u32, len: u32) f64 {
    // Encode a pair of a u32 pointer and length as a single f64.
    // This is a handy way to pass back (ptr, len) to WASM, which
    // natively supports only basic number types right now.
    return @bitCast(f64, [2]u32{ ptr, len });
}

const expect = @import("std").testing.expect;
const expectEqual = @import("std").testing.expectEqual;
const gotWant = std.testing.gotWant;

test "Vec.initMinMax" {
    var data = [_]f64{ 5, 1, 10 };
    const v = try Vec.initMinMax(3, data[0..], 3, 7);
    try expectEqual(v.min, 3);
    try expectEqual(v.max, 7);
}

test "Vec.initMin" {
    var data = [_]f64{ 5, 1, 10 };
    const v = try Vec.initMin(3, data[0..], 3);
    try expectEqual(v.min, 3);
    try expectEqual(v.max, 10);
}

test "Vec.initMax" {
    var data = [_]f64{ 5, 1, 10 };
    const v = try Vec.initMax(3, data[0..], 7);
    try expectEqual(v.min, 1);
    try expectEqual(v.max, 7);
}

test "Vec.init" {
    var data = [_]f64{ 5, 1, 10 };
    const v = try Vec.init(3, data[0..]);
    try expectEqual(v.min, 1);
    try expectEqual(v.max, 10);
}

test "Vec.initAuto" {
    var data = [_]f64{ 5, 1, 10 };
    {
        const v = try Vec.initAuto(3, data[0..], math.nan(f64), math.nan(f64));
        try expectEqual(v.min, 1);
        try expectEqual(v.max, 10);
    }
    {
        const v = try Vec.initAuto(3, data[0..], 3, math.nan(f64));
        try expectEqual(v.min, 3);
        try expectEqual(v.max, 10);
    }
    {
        const v = try Vec.initAuto(3, data[0..], math.nan(f64), 7);
        try expectEqual(v.min, 1);
        try expectEqual(v.max, 7);
    }
}

// todo:
// - test data with inifities and only nans; these cases should throw an error.
// - test cases where max > min, now that we support it