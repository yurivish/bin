const std = @import("std");
const math = std.math;
const assert = std.debug.assert;

pub fn NextVec(comptime T: type) type {
    assert(T == f32 or T == f64);

    // todo: rename trueMin, trueMax to
    return struct {
        data: []const T,
        // rename to extentLo, extentHi? loBound, hiBound?
        min: T,
        max: T,
        scale: T,
        nBins: T,

        // note *still* not the tru min or max; just the min and max of hte range;
        // rename to extentMin and extentMax, or minBound, maxBound?
        trueMin: T,
        trueMax: T,

        const Vec = @This();

        pub noinline fn initAuto(nBins: usize, data: []const T, min: T, max: T) !Vec {
            const nanMin = math.isNan(min);
            const nanMax = math.isNan(max);
            if (nanMin and nanMax) return Vec.init(nBins, data);
            if (nanMin) return Vec.initMax(nBins, data, max);
            if (nanMax) return Vec.initMin(nBins, data, min);
            return Vec.initMinMax(nBins, data, min, max);
        }

        pub fn init(nBins: usize, data: []const T) !Vec {
            var min = math.inf(T);
            var max = -math.inf(T);
            for (data[0..]) |d| {
                min = math.min(d, min);
                max = math.max(d, max);
            }
            return Vec.initMinMax(nBins, data, min, max);
        }

        pub fn initMin(nBins: usize, data: []const T, min: T) !Vec {
            var max = -math.inf(T);
            for (data) |d| max = math.max(d, max);
            return Vec.initMinMax(nBins, data, min, max);
        }

        pub fn initMax(nBins: usize, data: []const T, max: T) !Vec {
            var min = math.inf(T);
            for (data) |d| min = math.min(d, min);
            return Vec.initMinMax(nBins, data, min, max);
        }

        pub noinline fn initMinMax(nBins: usize, data: []const T, min: T, max: T) !Vec {
            if (!(math.isFinite(min) and math.isFinite(max))) return error.NonFiniteMinOrMax;

            const eps = switch (T) {
                f64 => blk: {
                    // f64 can exactly represent inegers up to  = 2^53 = 9,007,199,254,740,993
                    assert(nBins <= 32768 * 32768);
                    break :blk 1e-7;
                },
                f32 => blk: {
                    // f32 can exactly represent integers up to and including 4096^2 = 2^24 = 16,777,216
                    // !! so we can use f32 to compute intermediate array indices in 2d, but not 3d.
                    assert(nBins <= 4096);
                    break :blk 0.000244140625;
                },
                else => unreachable,
            };

            // can be subtracted from f32 (2^12 = 4096) to decrease it, so we can have up to 4k bins in f32 mode?
            // const eps32 =
            // Scale factor from [min, max] to [0..bins.len - eps]. The epsilon ensures
            // that maximum point gets mapped into the largest bin rather than beyond it.
            // It does squish all of the projected values slightly, but I haven't found this to be an issue for my use cases.
            // If max is equal to min, we want to map all value to the zeroth bin.
            const denom = if (max - min == 0) math.inf(T) else max - min;
            const scale = (@intToFloat(T, nBins) - eps) / denom;
            return Vec{
                .data = data,
                .min = min,
                .max = max,
                .scale = scale,
                .nBins = @intToFloat(T, nBins),
                .trueMin = math.min(min, max),
                .trueMax = math.max(min, max),
            };
        }

        pub fn bin(vec: Vec, value: T) usize {
            // Compute the bin index for the given value.
            const pc = (value - vec.min) * vec.scale;
            return @floatToInt(usize, pc);
        }

        pub fn binPc(vec: Vec, index: usize) T {
            const value = vec.data[index];
            if (!vec.inBounds(value)) return math.nan(T);
            return @floor((value - vec.min) * vec.scale);
        }

        const foo = struct { inBounds: bool, index: T };
        pub fn binPcBounds(vec: Vec, index: usize) foo {
            const value = vec.data[index];
            return .{
                .inBounds = vec.inBounds(value),
                .index = @floor((value - vec.min) * vec.scale),
            };
        }

        pub fn inBounds(vec: Vec, value: T) bool {
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
    };
}

fn ptrLenF64(ptr: u32, len: u32) f64 {
    // Encode a pair of a u32 pointer and length as a single T.
    // This is a handy way to pass back (ptr, len) to WASM, which
    // natively supports only basic number types right now.
    return @bitCast(f64, [2]u32{ ptr, len });
}
