// A convenient aspect of this design is that the previous buffer
// subarrays can continue to be used if we are resized during an alloc.
// note: could have a resetAndAlloc function like we did with the wasm allocator;
// that would match better to a C implementation since we could free the previous
// alloc, and allocate all currently-needed buffers from the new allocation, whereas
// the current design relies on the GC to keep the previous buffers around.
export class ScratchSpace {
  constructor(initialLength = 1024) {
    this.buf = new Uint32Array(initialLength);
    this.index = 0;
  }
  allocU32(length) {
    if (this.index + length > this.buf.length) {
      // Ensure we at least double every time we're resized
      this.buf = new Uint32Array(Math.max(length, 2 * this.buf.length));
      this.index = 0;
    }
    const subarray = this.buf.subarray(this.index, this.index + length);
    this.index += length;
    return subarray;
  }
  allocF64(length) {
    // quick hack version; we could instead have this.buf be a Buffer
    // and not need to allocate the intermediate Uint32Array object
    // here, using an allocBytes method that uses/resizes the buffer.
    const u32 = this.allocU32(2 * length)
    return new Float64Array(u32.buffer, 0, length)

  }
  reset() {
    this.index = 0;
  }
}