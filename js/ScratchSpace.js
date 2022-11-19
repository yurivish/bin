// A convenient aspect of this design is that the previous buffer
// subarrays can continue to be used if we are resized during an alloc.
export class ScratchSpace {
  constructor(initialLength = 1024) {
    this.buf = new Uint32Array(initialLength);
    this.index = 0;
  }
  alloc(length) {
    if (this.index + length > this.buf.length) {
      // Ensure we at least double every time we're resized
      this.buf = new Uint32Array(Math.max(length, 2 * this.buf.length));
      this.index = 0;
    }
    const sub = this.buf.subarray(this.index, this.index + length);
    this.index += length;
    return sub;
  }
  reset() {
    this.index = 0;
  }
}