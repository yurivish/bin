// note: top-level await is only for --target = es2022 and above:
// esbuild bin.js --bundle --loader:.wasm=binary --target=es2022 --format=esm --outfile=bin.bundle.js --watch
//
// import wasm from './dist/bitvector.wasm';
// const c = await WebAssembly.instantiate(wasm).then((r) => r.instance.exports);

let c = null;

export class CBitVector {
  constructor(length) {
    this.v = c.bitvector_init(length);
  }
  destroy() {
    c.bitvector_free(this.v);
  }
  one(i) {
    return c.bitvector_one(this.v, i | 0);
  }
  finish() {
    return c.bitvector_finish(this.v);
  }
  rank1(i) {
    // c.add1(i); // testing the effects of wasm call overhead
    return c.bitvector_rank1(this.v, i | 0);
  }
  rank0(i) {
    return c.bitvector_rank0(this.v, i | 0);
  }
  select1(i) {
    return c.bitvector_select1(this.v, i | 0);
  }
  select0(i) {
    return c.bitvector_select0(this.v, i | 0);
  }
  access(i) {
    return c.bitvector_access(this.v, i | 0);
  }
  approxSizeInBits() {
    return -1; // for now
  }
}
