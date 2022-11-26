// note: top-level await is only for --target = es2022 and above:
// esbuild js/bin.js --bundle --loader:.wasm=binary --target=es2022 --format=esm --outfile=bin.bundle.js --watch
//
// import wasm from './dist/bitvector.wasm';
// const C = await WebAssembly.instantiate(wasm).then((r) => r.instance.exports);

// let C = null;

export class CBitVector {
  constructor(length, opts) {
    const { C } = opts;
    this.C = C;
    this.v = C.bitvector_init(length);
    if (this.v === 0) throw new Error('c_error: bitvector_init');
    this.length = length;
  }
  destroy() {
    const ret = this.C.bitvector_free(this.v);
  }
  one(i) {
    const ret = this.C.bitvector_one(this.v, i | 0);
    if (ret === -1) throw new Error('c_error: bitvector_one');
    return ret;
  }
  finish() {
    const ret = this.C.bitvector_finish(this.v);
    if (ret === -1) throw new Error('c_error: bitvector_finish');
    return this;
  }
  rank1(i) {
    // this.C.add1(i); // testing the effects of wasm call overhead
    const ret = this.C.bitvector_rank1(this.v, i | 0);
    return ret;
  }
  rank0(i) {
    const ret = this.C.bitvector_rank0(this.v, i | 0);
    return ret;
  }
  select1(i) {
    const ret = this.C.bitvector_select1(this.v, i | 0);
    return ret;
  }
  select0(i) {
    const ret = this.C.bitvector_select0(this.v, i | 0);
    return ret;
  }
  access(i) {
    const ret = this.C.bitvector_access(this.v, i | 0);
    if (ret == -1) throw new Error('c_error: bitvector_access');
    return ret;
  }
  approxSizeInBits() {
    return -1; // for now
  }
}
