let c = null

// note: top-level await is only for --target = es2022 and above:
// esbuild bin.js --bundle --target=es2022 --format=esm --outfile=/Users/yurivish/Dropbox/Projects/bin-c/bin.bundle.js --watch
//
// const c = await WebAssembly.instantiateStreaming(
//   fetch("https://compy.local:8000/hello.wasm?" + Math.random()),
//   { env: { foo: console.log } }
// ).then((response) =>response.instance.exports)

export class CBitVector {
  constructor(length) {
    this.v = c.bitvector_init(length);
  }

  one(i) {
    return c.bitvector_one(this.v, i|0);
  }

  finish() {
    return c.bitvector_finish(this.v);
  }

  rank1(i) {
    return c.bitvector_rank1(this.v, i|0);
  }

  rank0(i) {
    return c.bitvector_rank0(this.v, i|0);
  }

  select1(i) {
    return c.bitvector_select1(this.v, i|0);
  }

  select0(i) {
    return c.bitvector_select0(this.v, i|0);
  }

  access(i) {
    return c.bitvector_access(this.v, i|0);
  }

  approxSizeInBits() {
    return -1; // for now
  }
}
