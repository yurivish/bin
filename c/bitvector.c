// https://observablehq.com/d/332a94f6f4b649cf
// https://gcc.gnu.org/onlinedocs/gcc/Other-Builtins.html
#include <stdint.h>

// walloc
typedef __SIZE_TYPE__ size_t;
#define WASM_EXPORT(name) __attribute__((export_name(#name))) name
void *malloc(size_t size);
void free(void *p);
void* WASM_EXPORT(walloc)(size_t size) { return malloc(size); }
void WASM_EXPORT(wfree)(void* ptr) { free(ptr); }

// common
#define NULL 0
#define ERROR -1

// bitvector blocks
#define BLOCK_TYPE uint64_t
#define BLOCK_SIZE (8 * sizeof(BLOCK_TYPE))
#define BLOCK_POPCOUNT(x) (BLOCK_SIZE == 32 ?  __builtin_popcount(x) : __builtin_popcountll(x))


// ideas:
// [] interleave 32-bit count with 32-bit rank, retrieve as a i64?
// [] implement a batched rank1 function that makes two rank1 calls in order to
//    reduce the number of cross-wasm calls


// external functions
// extern int foo(int arg1, int arg2);

typedef struct {
  BLOCK_TYPE* blocks;
  BLOCK_TYPE* rankSuperblocks;
  int numOnes;
  int maxOneIndex;
  int length;
  int storedLength;
} bitvector;

int ceilDivide(int a, int b) {
  // return ceil(a / b), avoiding overflow
  if (a == 0) return 0;
  return 1 + (a - 1) / b;
}

bitvector* WASM_EXPORT(bitvector_init)(int length) {
  bitvector* v = malloc(sizeof(bitvector));
  if (v == NULL) return NULL;
  
  // allocate blocks (even if length == 0)
  int numBlocks = ceilDivide(length, BLOCK_SIZE);
  BLOCK_TYPE *blocks = malloc(numBlocks * sizeof(BLOCK_TYPE));
  if (blocks == NULL) {
    free(v);
    return NULL;
  }
  for (int i = 0; i < numBlocks; i++) blocks[i] = 0;
  v->blocks = blocks;

  v->rankSuperblocks = NULL;
  v->numOnes = 0;
  v->maxOneIndex = -1;
  v->length = length;
  v->storedLength = 0;
  return v;
}

int WASM_EXPORT(bitvector_one)(bitvector* v, int i) {
  if (i >= v->length) return ERROR; // i must be < length
  int blockIndex = i / BLOCK_SIZE; 
  int bitOffset = i % BLOCK_SIZE; // i & (BLOCK_SIZE - 1);
  v->blocks[blockIndex] |= 1ULL << bitOffset;
  v->numOnes += 1;
  if (i > v->maxOneIndex) v->maxOneIndex = i;
  return 0;
}

int WASM_EXPORT(bitvector_finish)(bitvector* v) {
  v->storedLength = v->maxOneIndex + 1;
  // reallocate blocks if we oversized initially
  int origNumBlocks = ceilDivide(v->length, BLOCK_SIZE); 
  int numBlocks = ceilDivide(v->storedLength, BLOCK_SIZE); 
  if (numBlocks < origNumBlocks) {
    BLOCK_TYPE *origBlocks = v->blocks;
    BLOCK_TYPE *blocks = malloc(numBlocks * sizeof(BLOCK_TYPE));
    if (blocks == NULL) return ERROR;
    for (int i = 0; i < numBlocks; i++) blocks[i] = origBlocks[i];
    free(origBlocks);
    v->blocks = blocks;
  }

  BLOCK_TYPE *rankSuperblocks  = malloc(numBlocks * sizeof(BLOCK_TYPE));
  if (rankSuperblocks == NULL) return ERROR;

  BLOCK_TYPE* blocks = v->blocks;
  rankSuperblocks[0] = BLOCK_POPCOUNT(blocks[0]);
  for (int i = 1; i < numBlocks; i++) {
    rankSuperblocks[i] = rankSuperblocks[i - 1] + BLOCK_POPCOUNT(blocks[i]);
  }
  v->rankSuperblocks = rankSuperblocks;
  return 0;
}

// If we store blocks and rank superblocks interleaved as 
// [block] [superblock] [block] [superblock]
// ^^^^^^^^^^^^^^^^^^^^ 64 bits
// and retrieve 64 bits at a time, the below can become
// (B & 0x0000ffff) - BLOCK_POPCOUNT32((B >> 32) & mask);
// if mask is 32 bits. We can do one shift instead of two
// Need to be careful with right shifts; see Seacord, pp. 71
// Might be better to just write it all in Zig.
int WASM_EXPORT(bitvector_rank1)(bitvector* v, int i) {
  if (i < 0) return 0;
  if (i >= v->storedLength) return v->numOnes;
  int blockIndex = i / BLOCK_SIZE;
  int lowBitIndex = i % BLOCK_SIZE; // i & (BLOCK_SIZE - 1);
  BLOCK_TYPE rankSuperblock = v->rankSuperblocks[blockIndex];
  BLOCK_TYPE block = v->blocks[blockIndex];
  BLOCK_TYPE mask = 0xfffffffffffffffeULL << lowBitIndex;
  return rankSuperblock - BLOCK_POPCOUNT(block & mask);
}

int WASM_EXPORT(bitvector_rank0)(bitvector* v, int i) {
  if (i < 0) return 0;
  // We check against length here (rather than storedLength)
  // since we need to count the implicitly-represented zeros.
  // note: the final block is padded with zeros so rank0 will return
  // incorrect results if called with an out-of-bounds index that is
  // within the final block. So we do the bounds checks here too.
  if (i >= v->length) return v->length - v->numOnes;
  return i - bitvector_rank1(v, i) + 1;
}

// unhinted select (hinted: allow specifying L and R as arguments)
int WASM_EXPORT(bitvector_select1)(bitvector* v, int i) {
  if (i < 1) return ERROR; // out of bounds: i < 1
  if (i > v->numOnes) return ERROR; // out of bounds: i > numOnes
  // Search based on the structure of binarySearchBefore
  int L = 0;
  int R = v->length;
  while (L < R) {
    int m = (L + R) >> 1;
    if (bitvector_rank1(v, m) < i) L = m + 1;
    else R = m;
  }
  return L;
}

int WASM_EXPORT(bitvector_select0)(bitvector* v, int i) {
  if (i < 1) return ERROR; // out of bounds: i < 1
  int numZeros = v->length - v->numOnes;
  if (i > numZeros) return ERROR; // out of bounds: i > numZeros
  // Search based on the structure of binarySearchBefore
  int L = 0;
  int R = v->length;
  while (L < R) {
    int m = (L + R) >> 1;
    if (bitvector_rank0(v, m) < i) L = m + 1;
    else R = m;
  }
  return L;
}

int WASM_EXPORT(bitvector_access)(bitvector* v, int i) {
  if (i < 0 || i >= v->length) return ERROR; // access: out of bounds
  if (i >= v->storedLength) return 0;
  int blockIndex = i / BLOCK_SIZE;
  int bitOffset = i % BLOCK_SIZE; // i & (BLOCK_SIZE - 1);
  BLOCK_TYPE block = v->blocks[blockIndex];
  BLOCK_TYPE targetMask = 1ULL << bitOffset; // mask out the target bit
  return (block & targetMask) != 0;
}

void WASM_EXPORT(bitvector_free)(bitvector* v) {
  free(v->blocks);
  free(v->rankSuperblocks);
  free(v);
}

int WASM_EXPORT(add1)(int i) {
  return i + 1;
}
