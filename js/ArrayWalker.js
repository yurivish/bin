
// helper for turning recursion into iteration, encapsulating the logic
// of walking an array zero to two outputs for each input.
// left elements are filled in from the front, while right elements are filled
// in from the back, and reset() will preserve all left children in order at
// the front of the array, and move the right children in order immediately
// after them. to interleave left and right in the order they are generated,
// use nextBackIndex() for both the left and right outputs (we cannot use nextFrontIndex()
// for both, since that would overwrite elements as they are being processed
// in left-to-right order).
class ArrayWalker {
  constructor(length, cap) {
    this.length = length; // length taken up by existing elements
    this.cap = cap; // capacity for additional elements
    this.frontIndex = 0;
    this.backIndex = cap; // nextIndex + 1
  }

  // Return the next index from the front (left side) of the array
  // note: ensure that the current value has been retrieved before
  // writing to the left, since it will overwrite the current value.
  nextFrontIndex() {
    const index = this.frontIndex;
    this.frontIndex += 1;
    return index;
  }
  // Return the next index from the back (right side) of the array
  nextBackIndex() {
    // return the next index at which we can append an element
    // (as we fill the array in backwards from arr[cap - 1])
    return (this.backIndex -= 1);
  }

  // `reverse`: Whether to reverse the [backIndex...cap] subarray.
  // In practice, we want to reverse it if we've been iterating
  // the source array from left-to-right, and want to not reverse
  // if we've been iterating right-to-left.
  reset(reverse, ...arrays) {
    if (this.backIndex < this.cap) {
      // move the filled-in elements from the end
      // to the front of the array
      for (let i = 0; i < arrays.length; i++) {
        const arr = arrays[i];
        const sub = arr.subarray(this.backIndex, this.cap);
        // reverse (if needed)
        if (reverse) sub.reverse();
        // move right elements to follow the left element directly
        arr.set(sub, this.frontIndex);
      }
    }
    // apply the same logical change to the last and length markers
    this.length = this.frontIndex + (this.cap - this.backIndex);
    this.frontIndex = 0;
    this.backIndex = this.cap;
  }
}
