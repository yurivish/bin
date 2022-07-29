# Bin

This is a small project created to see what it's like to use WASM-compiled [Zig](https://ziglang.org) in a browser environment. 

The library implements simple binning into linearly-spaced bins in 1 and 2 dimensions, as well as the Euclidean [Distance Transform](https://en.wikipedia.org/wiki/Distance_transform). 

For details on the algorithm used to compute the EDT, see the paper [Distance Transforms of Sampled Functions](
http://people.cs.uchicago.edu) by Felzenszwalb and Huttenlocher.

For improved efficiency, the EDT implementation also uses an idea from the 1996 paper [An Efficient Algorithm for the Euclidean Distance Transformation](https://onlinelibrary.wiley.com/doi/abs/10.1002/scj.4690270702) in the common case when you're computing the EDT of a binary image.

See the companion [Observable notebook](https://observablehq.com/@yurivish/bin) for the JavaScript API.

For a JavaScript implementation of the EDT, see [this notebook](https://observablehq.com/@mourner/fast-distance-transform) and [library](https://github.com/mapbox/tiny-sdf) by Vladimir Agafonkin.

Build with:

```sh
zig build-lib main.zig --name bin -target wasm32-freestanding -dynamic -OReleaseSafe
````
