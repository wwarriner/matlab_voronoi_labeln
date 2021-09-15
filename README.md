# Voronoi Label N
Tool for drawing a Voronoi label matrix given a set of centroids. Creates an N-D Voronoi diagram label matrix from input centroid coordinates. May be used directly on pixel coordinate centroids, or used on arbitrary numeric coordinate centroids with args 2-4.

## Two methods are available:
1. `"coordinates"` - uses a compute- and memory-intensive direct computation of nearest centroid to each element. Labels all elements, biased towards smaller valued labels when there is a tie. Complexity is `O(M*N*P)` where `M` is the number of centroids, `N` is the number of dimensions, `P` is the number of elements in the output image. Suitable if Image Processing Toolbox is not available.
2. `"watershed"` - uses a faster, less-memory-intensive computation using optimized image processing algorithms `bwdist()` and `watershed()`. Prefer this method if Image Processing Toolbox is available, it is about twice as fast.

## Usage:
1. Pixel coordinates, call with 1st arg only, or with 1st and 5th (args 2-4 empty arrays).
2. Arbitrary coordinates, call with 2nd arg to give shape, 3rd to give origin and 4th to scale pixels/voxels.

```matlab
% direct
labels = voronoi_labeln(centroids) % implied shape
labels = voronoi_labeln(centroids, shape) % explicit shape
labels = voronoi_labeln(centroids, shape, [], px_len) % scale
labels = voronoi_labeln(centroids, shape, origin, []) % translate
labels = voronoi_labeln(centroids, shape, origin, px_len) % both

% watershed
labels = voronoi_labels(centroids, [], [], [], "watershed")
labels = voronoi_labels(centroids, shape, [], [], "watershed")
% ... etc.
```