function labels = voronoi_labeln(centroids, shape, origin, element_length, method)
%{
Creates an N-D Voronoi diagram label matrix from input coordinates. May be used directly
on pixel coordinate centroids, or used on arbitrary numeric coordinate
centroids.

Inputs:
1. centroids - A M-by-N numeric array of centroid coordinates. M is centroids, N
   is spatial dimensions.
2. shape - A numeric N-vector giving the desired output shape in pixels.
   Optional, default max of input centroids in pixel units, use []. If a shape
   is provided, all centroids MUST lie within the output image. This is checked
   by converting centroids to pixel coordinates using origin and element_length,
   then checking 1 <= centroids <= shape.
3. origin - A numeric N-vector giving the origin of the images in centroid
   coordinates. Optional, defaults to zero vector, use [].
4. element_length - A numeric scalar giving the side length of a pixel in
   centroid coordinates. Optional, defaults to 1.0, use [].
5. method - A scalar string-like denoting which method to use. Must be one of
   "coordinates" or "watershed". Optional, defaults to "coordinates" for
   accuracy.

Two methods are available:
1. "coordinates" - uses a compute- and memory-intensive direct computation of
   nearest centroid to each element. Labels all elements, biased towards smaller
   valued labels when there is a tie. Complexity is O(M*N*P) where M is the
   number of centroids, N is the number of dimensions, P is the number of
   elements in the output image.
2. "watershed" - uses a faster, less-memory-intensive computation using
   optimized image processing algorithms bwdist() and watershed(). Boundaries
   between neighboring Voronoi cells are unlabeled, even when there is no
   ambiguity about distance. Prefer this if boundary cells are unimportant.

Output:
1. labels - A uint32 array with the same size as input "shape".

Usage:
1. Pixel coordinates, call with 1st arg only, or with 1st and 5th (args 2-4
   empty arrays).
2. Arbitrary coordinates, call with 2nd arg to give shape, 3rd to give origin
   and 4th to scale pixels/voxels.
%}

assert(ismatrix(centroids));
assert(isnumeric(centroids));
assert(isreal(centroids));
assert(all(isfinite(centroids), "all"));

CENTROIDS = 1;
DIMENSIONS = 2;

centroid_count = size(centroids, CENTROIDS);
dimension_count = size(centroids, DIMENSIONS);

if nargin < 2 || isempty(shape)
    shape = max(centroids, [], CENTROIDS);
end
assert(isvector(shape));
assert(isnumeric(shape));
assert(isreal(shape));
assert(all(shape == fix(shape)));

if nargin < 3 || isempty(origin)
    origin = zeros(1, dimension_count);
end
assert(isvector(origin));
assert(isnumeric(origin));
assert(isreal(origin));
assert(all(isfinite(origin)));

if nargin < 4 || isempty(element_length)
    element_length = 1;
end
assert(isscalar(element_length));
assert(isnumeric(element_length));
assert(isreal(element_length));
assert(isfinite(element_length));
assert(0.0 < element_length);

if nargin < 5
    method = "coordinates";
end
method = string(method);
assert(isscalar(method));
assert(any(method == ["coordinates", "watershed"]));

centroids = centroids - origin;
centroids = centroids .* element_length;

% sanity checks
assert(all(1 <= centroids, "all"));
assert(all(centroids <= shape, "all"));

if strcmpi(method, "coordinates")
    labels = direct_compute_labels(centroids, shape);
elseif strcmpi(method, "watershed")
    labels = compute_labels_watershed(centroids, shape);
    indices = find(labels == 0);
    labels(indices) = direct_compute_labels(centroids, shape, indices);
else
    assert(false);
end

labels = uint32(labels);
labels = squeeze(labels);

assert(all(size(labels) == shape));
assert(isnumeric(labels));
assert(isreal(labels));
assert(all(labels == fix(labels), "all"));
assert(isempty(setdiff(1 : centroid_count, unique(labels))));

end


function labels = direct_compute_labels(centroids, shape, flat_indices)
%{
Here we go on a wild ride through layers of dimensionality!

For N-D centroids, we end up doing calculations in N+2 dimensions. One through N
are spatial dimensions. N+1 are coordinates at a given spatial location. N+2 are
the list of centroids.

So we end up with variable X as the actual coordinates of each point in space
with size (X1,...,XN,N,1) for dimension N. These are implicitly replicated over
all centroids as the trailing size of 1. Similarly, C is the centroid coordinates, replicated over each point in space with size
(1,...,1,N,M) for M centroids of dimension N, where there are N leading sizes of
1. The replication over space is implicit via the leading 1's. When we subtract
C from X, MATLAB's implicit expansion fills in all the 1's correctly.

If flat_indices is a vector of indices into the image implied by shape, then
computation is performed only for those indices.
%}

if nargin < 3
    flat_indices = [];
end

if ~isempty(flat_indices)
    assert(isvector(flat_indices));
    assert(isnumeric(flat_indices));
    assert(all(flat_indices == fix(flat_indices)));
    assert(all(1 <= flat_indices));
    assert(all(flat_indices <= prod(shape)));
end

centroid_count = size(centroids, CENTROIDS);
dimension_count = size(centroids, DIMENSIONS);

X = cell(1, dimension_count);
if ~isempty(flat_indices)
    [X{:}] = ind2sub(shape, flat_indices);
    X = cat(DIMENSIONS, X{:});
    C = reshape(centroids.', [1 dimension_count centroid_count]);
else
    S = arrayfun(@(x)1:x, shape, "uniformoutput", false);
    [X{:}] = ndgrid(S{:});
    X = cat(dimension_count + 1, X{:});
    C = reshape(centroids.', [ones(1, dimension_count) dimension_count centroid_count]);
end
labels = find_min_indices(X, C);

end


function labels = compute_labels_watershed(centroids, shape)
%{
Computes labels using a distance transform followed by watershed. Note the ridge
lines are assigned zero, so we need to directly compute those boundaries.
%}

centroid_count = size(centroids, CENTROIDS);

im = false(shape);
for i = 1 : centroid_count
    indices = num2cell(centroids(i, :));
    im(indices{:}) = true;
end
im = bwdist(im);
labels = watershed(im);
    
end


function indices = find_min_indices(X, C)
%{
X has size (S,N,1) where S are spatial dimensions (one or more), N are coordinates, 1 may be
implicit.

C has size (1's,N,M) where there are as many 1's as S above, N as above, M
arbitrary.
%}

assert(2 <= ndims(X));
dimension_count = ndims(X) - 1;
assert(ndims(C) - 2 == dimension_count);

d = X - C; % difference of coordinates by implicit expansion
d = d .^ 2; % euclidean distance (we don't need sqrt when taking min)
d = sum(d, dimension_count + 1); % sum along coordinates
[~, indices] = min(d, [], dimension_count + 2); % min along centroids

space_count = prod(size(X, 1 : dimension_count));
assert(numel(indices) == space_count);

end


function n = CENTROIDS()

n = 1;

end


function n = DIMENSIONS()

n = 2;

end
