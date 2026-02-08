[![](https://jitpack.io/v/micycle1/GeoBlitz.svg)](https://jitpack.io/#micycle1/GeoBlitz)

# GeoBlitz

GeoBlitz is an experimental collection of very fast, JTS-inspired spatial indexes and geometry utilities.

Designed for **build-once, query-many** scenarios, many classes serve as drop-in (or near drop-in) replacements for common [JTS](https://github.com/locationtech/jts/) operations, but with significantly better performance.

> **Important Notes**  
> - **Experimental (for now)**: APIs and implementations may change.
> - **Single-package**: All classes live under `com.github.micycle1.geoblitz` during early development.
> - **Target Use**: Repeated spatial queries against fixed geometries (build once, query many).

---

### Core Indexes / Helpers

#### `FastConvexHull`
A compact, low-allocation convex hull implementation based on Andrew's monotone chain optimised for moderate-size inputs.
- Performs in-place sorting and stack-style upper/lower hull construction to reduce memory churn.
- Focuses on minimising allocations for faster execution on moderate-sized point sets.
- Note: slower than JTS ConvexHull on extremely large inputs (≈1M+ points).

#### `HilbertParallelPolygonUnion`
High-performance polygon union using Hilbert curve ordering and parallel reduction.
- Sorts geometries by a Hilbert-encoded key derived from envelopes to improve spatial locality.
- Executes the union as a parallel reduction so many union tasks run concurrently.
- Discards non-polygonal artifacts; returns clean `Polygon` or `MultiPolygon`.

#### `DiskUnion`
Fast union of disks/circles by preserving circular arcs in the boundary representation instead of linearizing circles upfront.
- Significantly faster than polygonizing circles and using JTS `CascadedPolygonUnion`, especially at high circle resolution.

#### `FastVariableBuffer`
Identical to JTS VariableBuffer but uses HilbertParallelPolygonUnion for faster unioning of buffer components.

#### `HPRtreeX`
A typed Hilbert-packed R-tree with extended query capabilities.
- Generic type support (`<T>`) for storing arbitrary user objects.
- Early-exit item visitors to allow immediate termination of queries.
- Efficient nearest-neighbor (best-first) and range (depth-first) queries with bounding-box pruning.

#### `IndexedLengthIndexedLine`
Length-based linear-referencing with O(log n) length→location lookups for repeated queries.
- Length-based linear referencing with O(log n) index for repeated `extractPoint()`/`extractLine()` calls.
- Precomputes cumulative segment lengths.  
- Ideal for large linestrings with many queries – much faster than JTS’s linear-scan `LengthIndexedLine`.

#### `FastLineIntersector`
High-throughput line-segment intersection that trades the absolute robustness of the JTS robust intersector for speed.
- Skips expensive DD (double-double) arithmetic and full determinant expansion in the common case; falls back to a robust determinant predicate only when a tiny relative epsilon indicates near-degeneracy.
- Skips coordinate z-interpolation.
- Use case: many short-lived segment intersection computations where throughput matters and the input is not adversarial. Not recommended when absolute bitwise-robustness for degenerate, adversarial or extremely large-scale inputs is required.

#### `SegmentVoronoiIndex`
Approximate nearest-segment index built from Voronoi cells of sampled points along input segments.
- Samples each segment at configurable spacing to trade accuracy for build time and memory.
- Builds a Voronoi diagram from samples and unions cells belonging to the same segment.
- Answers nearest-segment queries via point-in-cell tests using `YStripesPointInAreaLocator`.

#### `YStripesPointInAreaLocator`
Fast point-in-area locator using per-polygon Y-stripe indexes and an HPR-tree for candidate selection.
- Single-polygon fast path; HPR-tree of per-polygon locators for multi-polygons.
- Immutable and safe for concurrent repeated point-in-area queries.
- Inspired by the [tg library’s YStripes](https://github.com/tidwall/tg/blob/main/docs/POLYGON_INDEXING.md#ystripes).

#### `ProHausdorffDistance`
Projection-based Hausdorff Distance approximation. Implements [ProHD](https://www.arxiv.org/abs/2511.18207).
- Projects geometry points onto informative directions (centroid axis, PCA axis) and computes Hausdorff distance on a small subset of extreme points.
- Extremely fast estimate of the true Hausdorff distance.
- Adjustable `alpha` parameter controls accuracy vs speed.

#### `PointDistanceIndex`
Fast exact distance queries from points to a boundary and optional obstacles, with optional signed distance determined by the boundary.

- Uses `HPRtree` to index line “facets” (small multi-segment chunks) rather than individual segments, reducing index size and depth.
- Targets include boundary rings plus obstacle lines and points.
- Supports unsigned distance (always non-negative) and optional signed distance (negative if outside the boundary).


#### `EndpointSnapper`
Endpoint-only snapping for near-coverage linework (with optional polygon-vertex anchoring).
- Snaps only LineString *endpoints* together within a tolerance; interior line vertices and polygon vertices are never moved.
- Uses transitive clustering to close tiny endpoint gaps and near-misses.
- Preserves structure: returns snapped LineStrings; Polygons/MultiPolygons (and other types) are passed through unchanged.
- Intended for line arrangements that should form a valid coverage but don’t quite join; typically followed by noding and polygonization to build faces. Not a general vertex/edge snapper.

### Comparison table
*For classes with a JTS counterpart*
| GeoBlitz Class | JTS Counterpart | Speedup |
|---|---|---:|
| FastConvexHull | ConvexHull | TBD |
| HilbertParallelPolygonUnion | CascadedPolygonUnion | TBD |
| DiskUnion | CascadedPolygonUnion | ~30x |
| HPRtreeX | HPRtree | *provides NN/range/early-exit features* |
| IndexedLengthIndexedLine | LengthIndexedLine | O(log n) queries vs O(n) scan – large speedups for repeated queries |
| FastLineIntersector | RobustLineIntersector | ~2x |
| SegmentVoronoiIndex | IndexedFacetDistance | ~4x (dataset & sampling dependent) |
| YStripesPointInAreaLocator | IndexedPointInAreaLocator | ~4x |
| ProHausdorffDistance | DiscreteHausdorffDistance | 10x+ |
| PointDistanceIndex | IndexedFacetDistance | ~2x |


### Benchmarks
Include `jmh:benchmark` as a Maven goal to run the benchmarks (see pom.xml).

Reported speedups depend on dataset, geometry complexity, and parameters (e.g., sampling spacing); please run the included JMH benchmarks on your data for authoritative numbers.
