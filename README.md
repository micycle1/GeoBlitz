# GeoBlitz

GeoBlitz is an experimental collection of very fast, JTS-inspired spatial indexes and geometry functionality.

### Important notes
- Experimental: APIs and implementations are subject to change.
- Single-package for now: all classes live under `com.github.micycle1.geoblitz` to keep the project compact while iterating.
- Intended use: generally, repeated spatial queries against fixed geometries (build-once, query-many).
- Many GeoBlitz classes are designed as drop-in (or near drop-in) replacements for common JTS spatial operations.

### Quick overview
- Y-stripe point-in-area locator for fast point-in-area tests on complex polygons.
- A Voronoi-derived segment index that maps plane regions to their nearest line segment for fast nearest-segment queries.
- A modernised, typed HPRtree (HHPRtree) with early-exit item visitor and nearest-neighbor search support.

### Comparison vs JTS
| GeoBlitz | JTS Counterpart | Use | Notes | Speedup |  |
|---|---|---|---|---|---|
| `YStripesPointInAreaLocator` | `IndexedPointInAreaLocator` | Point-in-polygon/area test (`contains`, `covers`, etc.) | Per-polygon Y-stripe locators + STRtree for multi-polygons. Reduces candidate segments per query. | ~4x |  |
| `HHPRtree` | `HPRtree` | Spatial index for envelope queries | HHPRtree is typed (`<T>`) variant of HPRtree, and offers a `nearestNeighbor()` method based on "best-first" search | ... |  |
| `SegmentVoronoiIndex` | `IndexedFacetDistance` | Nearest line segment to a coordinate / distance to a polygon | Fast spatial index for approximate nearest-segment queries using a Voronoi-based partitioning of the plane. O(log n) average query. Note: approximates true nearest-segment queries by sampling; accuracy depends on sample spacing. | ~4x |  |
| `OmniUnion` | `CascadedPolygonUnion` | Union of many polygons | Faster union of many polygons by using a spatial index to find likely intersecting polygons. | TBD |  |
| `ConvexHull` | `ConvexHull` | Convex hull of geometry | TBD |  |