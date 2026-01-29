package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.locationtech.jts.algorithm.Distance;
import org.locationtech.jts.algorithm.locate.PointOnGeometryLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.overlayng.CoverageUnion;
import org.locationtech.jts.triangulate.VoronoiDiagramBuilder;

/**
 * Approximate nearest-line-segment spatial index built from sample points taken
 * along each segment and using the resulting Voronoi partition.
 * <p>
 * Build the index once (from a Polygon or a List&lt;LineSegment&gt;) and then
 * answer nearest-segment queries quickly. The sampling spacing controls the
 * trade-off between accuracy and build time: smaller spacing = higher accuracy,
 * larger spacing = faster build and lower memory use.
 * </p>
 * <p>
 * Queries: {@link #nearestSegment(Coordinate)} returns the original LineSegment
 * whose Voronoi cell contains the point, or {@code null} if the point lies
 * outside the configured clip envelope.
 * {@link #distanceToNearestSegment(Coordinate)} returns the Euclidean distance
 * to that segment (or +Inf if none).
 * </p>
 * <p>
 * Limitations: this is an approximate index — results near boundaries may be
 * affected by sampling density. The index is intended for read-only concurrent
 * queries after construction; it is not designed for incremental modification.
 * </p>
 *
 * @author Michael Carleton
 */
public class SegmentVoronoiIndex {

	private final GeometryFactory gf;
	private Envelope clipEnvelope;
	private final List<LineSegment> segments;
	/**
	 * The spatial index of Voronoi cells.
	 */
	public final HPRtreeX<CellRecord> cellIndex; // indexes Voronoi cells

	/**
	 * Constructs a SegmentVoronoiIndex from the given polygon. The index is built
	 * during construction and is ready for concurrent read-only queries after the
	 * constructor returns.
	 * <p>
	 * The clip envelope used to clip the Voronoi diagram is computed from the
	 * polygon's envelope and expanded by one quarter of its diameter to reduce edge
	 * effects.
	 *
	 * @param polygon       the polygon whose exterior ring will be converted to
	 *                      line segments and used to build the index (must not be
	 *                      null)
	 * @param sampleSpacing maximum spacing (in coordinate units) between sample
	 *                      points along segments; smaller values increase accuracy
	 *                      and memory/time cost, larger values make the build
	 *                      faster and use less memory
	 */
	public SegmentVoronoiIndex(Polygon polygon, double sampleSpacing) {
		this.gf = new GeometryFactory();
		this.clipEnvelope = polygon.getEnvelopeInternal();
		this.clipEnvelope.expandBy(clipEnvelope.getDiameter() / 4);
		var segments = segmentsFromPolygon(polygon, true);
		this.segments = new ArrayList<>(segments);
		this.cellIndex = new HPRtreeX<>();

		build(sampleSpacing);
		this.cellIndex.build();

	}

	/**
	 * Constructs a SegmentVoronoiIndex from the given polygon. The index is built
	 * during construction and is ready for concurrent read-only queries after the
	 * constructor returns.
	 * <p>
	 * The provided clipEnvelope is used to clip the Voronoi diagram instead of
	 * deriving one from the polygon.
	 *
	 * @param polygon       the polygon whose rings will be converted to line
	 *                      segments and used to build the index (must not be null)
	 * @param clipEnvelope  the envelope that defines the region to cover; may be
	 *                      null to indicate no clipping
	 * @param sampleSpacing maximum spacing (in coordinate units) between sample
	 *                      points along segments; smaller values increase accuracy
	 *                      and memory/time cost, larger values make the build
	 *                      faster and use less memory
	 */
	public SegmentVoronoiIndex(Polygon polygon, Envelope clipEnvelope, double sampleSpacing) {
		this(segmentsFromPolygon(polygon, true), clipEnvelope, sampleSpacing);
	}

	/**
	 * Constructs a SegmentVoronoiIndex from an explicit list of line segments and
	 * builds a Voronoi-based approximate nearest-segment index. The index is built
	 * during construction and is ready for concurrent read-only queries after the
	 * constructor returns.
	 * <p>
	 * The supplied list of segments is copied; segments may be provided in any
	 * order. If clipEnvelope is null, no clipping envelope is applied when
	 * constructing the Voronoi diagram.
	 *
	 * @param segments      the list of LineSegment objects to index (must not be
	 *                      null; may be empty)
	 * @param clipEnvelope  the envelope that defines the region to cover; may be
	 *                      null to indicate no clipping
	 * @param sampleSpacing maximum spacing (in coordinate units) between sample
	 *                      points along segments; smaller values increase accuracy
	 *                      and memory/time cost, larger values make the build
	 *                      faster and use less memory
	 */
	public SegmentVoronoiIndex(List<LineSegment> segments, Envelope clipEnvelope, double sampleSpacing) {
		this.gf = new GeometryFactory();
		this.clipEnvelope = clipEnvelope == null ? null : new Envelope(clipEnvelope);
		this.segments = new ArrayList<>(segments);
		this.cellIndex = new HPRtreeX<>();

		build(sampleSpacing);
		this.cellIndex.build();
	}

	/**
	 * Returns the original LineSegment whose Voronoi cell contains the given point.
	 * If the point lies outside the configured clip envelope, {@code null} is
	 * returned.
	 * <p>
	 * Note: this is an approximate index built from discrete samples along the
	 * segments; results for points near cell boundaries may depend on sampling
	 * density. The method is thread-safe for concurrent readers.
	 *
	 * @param p the query coordinate (must not be null)
	 * @return the LineSegment whose Voronoi cell contains the point, or
	 *         {@code null} if the point lies outside the clip envelope or no
	 *         containing cell was found
	 */
	public LineSegment nearestSegment(Coordinate p) {
		if (!clipEnvelope.contains(p)) {
			// TODO if outside convex hull, nearest lies on convex hull
			return null;
		}

		var e = new Envelope(p);
		final CellRecord[] hitBox = new CellRecord[1];
		cellIndex.query(e, cr -> {
			if (cr.contains(p)) {
				hitBox[0] = cr;
				return false; // terminate query (continue=false)
			}
			return true; // continue
		});

		return hitBox[0].segment; // shouldn't be null?!
	}

	/**
	 * Returns the Euclidean distance from the given point to the nearest indexed
	 * line segment as determined by the Voronoi-based index. If no nearest segment
	 * can be determined (for example, the point is outside the clip envelope),
	 * {@code Double.POSITIVE_INFINITY} is returned.
	 * <p>
	 * This method delegates to {@link #nearestSegment(Coordinate)} and then
	 * computes the distance to the returned segment.
	 *
	 * @param p the query coordinate (must not be null)
	 * @return the Euclidean distance to the nearest indexed LineSegment, or
	 *         {@code Double.POSITIVE_INFINITY} if no segment was found
	 */
	public double distanceToNearestSegment(Coordinate p) {
		LineSegment seg = nearestSegment(p);
		return seg == null ? Double.POSITIVE_INFINITY : Distance.pointToSegment(p, seg.p0, seg.p1);
	}

	private void build(double sampleSpacing) {
		if (segments.isEmpty()) {
			return;
		}

		// 1) Sample each segment
		List<Coordinate> sampleCoords = new ArrayList<>();
		Map<Coordinate, Integer> siteKeyToSegIdx = new HashMap<>();

		for (int si = 0; si < segments.size(); si++) {
			LineSegment seg = segments.get(si);
			List<Coordinate> samples = sampleSegment(seg, sampleSpacing);
			for (Coordinate c : samples) {
				sampleCoords.add(c);
				siteKeyToSegIdx.put(c, si);
			}
		}

		if (sampleCoords.isEmpty()) {
			return;
		}

		// 2) Build Voronoi
		VoronoiDiagramBuilder vdb = new VoronoiDiagramBuilder();
		var sites = gf.createMultiPointFromCoords(sampleCoords.toArray(new Coordinate[0]));
		// setSites de-dedupes points
		vdb.setSites(sites);
		if (clipEnvelope == null) {
			clipEnvelope = sites.getEnvelopeInternal();
			clipEnvelope.expandBy(1e-6);
		}
		vdb.setClipEnvelope(clipEnvelope);

		Geometry diagram = vdb.getDiagram(gf);

		// 3) Collect cells per segment
		Map<Integer, List<Polygon>> segToCells = new HashMap<>(segments.size());

		for (int i = 0; i < diagram.getNumGeometries(); i++) {
			Geometry g = diagram.getGeometryN(i);
			if (!(g instanceof Polygon)) {
				continue;
			}
			Polygon cell = (Polygon) g;

			Coordinate site = null;
			Object ud = cell.getUserData();
			if (ud instanceof Coordinate) {
				site = (Coordinate) ud;
			} else {
				site = cell.getInteriorPoint().getCoordinate();
			}

			Integer segIdx = siteKeyToSegIdx.get(site);

			segToCells.computeIfAbsent(segIdx, k -> new ArrayList<>()).add(cell);
		}

		cellsList = new CopyOnWriteArrayList<>();
		// 4) Union cells per segment to one blob (Polygon), index once per segment
		var cellBlobs = segToCells.entrySet().parallelStream().flatMap(e -> {
			List<Polygon> cells = e.getValue();

			Geometry union = CoverageUnion.union(gf.buildGeometry(cells));
			if (union == null || union.isEmpty()) {
				return Stream.empty();
			}

			LineSegment seg = segments.get(e.getKey());
			List<CellRecord> results = new ArrayList<>();

			// exceedingly rare, but handle
			if (union instanceof MultiPolygon || union instanceof GeometryCollection) {
				int n = union.getNumGeometries();
				for (int i = 0; i < n; i++) {
					Geometry g = union.getGeometryN(i);
					if (g instanceof Polygon && !g.isEmpty()) {
						results.add(new CellRecord((Polygon) g, seg));
					}
				}
			} else {
				Polygon poly = (Polygon) union;
				if (!poly.isEmpty()) {
					results.add(new CellRecord(poly, seg));
				}
			}

			return results.stream();
		}).collect(Collectors.toList());

		cellBlobs.forEach(cb -> cellIndex.insert(cb.e, cb));
	}

	/**
	 * List of lists of polygons representing the cells.
	 */
	public List<List<Polygon>> cellsList;

	private static List<Coordinate> sampleSegment(LineSegment seg, double spacing) {
		double len = seg.getLength();
		if (!(spacing > 0) || len == 0) {
			return Arrays.asList(new Coordinate(seg.p0), new Coordinate(seg.p1));
		}
		int nParts = Math.max(1, (int) Math.ceil(len / spacing));
		List<Coordinate> out = new ArrayList<>(nParts + 1);
		for (int i = 0; i <= nParts; i++) {
			double frac = (double) i / (double) nParts;
			out.add(new Coordinate(seg.pointAlong(frac)));
		}
		return out;
	}

	private static List<LineSegment> segmentsFromPolygon(Polygon poly, boolean includeHoles) {
		List<LineSegment> segs = new ArrayList<>();
		addSegmentsFromLineString(poly.getExteriorRing(), segs);
		if (includeHoles) {
			for (int i = 0; i < poly.getNumInteriorRing(); i++) {
				addSegmentsFromLineString(poly.getInteriorRingN(i), segs);
			}
		}
		return segs;
	}

	private static void addSegmentsFromLineString(LineString ls, List<LineSegment> out) {
		Coordinate[] cs = ls.getCoordinates();
		for (int i = 1; i < cs.length; i++) {
			out.add(new LineSegment(cs[i - 1], cs[i]));
		}
	}

	/**
	 * A record containing the Voronoi cell polygon and its associated segment.
	 */
	public static final class CellRecord {
		/**
		 * The Voronoi cell polygon.
		 */
		public final Polygon cell;
		private final LineSegment segment;
		PointOnGeometryLocator areaLocator;
		Envelope e;

		CellRecord(Polygon cell, LineSegment segment) {
			this.cell = cell;
			this.segment = segment;
			areaLocator = new YStripesPointInAreaLocator(cell);
			e = cell.getEnvelopeInternal();
		}

		boolean contains(Coordinate c) {
			return areaLocator.locate(c) != Location.EXTERIOR;
		}
	}
}