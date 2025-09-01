package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import org.locationtech.jts.algorithm.Distance;
import org.locationtech.jts.algorithm.locate.PointOnGeometryLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineSegment;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.overlayng.CoverageUnion;
import org.locationtech.jts.triangulate.VoronoiDiagramBuilder;

/**
 * Approximate nearest-line-segment spatial index using a Voronoi diagram
 * constructed from sample points along each segment.
 *
 * Build once, then answer point queries in O(log n) average time.
 * 
 * @author Michael Carleton
 */
public class SegmentVoronoiIndex {

	private final GeometryFactory gf;
	private final Envelope clipEnvelope;
	private final List<LineSegment> segments;
	private final HHPRtree<CellRecord> cellIndex; // indexes Voronoi cells

	public SegmentVoronoiIndex(Polygon polygon, Envelope clipEnvelope, double sampleSpacing) {
		this(segmentsFromPolygon(polygon, false), clipEnvelope, sampleSpacing);
	}

	/**
	 * Build the index.
	 *
	 * @param segments      segments to index (in any order)
	 * @param clipEnvelope  envelope of the plane region to cover
	 * @param sampleSpacing maximum spacing (in coordinate units) between sample
	 *                      points along segments smaller values = more accurate
	 *                      Voronoi; larger values = faster build
	 * @param containsEps   tolerance (in coordinate units) when testing if a point
	 *                      lies in a cell (useful on boundaries)
	 */
	public SegmentVoronoiIndex(List<LineSegment> segments, Envelope clipEnvelope, double sampleSpacing) {
		this.gf = new GeometryFactory();
		this.clipEnvelope = new Envelope(clipEnvelope);
		this.segments = new ArrayList<>(segments);
		this.cellIndex = new HHPRtree<>();

		build(sampleSpacing);
		this.cellIndex.build();
	}

	/**
	 * Returns the nearest line segment to the given point, or null if the point
	 * lies outside the clip envelope.
	 */
	public LineSegment nearestSegment(Coordinate p) {
		if (!clipEnvelope.contains(p)) {
			return null;
		}

		var e = new Envelope(p);
		final CellRecord[] hitBox = new CellRecord[1];
		cellIndex.query(e, cr -> {
			if (cr.contains(p)) {
				hitBox[0] = cr;
				return false;
			}
			return true;
		});

		if (hitBox[0] != null) {
			return hitBox[0].segment;
		}

		// Fallback: scan nearest cell by distance if not strictly contained (extremely
		// rare)
		Point pt = gf.createPoint(p);
		CellRecord best = null;
		double bestDist = Double.POSITIVE_INFINITY;
		for (CellRecord cr : cellIndex.query(e)) {
			double d = cr.cell.distance(pt);
			if (d < bestDist) {
				bestDist = d;
				best = cr;
			}
		}
		return best != null ? best.segment : null;
	}

	/**
	 * Returns the distance from the point to the nearest line segment, or +Inf if
	 * not found.
	 */
	public double distanceToNearestSegment(Coordinate p) {
		LineSegment seg = nearestSegment(p);
		return seg == null ? Double.POSITIVE_INFINITY : Distance.pointToSegment(p, seg.p0, seg.p1);
	}

	private void build(double sampleSpacing) {
		if (segments.isEmpty())
			return;

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

		if (sampleCoords.isEmpty())
			return;

		// 2) Build Voronoi
		VoronoiDiagramBuilder vdb = new VoronoiDiagramBuilder();
		vdb.setSites(gf.createMultiPointFromCoords(sampleCoords.toArray(new Coordinate[0])));
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
				System.out.println("?");
				site = cell.getInteriorPoint().getCoordinate();
			}

			Integer segIdx = siteKeyToSegIdx.get(site);
//			if (segIdx == null) {
//				segIdx = nearestByKey(site, siteKeyToSegIdx);
//			}
			if (segIdx == null) {
				System.out.println("?");
				continue;
			}

			segToCells.computeIfAbsent(segIdx, k -> new ArrayList<>()).add(cell);
		}

		// 4) Union cells per segment to one blob (Polygon), index once per segment
		var cellBlobs = segToCells.entrySet().parallelStream().map(e -> {
			List<Polygon> cells = e.getValue();
			if (cells.isEmpty()) {
				return null;
			}

			Polygon single = (Polygon) CoverageUnion.union(gf.buildGeometry(cells));
			if (single == null || single.isEmpty()) {
				return null;
			}
			LineSegment seg = segments.get(e.getKey());
			CellRecord cr = new CellRecord(single, seg);
			return cr;
		}).filter(Objects::nonNull).toList();

		cellBlobs.forEach(cb -> cellIndex.insert(cb.e, cb));

	}

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

	private static final class CellRecord {
		private final Polygon cell;
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