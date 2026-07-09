package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.index.strtree.STRtree;
import org.locationtech.jts.noding.MCIndexNoder;
import org.locationtech.jts.noding.NodedSegmentString;
import org.locationtech.jts.noding.Noder;
import org.locationtech.jts.noding.snapround.SnapRoundingNoder;
import org.locationtech.jts.operation.polygonize.Polygonizer;

public final class NAryUnion {

	/*-
	 * EdgeSetIntersector
	 */

	private NAryUnion() {
	}

	// Convenience: union of any polygonal inputs. Non-polygonal are ignored.
//	public static Geometry union(Collection<Geometry> inputs) {
//		if (inputs == null || inputs.isEmpty()) {
//			return new GeometryFactory().createEmpty(2);
//		}
//		GeometryFactory gf = inputs.iterator().next().getFactory();
//		List<Polygon> polys = new ArrayList<>();
//		PolygonExtracter.getPolygons(inputs, polys);
//		if (polys.isEmpty())
//			return gf.createEmpty(2);
//		if (polys.size() == 1)
//			return polys.get(0);
//
//		// If your coordinates are noisy (binary FP with excess decimals),
//		// set a fixed precision model (e.g., scale 1e8) and pass it below.
//		return unionPolygons(polys, null);
//	}

	// Main entry: pass Polygon list and optional fixed precision (recommended).
	public static Geometry union(List<Polygon> polys, PrecisionModel pm) {
		if (polys.isEmpty()) {
			return new GeometryFactory().createEmpty(2);
		}
		GeometryFactory gf = polys.get(0).getFactory();

		// Optional defensive fix for rare invalid inputs (fast-skip valid ones)
//		polys = fixInvalid(polys);

		// 1) Node all rings once (snap-round if pm != null; else robust index noder)
		List<NodedSegmentString> noded = nodeAllRings(polys, pm);

		// 2) Polygonize the arrangement edges into bounded faces
		List<LineString> nodedLines = toLineStrings(noded, gf);
		Polygonizer polyz = new Polygonizer();
		polyz.add(nodedLines);
		@SuppressWarnings("unchecked")
		Collection<Polygon> facesAll = polyz.getPolygons();

		if (facesAll.isEmpty()) {
			return gf.createEmpty(2);
		}

		// 3) Keep only faces covered by at least one input polygon
		// Build a spatial index of prepared polygons for fast point-in-polygon
		STRtree index = new STRtree();
		List<PreparedGeometry> prepared = new ArrayList<>(polys.size());
		for (Polygon p : polys) {
			PreparedGeometry prep = PreparedGeometryFactory.prepare(p);
			prepared.add(prep);
			index.insert(p.getEnvelopeInternal(), prep);
		}

		List<Polygon> keptFaces = new ArrayList<>();
		for (Polygon f : facesAll) {
			// representative interior point
			Point ip = f.getInteriorPoint();
			@SuppressWarnings("unchecked")
			List<PreparedGeometry> cands = index.query(new Envelope(ip.getX(), ip.getX(), ip.getY(), ip.getY()));
			boolean covered = false;
			for (PreparedGeometry prep : cands) {
				if (prep.covers(ip)) { // covers handles boundary cases
					covered = true;
					break;
				}
			}
			if (covered) {
				keptFaces.add(f);
			}
		}
		if (keptFaces.isEmpty()) {
			return gf.createEmpty(2);
		}
		if (keptFaces.size() == 1) {
			return keptFaces.get(0);
		}

		// 4) Dissolve internal edges by canceling doubled segments
		// Count undirected segments over all kept faces
		Map<EdgeKey, Integer> edgeCount = new HashMap<>(keptFaces.size() * 8);
		for (Polygon f : keptFaces) {
			addRingEdges(f.getExteriorRing(), edgeCount);
			int nh = f.getNumInteriorRing();
			for (int i = 0; i < nh; i++) {
				// faces from polygonizer should be hole-free; safeguard:
				addRingEdges(f.getInteriorRingN(i), edgeCount);
			}
		}

		// Keep only boundary edges (odd count == 1 after pairing)
		List<LineString> boundarySegments = new ArrayList<>();
		for (Map.Entry<EdgeKey, Integer> e : edgeCount.entrySet()) {
			if ((e.getValue() & 1) == 1) { // count == 1
				boundarySegments.add(e.getKey().toLineString(gf));
			}
		}
		if (boundarySegments.isEmpty()) {
			return gf.createEmpty(2);
		}

		// 5) Polygonize boundary to get dissolved union polygons
		Polygonizer boundaryPolyz = new Polygonizer();
		boundaryPolyz.add(boundarySegments);
		@SuppressWarnings("unchecked")
		Collection<Polygon> dissolved = boundaryPolyz.getPolygons();

		if (dissolved.isEmpty()) {
			return gf.createEmpty(2);
		}
		if (dissolved.size() == 1) {
			return dissolved.iterator().next();
		}
		return gf.buildGeometry(new ArrayList<>(dissolved));
	}

	// Node all rings into segment strings
	private static List<NodedSegmentString> nodeAllRings(List<Polygon> polys, PrecisionModel pm) {
		List<NodedSegmentString> segs = new ArrayList<>();
		for (Polygon p : polys) {
			addRing(p.getExteriorRing(), segs);
			int nh = p.getNumInteriorRing();
			for (int i = 0; i < nh; i++) {
				addRing(p.getInteriorRingN(i), segs);
			}
		}
		Noder noder = (pm != null) ? new SnapRoundingNoder(pm) : new MCIndexNoder();
		noder.computeNodes(segs);
		@SuppressWarnings("unchecked")
		List<NodedSegmentString> noded = (List<NodedSegmentString>) noder.getNodedSubstrings();
		return noded;
	}

	private static void addRing(LineString ring, List<NodedSegmentString> out) {
		Coordinate[] pts = ring.getCoordinates();
		if (pts.length < 2) {
			return;
		}
		// Ensure closed: JTS rings are closed; just be defensive
		if (!pts[0].equals2D(pts[pts.length - 1])) {
			Coordinate[] closed = Arrays.copyOf(pts, pts.length + 1);
			closed[closed.length - 1] = closed[0].copy();
			pts = closed;
		}
		out.add(new NodedSegmentString(pts, null));
	}

	private static List<LineString> toLineStrings(List<NodedSegmentString> segs, GeometryFactory gf) {
		List<LineString> ls = new ArrayList<>(segs.size());
		for (NodedSegmentString ss : segs) {
			Coordinate[] c = ss.getCoordinates();
			if (c.length >= 2) {
				ls.add(gf.createLineString(c));
			}
		}
		return ls;
	}

	// Count edges of a ring as undirected segments
	private static void addRingEdges(LineString ring, Map<EdgeKey, Integer> count) {
		Coordinate[] c = ring.getCoordinates();
		for (int i = 1; i < c.length; i++) {
			EdgeKey k = new EdgeKey(c[i - 1], c[i]);
			count.merge(k, 1, Integer::sum);
		}
	}

	// Undirected edge key, robust under snap rounding (exact doubles on a fixed
	// grid)
	private static final class EdgeKey {
		final double ax, ay, bx, by;

		EdgeKey(Coordinate a, Coordinate b) {
			// Normalize undirected: (min, max) lexicographically
			int cmp = cmp(a, b);
			if (cmp <= 0) {
				this.ax = a.x;
				this.ay = a.y;
				this.bx = b.x;
				this.by = b.y;
			} else {
				this.ax = b.x;
				this.ay = b.y;
				this.bx = a.x;
				this.by = a.y;
			}
		}

		static int cmp(Coordinate p, Coordinate q) {
			if (p.x < q.x) {
				return -1;
			}
			if (p.x > q.x) {
				return 1;
			}
			return Double.compare(p.y, q.y);
		}

		@Override
		public boolean equals(Object o) {
			if (this == o) {
				return true;
			}
			if (!(o instanceof EdgeKey)) {
				return false;
			}
			EdgeKey ek = (EdgeKey) o;
			return ax == ek.ax && ay == ek.ay && bx == ek.bx && by == ek.by;
		}

		@Override
		public int hashCode() {
			long h = 1469598103934665603L;
			h = mix(h, ax);
			h = mix(h, ay);
			h = mix(h, bx);
			h = mix(h, by);
			return (int) (h ^ (h >>> 32));
		}

		private static long mix(long h, double v) {
			long x = Double.doubleToLongBits(v);
			h ^= x;
			h *= 1099511628211L;
			return h;
		}

		LineString toLineString(GeometryFactory gf) {
			return gf.createLineString(new Coordinate[] { new Coordinate(ax, ay), new Coordinate(bx, by) });
		}
	}

	// Optional: fix only the invalid ones (cheap guard)
//	private static List<Polygon> fixInvalid(List<Polygon> polys) {
//		List<Polygon> out = new ArrayList<>(polys.size());
//		for (Polygon p : polys) {
//			if (p.isEmpty())
//				continue;
//			if (p.isValid()) {
//				out.add(p);
//			} else {
//				// GeometryFixer produces polygonal; safer than buffer(0) for some cases
//				Geometry fixed = org.locationtech.jts.geom.util.GeometryFixer.fix(p);
//				List<Polygon> extracted = new ArrayList<>();
//				PolygonExtracter.getPolygons(Collections.singletonList(fixed), extracted);
//				out.addAll(extracted);
//			}
//		}
//		return out;
//	}
}