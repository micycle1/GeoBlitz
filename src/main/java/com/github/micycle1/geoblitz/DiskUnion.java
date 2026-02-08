package com.github.micycle1.geoblitz;

import org.apache.commons.math3.util.FastMath;
import org.locationtech.jts.geom.*;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Computes the union of circular disks by preserving circular arcs rather than
 * linearizing circles upfront.
 * 
 * <p>
 * This approach is significantly faster than creating polygon approximations of
 * each circle and using JTS {@code CascadedPolygonUnion}. By working with exact
 * circular arcs throughout the computation, it avoids expensive polygon
 * operations on high-vertex-count geometries.
 */
public final class DiskUnion {

	private static final class Disk {
		public final int id;
		public final Coordinate c;
		public final double r;

		/**
		 * @param id     disk identifier
		 * @param circle coordinate where {@code x,y} are center and {@code z} is radius
		 */
		public Disk(int id, Coordinate circle) {
			this.id = id;
			Objects.requireNonNull(circle, "circle");
			this.c = new Coordinate(circle.x, circle.y);
			this.r = circle.getZ();
			if (Double.isNaN(this.r)) {
				throw new IllegalArgumentException("circle.z (radius) must be set");
			}
		}

		public Envelope envelope() {
			return new Envelope(c.x - r, c.x + r, c.y - r, c.y + r);
		}
	}

	/**
	 * A boundary arc of the union, oriented so that union interior is on the left.
	 */
	private static final class Arc {
		public final Disk circle;
		/** start angle in [0, 2π) */
		public final double a0;
		/** end angle in (a0, a0+2π], i.e., CCW from a0 to a1 */
		public final double a1;

		// Filled after snapping:
		SnapVertex start;
		SnapVertex end;

		Arc(Disk circle, double a0, double a1) {
			this.circle = circle;
			this.a0 = a0;
			this.a1 = a1;
		}

		public Coordinate startPoint() {
			return pointOn(circle, a0);
		}

		public Coordinate endPoint() {
			return pointOn(circle, normalizeAngleTo2Pi(a1));
		}

		/** Tangent direction angle at start, for CCW traversal. */
		public double tangentAngleAtStart() {
			// point angle is a0; tangent for CCW is a0 + pi/2
			return normalizeAngleTo2Pi(a0 + Math.PI * 0.5);
		}

		/** Tangent direction angle at end, for CCW traversal. */
		public double tangentAngleAtEnd() {
			double a = normalizeAngleTo2Pi(a1);
			return normalizeAngleTo2Pi(a + Math.PI * 0.5);
		}

		@Override
		public String toString() {
			return "Arc{circle=" + circle.id + ", a0=" + a0 + ", a1=" + a1 + "}";
		}
	}

	/** One closed boundary component: a cyclic list of arcs. */
	private static final class ArcCycle {
		public final List<Arc> arcs;

		ArcCycle(List<Arc> arcs) {
			this.arcs = arcs;
		}
	}

	/** Result: boundary as arc-cycles (outer shells + holes). */
	public static final class ArcBoundary {
		public final List<ArcCycle> cycles;

		ArcBoundary(List<ArcCycle> cycles) {
			this.cycles = cycles;
		}
	}

	/**
	 * Compute union boundary as circular arcs.
	 *
	 * @param circles input circles as Coordinates ({@code x,y} center and {@code z}
	 *                radius)
	 * @param eps     geometric tolerance for snapping / dedup, in world units
	 *                (e.g., 1e-9..1e-6)
	 */
	public static ArcBoundary computeBoundaryArcs(List<Coordinate> circles, double eps) {
		Objects.requireNonNull(circles, "circles");
		if (circles.isEmpty())
			return new ArcBoundary(Collections.emptyList());

		List<Disk> disks = toDisks(circles);

		// 1) Index circles by envelope
		var index = new HPRtreeX<Disk>();
		Map<Integer, Disk> byId = new HashMap<>();
		for (Disk c : disks) {
			byId.put(c.id, c);
			index.insert(c.envelope(), c);
		}
		index.build();

		// 2) For each circle, gather intersection angles
		Map<Integer, List<Double>> anglesByCircle = new HashMap<>();
		for (Disk c : disks)
			anglesByCircle.put(c.id, new ArrayList<>());

		// Also track "definitely contained" circles (fully covered by a single other
		// disk)
		Set<Integer> contained = new HashSet<>();

		// Pair generation using tree pruning:
		// For each circle i, only consider candidate j with j.id > i.id to avoid
		// duplicates.
		List<Disk> sorted = new ArrayList<>(disks);
		sorted.sort(Comparator.comparingInt(cc -> cc.id));

		for (Disk ci : sorted) {
			List<Disk> candidates = index.query(ci.envelope());
			for (Disk cj : candidates) {
				if (cj.id <= ci.id)
					continue;
				// quick bbox prune already done; now do actual circle relation
				CircleRelation rel = circleCircleRelation(ci, cj, eps);
				if (rel.type == RelationType.DISJOINT)
					continue;

				if (rel.type == RelationType.A_INSIDE)
					contained.add(ci.id);
				if (rel.type == RelationType.B_INSIDE)
					contained.add(cj.id);

				if (rel.type == RelationType.SECANT || rel.type == RelationType.TANGENT) {
					// For tangency we still add the same angle once (dedupe later)
					for (Coordinate p : rel.points) {
						if (p == null)
							continue;
						double ai = angleAt(ci, p);
						double aj = angleAt(cj, p);
						anglesByCircle.get(ci.id).add(ai);
						anglesByCircle.get(cj.id).add(aj);
					}
				}
			}
		}

		// 3) Build candidate arcs per circle by splitting at intersection angles
		List<Arc> keptArcs = new ArrayList<>();

		for (Disk c : sorted) {
			if (contained.contains(c.id)) {
				// Contained circles cannot contribute to boundary.
				continue;
			}

			List<Double> angles = anglesByCircle.get(c.id);
			angles = dedupeAndSortAngles(angles, eps);

			if (angles.isEmpty()) {
				// No intersections: either fully covered (contained) or disjoint -> whole
				// circle boundary
				if (!isDiskFullyCovered(c, index, eps)) {
					keptArcs.add(new Arc(c, 0.0, TWO_PI)); // full circle
				}
				continue;
			}

			// Create CCW arcs between consecutive angles, including wrap-around
			for (int k = 0; k < angles.size(); k++) {
				double a0 = angles.get(k);
				double a1 = (k + 1 < angles.size()) ? angles.get(k + 1) : (angles.get(0) + TWO_PI);
				if (a1 <= a0 + angleEps(eps, c.r))
					continue; // degenerate

				Arc arc = new Arc(c, a0, a1);

				// Keep arc if its midpoint is NOT covered by any other disk
				double amid = 0.5 * (a0 + a1);
				Coordinate pmid = pointOn(c, amid);
				if (!isPointCoveredByOtherDisk(c, pmid, index, eps)) {
					keptArcs.add(arc);
				}
			}
		}

		// 4) Snap arc endpoints -> Nodes, and build outgoing adjacency
		VertexSnapper snapper = new VertexSnapper(eps);
		Map<SnapVertex, List<Arc>> outgoing = new HashMap<>();

		for (Arc a : keptArcs) {
			SnapVertex s = snapper.node(a.startPoint());
			SnapVertex e = snapper.node(a.endPoint());
			a.start = s;
			a.end = e;

			outgoing.computeIfAbsent(s, __ -> new ArrayList<>()).add(a);
		}

		// 5) Stitch arcs into cycles by face-walk (keep interior on left)
		Set<Arc> used = Collections.newSetFromMap(new IdentityHashMap<>());
		List<ArcCycle> cycles = new ArrayList<>();

		// Optional: sort outgoing arcs at each node by start tangent angle for faster
		// selection
		for (Map.Entry<SnapVertex, List<Arc>> en : outgoing.entrySet()) {
			en.getValue().sort(Comparator.comparingDouble(Arc::tangentAngleAtStart));
		}

		for (Arc startArc : keptArcs) {
			if (used.contains(startArc))
				continue;
			if (startArc.start == null || startArc.end == null)
				continue;
			// Ignore self-loop full circle arcs for stitching; treat as its own cycle
			if (startArc.a1 - startArc.a0 >= TWO_PI - 1e-12 && startArc.start.equals(startArc.end)) {
				used.add(startArc);
				cycles.add(new ArcCycle(List.of(startArc)));
				continue;
			}

			List<Arc> cycle = new ArrayList<>();
			Arc cur = startArc;
			SnapVertex startNode = startArc.start;

			int guard = 0;
			while (true) {
				if (guard++ > 1_000_000) {
					throw new IllegalStateException("Cycle stitching runaway; check eps/degeneracies.");
				}
				if (used.contains(cur)) {
					// If we re-enter an already-used arc mid-trace, abort this trace.
					break;
				}
				used.add(cur);
				cycle.add(cur);

				SnapVertex at = cur.end;
				if (at.equals(startNode)) {
					cycles.add(new ArcCycle(cycle));
					break;
				}

				List<Arc> outs = outgoing.getOrDefault(at, List.of());
				if (outs.isEmpty()) {
					// Open chain (shouldn't happen for clean union boundary); abandon.
					break;
				}

				double incomingDir = cur.tangentAngleAtEnd();
				Arc next = chooseNextArc(incomingDir, outs, used);

				if (next == null) {
					// No unused continuation; abandon.
					break;
				}
				cur = next;
			}
		}

		return new ArcBoundary(cycles);
	}

	private static List<Disk> toDisks(List<Coordinate> circles) {
		List<Disk> out = new ArrayList<>(circles.size());
		for (int i = 0; i < circles.size(); i++) {
			out.add(new Disk(i, circles.get(i)));
		}
		return out;
	}

	/**
	 * Optional: linearize arc cycles to a Geometry (MultiPolygon) with only final
	 * approximation. This is a lightweight conversion; depending on data, you may
	 * want more robust ring nesting.
	 */
	public static Geometry toJtsGeometry(ArcBoundary boundary, GeometryFactory gf, double maxSegmentLength) {
		List<LinearRing> rings = new ArrayList<>();
		for (ArcCycle cy : boundary.cycles) {
			List<Coordinate> coords = new ArrayList<>();
			for (Arc a : cy.arcs) {
				appendLinearizedArc(coords, a, maxSegmentLength);
			}
			// Close ring
			if (coords.isEmpty())
				continue;
			if (!coords.get(0).equals2D(coords.get(coords.size() - 1))) {
				coords.add(new Coordinate(coords.get(0)));
			}
			// Need at least 4 coords for a LinearRing (closed triangle)
			if (coords.size() < 4)
				continue;
			rings.add(gf.createLinearRing(coords.toArray(new Coordinate[0])));
		}

		// Classify rings into shells vs holes by signed area (CCW => shell, CW => hole)
		List<LinearRing> shells = new ArrayList<>();
		List<LinearRing> holes = new ArrayList<>();
		for (LinearRing r : rings) {
			double a = signedArea(r.getCoordinates());
			if (a >= 0)
				shells.add(r);
			else
				holes.add(r);
		}

		// Naive nesting: assign each hole to the smallest shell that contains its
		// centroid.
		// For complicated cases, use a more robust nesting algorithm.
		List<Polygon> polys = new ArrayList<>();
		List<List<LinearRing>> holesByShell = new ArrayList<>();
		for (int i = 0; i < shells.size(); i++)
			holesByShell.add(new ArrayList<>());

		for (LinearRing h : holes) {
			Point hp = gf.createPoint(h.getCentroid().getCoordinate());
			int bestShell = -1;
			double bestArea = Double.POSITIVE_INFINITY;
			for (int i = 0; i < shells.size(); i++) {
				Polygon shellPoly = gf.createPolygon(shells.get(i));
				if (shellPoly.contains(hp)) { // TODO speed up
					double area = Math.abs(shellPoly.getArea());
					if (area < bestArea) {
						bestArea = area;
						bestShell = i;
					}
				}
			}
			if (bestShell >= 0)
				holesByShell.get(bestShell).add(h);
		}

		for (int i = 0; i < shells.size(); i++) {
			LinearRing shell = shells.get(i);
			LinearRing[] hs = holesByShell.get(i).toArray(new LinearRing[0]);
			polys.add(gf.createPolygon(shell, hs));
		}

		if (polys.isEmpty())
			return gf.createGeometryCollection();
		if (polys.size() == 1)
			return polys.get(0);
		return gf.createMultiPolygon(polys.toArray(new Polygon[0]));
	}

	// ---------------------------- Stitching helpers ----------------------------

	private static Arc chooseNextArc(double incomingDir, List<Arc> outs, Set<Arc> used) {
		Arc best = null;
		double bestDelta = Double.POSITIVE_INFINITY;

		// Choose unused outgoing arc that makes the smallest CCW turn from incomingDir
		// (keeps face interior to the left in typical DCEL "next edge" rule).
		for (Arc cand : outs) {
			if (used.contains(cand))
				continue;
			double outDir = cand.tangentAngleAtStart();
			double delta = normalizeAngleTo2Pi(outDir - incomingDir);
			if (delta < bestDelta) {
				bestDelta = delta;
				best = cand;
			}
		}
		return best;
	}

	// ---------------------------- Coverage tests ----------------------------

	private static boolean isPointCoveredByOtherDisk(Disk self, Coordinate p, HPRtreeX<Disk> index, double eps) {
		Envelope pe = new Envelope(p);
		double px = p.x, py = p.y;
		final boolean[] covered = new boolean[1];
		index.query(pe, c -> {
			if (c.id == self.id) {
				return true;
			}
			double dx = px - c.c.x;
			double dy = py - c.c.y;
			double rr = c.r + eps;
			if (dx * dx + dy * dy <= rr * rr) {
				covered[0] = true;
				return false; // early exit
			}
			return true;
		});
		return covered[0];
	}

	private static boolean isDiskFullyCovered(Disk self, HPRtreeX<Disk> index, double eps) {
		// Disk i is fully covered if exists j such that dist(ci,cj)+ri <= rj (+eps)
		final boolean[] covered = new boolean[1];
		index.query(self.envelope(), c -> {
			if (c.id == self.id) {
				return true;
			}
			double d = dist(self.c, c.c);
			if (d + self.r <= c.r + eps) {
				covered[0] = true;
				return false; // early exit
			}
			return true;
		});
		return covered[0];
	}

	// ---------------------------- Circle-circle relation---------

	private enum RelationType {
		/** The circles do not touch or overlap. */
		DISJOINT,
		/** The circles touch at exactly one point (externally tangent). */
		TANGENT,
		/** The circles intersect at exactly two distinct points. */
		SECANT,
		/** The first circle is completely contained within the second circle. */
		A_INSIDE,
		/** The second circle is completely contained within the first circle. */
		B_INSIDE,
		/** The circles have coincident centers and equal (or nearly equal) radii. */
		COINCIDENT
	}

	private static final class CircleRelation {
		final RelationType type;
		final Coordinate[] points; // size 2, may contain nulls

		CircleRelation(RelationType type, Coordinate p1, Coordinate p2) {
			this.type = type;
			this.points = new Coordinate[] { p1, p2 };
		}

		static CircleRelation disjoint() {
			return new CircleRelation(RelationType.DISJOINT, null, null);
		}
	}

	private static CircleRelation circleCircleRelation(Disk a, Disk b, double eps) {
		double dx = b.c.x - a.c.x;
		double dy = b.c.y - a.c.y;
		double d2 = dx * dx + dy * dy;
		double d = Math.sqrt(d2);

		double ra = a.r, rb = b.r;

		// Coincident centers
		if (d <= eps) {
			if (Math.abs(ra - rb) <= eps)
				return new CircleRelation(RelationType.COINCIDENT, null, null);
			// One contains the other (no boundary intersections)
			if (ra < rb)
				return new CircleRelation(RelationType.A_INSIDE, null, null);
			else
				return new CircleRelation(RelationType.B_INSIDE, null, null);
		}

		// Containment without intersection
		if (d + Math.min(ra, rb) < Math.max(ra, rb) - eps) {
			if (ra < rb)
				return new CircleRelation(RelationType.A_INSIDE, null, null);
			else
				return new CircleRelation(RelationType.B_INSIDE, null, null);
		}

		// Too far apart
		if (d > ra + rb + eps)
			return CircleRelation.disjoint();

		// Compute intersections (including tangency)
		// Using standard formula for two circles.
		double aLen = (ra * ra - rb * rb + d2) / (2.0 * d);
		double h2 = ra * ra - aLen * aLen;

		// Numerical clamp
		if (h2 < 0 && h2 > -eps * eps)
			h2 = 0;

		if (h2 < 0) {
			// No real intersections (should be disjoint or contained; treat as disjoint)
			return CircleRelation.disjoint();
		}

		double xm = a.c.x + aLen * (dx / d);
		double ym = a.c.y + aLen * (dy / d);

		double h = Math.sqrt(h2);
		double rx = -dy * (h / d);
		double ry = dx * (h / d);

		Coordinate p1 = new Coordinate(xm + rx, ym + ry);
		Coordinate p2 = new Coordinate(xm - rx, ym - ry);

		if (h <= eps) {
			// tangent (one intersection)
			return new CircleRelation(RelationType.TANGENT, p1, null);
		}
		return new CircleRelation(RelationType.SECANT, p1, p2);
	}

	// ---------------------------- Angle/arc utilities ----------------------------

	private static final double TWO_PI = Math.PI * 2.0;

	private static double angleAt(Disk c, Coordinate p) {
		double a = FastMath.atan2(p.y - c.c.y, p.x - c.c.x);
		return normalizeAngleTo2Pi(a);
	}

	private static Coordinate pointOn(Disk c, double angle) {
		return pointOn(c, c, angle);
	}

	private static Coordinate pointOn(Disk base, Disk c, double angle) {
		double a = angle;
		return new Coordinate(c.c.x + c.r * Math.cos(a), c.c.y + c.r * Math.sin(a));
	}

	private static double normalizeAngleTo2Pi(double a) {
		double x = a % TWO_PI;
		if (x < 0)
			x += TWO_PI;
		return x;
	}

	private static double angleEps(double eps, double r) {
		// angle tolerance corresponding to linear eps on radius r (approx eps/r)
		if (r <= 0)
			return 1e-12;
		return Math.min(1e-6, Math.max(1e-12, eps / r));
	}

	private static List<Double> dedupeAndSortAngles(List<Double> angles, double eps) {
		if (angles.isEmpty())
			return angles;

		// Normalize & sort
		List<Double> a = angles.stream().map(DiskUnion::normalizeAngleTo2Pi).sorted().collect(Collectors.toCollection(ArrayList::new));

		// Dedupe with angular tolerance (approx)
		// Note: using a fixed tolerance is tricky; this works well in practice when eps
		// is reasonable.
		List<Double> out = new ArrayList<>();
		double tol = 1e-12; // angular tol baseline; we can't infer radius here
		double prev = Double.NaN;
		for (double v : a) {
			if (out.isEmpty() || Math.abs(v - prev) > tol) {
				out.add(v);
				prev = v;
			}
		}
		// Also merge wrap-around duplicates near 0/2π
		if (out.size() >= 2) {
			double first = out.get(0);
			double last = out.get(out.size() - 1);
			if (Math.abs((first + TWO_PI) - last) <= tol) {
				out.set(0, first); // keep first
				out.remove(out.size() - 1);
			}
		}
		return out;
	}

	// ---------------------------- Snapping nodes ----------------------------

	private static final class SnapVertex {
		final long qx, qy; // quantized grid key
		final Coordinate coord; // representative coordinate

		SnapVertex(long qx, long qy, Coordinate coord) {
			this.qx = qx;
			this.qy = qy;
			this.coord = coord;
		}

		@Override
		public boolean equals(Object o) {
			if (this == o)
				return true;
			if (!(o instanceof SnapVertex))
				return false;
			SnapVertex n = (SnapVertex) o;
			return qx == n.qx && qy == n.qy;
		}

		@Override
		public int hashCode() {
			return Objects.hash(qx, qy);
		}
	}

	private static final class VertexSnapper {
		private final double snapTol;
		private final Map<Long, SnapVertex> vertices = new HashMap<>();

		VertexSnapper(double eps) {
			this.snapTol = eps <= 0 ? 1e-12 : eps;
		}

		SnapVertex node(Coordinate p) {
			long qx = quantize(p.x);
			long qy = quantize(p.y);
			long key = interleave(qx, qy);
			SnapVertex n = vertices.get(key);
			if (n != null)
				return n;
			SnapVertex created = new SnapVertex(qx, qy, new Coordinate(p));
			vertices.put(key, created);
			return created;
		}

		private long quantize(double v) {
			return Math.round(v / snapTol);
		}

		// simple reversible mix; not a real Morton code, just stable keying
		private static long interleave(long a, long b) {
			return (a * 73856093L) ^ (b * 19349663L);
		}
	}

	// ---------------------------- Linearization ----------------------------

	private static void appendLinearizedArc(List<Coordinate> out, Arc arc, double maxSegLen) {
		Disk c = arc.circle;
		double a0 = arc.a0;
		double a1 = arc.a1;
		double sweep = a1 - a0;
		if (sweep <= 0)
			return;

		// If full circle
		if (sweep >= TWO_PI - 1e-12) {
			// approximate full circle with N segments
			int n = Math.max(16, (int) Math.ceil((TWO_PI * c.r) / maxSegLen));
			for (int i = 0; i < n; i++) {
				double a = (TWO_PI * i) / n;
				out.add(pointOn(c, a));
			}
			// do not close here; caller closes ring
			return;
		}

		// choose segment count so chord <= maxSegLen
		int n = Math.max(1, (int) Math.ceil((sweep * c.r) / maxSegLen));
		for (int i = 0; i <= n; i++) {
			double a = a0 + (sweep * i) / n;
			Coordinate p = pointOn(c, a);
			// avoid duplicating last point from previous arc
			if (out.isEmpty() || !out.get(out.size() - 1).equals2D(p)) {
				out.add(p);
			}
		}
		// Caller will handle closure and de-dup between arcs
	}

	// ---------------------------- Misc math ----------------------------

	private static double dist(Coordinate a, Coordinate b) {
		return a.distance(b);
	}

	private static double signedArea(Coordinate[] ring) {
		// Shoelace; assumes ring closed or not—works either way
		double s = 0;
		for (int i = 0; i < ring.length - 1; i++) {
			Coordinate p = ring[i];
			Coordinate q = ring[i + 1];
			s += (p.x * q.y - q.x * p.y);
		}
		return 0.5 * s;
	}
}
