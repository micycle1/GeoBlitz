package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

import org.apache.commons.math3.util.FastMath;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;

/**
 * Computes the union of circular disks, preserving exact circular arcs in the
 * boundary representation and deferring linearization until final geometry
 * construction.
 *
 * <p>
 * This is a planar (Euclidean) port of Volodymyr Agafonkin's
 * <a href="https://github.com/mourner/circle-union">circle-union</a>, which
 * operates on geodesic disks on a sphere. The pipeline is:
 * <ol>
 * <li><b>Build</b> — sort circles by radius descending (so the larger circle of
 * any pair is processed first) and index centers in a flat KD-tree.</li>
 * <li><b>Scan</b> — a single radius-descending sweep classifies every
 * overlapping pair: circles engulfed by a larger one are dropped globally; each
 * properly-intersecting pair solves for its two boundary points, stored once
 * under stable integer IDs shared by both circles, and the pair is unioned in a
 * disjoint-set forest (the overlap graph's connected components are the union's
 * components).</li>
 * <li><b>Arcs</b> — per circle, each intersecting neighbor covers one angular
 * interval of its boundary; a circular +1/−1 depth sweep unions the covered
 * intervals and emits the complement as boundary arcs, endpoints referenced by
 * point ID.</li>
 * <li><b>Stitch</b> — every intersection point is the end of exactly one arc
 * and the start of exactly one other, so rings are walked by exact integer ID
 * handoff (no geometric search). Rings are nested into polygons by
 * connectivity: each disjoint-set component contributes one shell (positive
 * signed area) plus its holes — no point-in-ring containment tests.</li>
 * </ol>
 *
 * <p>
 * Classification is exact floating-point, intersection points are computed
 * once and shared by ID rather than re-derived and snapped, and all hot loops
 * run over flat primitive arrays. Degenerate inputs that break the
 * exact-topology invariants (e.g. three circles passing through the same
 * point) throw {@link IllegalStateException} rather than producing corrupt
 * output. Exact duplicate circles (bit-identical centers) are deduplicated up
 * front; near-duplicates are handled by the regular engulf/intersect logic.
 *
 * <pre>{@code
 * CircleUnion cu = new CircleUnion();
 * for (Coordinate c : circles) {
 * 	cu.add(c.x, c.y, c.getZ());
 * }
 * Geometry union = cu.geometry(0.5);
 * }</pre>
 *
 * @author Michael Carleton
 * @author Volodymyr Agafonkin (original spherical implementation)
 */
public final class CircleUnion {

	private static final double TWO_PI = Math.PI * 2.0;

	// input buffers (insertion order)
	private double[] inX;
	private double[] inY;
	private double[] inR;
	private int pos;

	private static final int DEFAULT_CAPACITY = 8;

	/** Cached arc topology, invalidated by {@link #add}. */
	private Topology topology;

	/**
	 * Creates an empty builder.
	 */
	public CircleUnion() {
		this(DEFAULT_CAPACITY);
	}

	/**
	 * Creates an empty builder with the given initial capacity.
	 *
	 * <p>
	 * This is purely a performance hint to avoid repeated array growth when the
	 * final circle count is known in advance; {@link #add} grows the backing
	 * storage automatically, so any number of circles may be added regardless of
	 * the value passed here.
	 *
	 * @param initialCapacity initial capacity hint
	 * @throws IllegalArgumentException if {@code initialCapacity} is negative
	 */
	public CircleUnion(int initialCapacity) {
		if (initialCapacity < 0) {
			throw new IllegalArgumentException("initialCapacity must be non-negative");
		}
		inX = new double[initialCapacity];
		inY = new double[initialCapacity];
		inR = new double[initialCapacity];
	}

	/**
	 * Adds a circle and returns its insertion index.
	 *
	 * @param x center x
	 * @param y center y
	 * @param r radius (finite, non-negative; zero-radius circles contribute nothing
	 *          to the union)
	 * @return the index of the added circle
	 * @throws IllegalArgumentException if the center is not finite or the radius is
	 *                                  NaN, infinite or negative
	 */
	public int add(double x, double y, double r) {
		if (!Double.isFinite(x) || !Double.isFinite(y)) {
			throw new IllegalArgumentException("circle center must be finite");
		}
		if (!(r >= 0) || Double.isInfinite(r)) {
			throw new IllegalArgumentException("radius must be finite and non-negative");
		}
		if (pos == inX.length) {
			grow();
		}
		int i = pos++;
		inX[i] = x;
		inY[i] = y;
		inR[i] = r;
		topology = null;
		return i;
	}

	private void grow() {
		int newCapacity = Math.max(DEFAULT_CAPACITY, inX.length * 2);
		inX = Arrays.copyOf(inX, newCapacity);
		inY = Arrays.copyOf(inY, newCapacity);
		inR = Arrays.copyOf(inR, newCapacity);
	}

	/**
	 * Returns the number of circles added so far.
	 *
	 * @return circle count
	 */
	public int size() {
		return pos;
	}

	/**
	 * Computes (or reuses the cached) union topology and linearizes it using a
	 * default {@link GeometryFactory}.
	 *
	 * @param maxSegmentLength maximum chord length for linearized arc segments;
	 *                         smaller values produce smoother output
	 * @return {@code Polygon} or {@code MultiPolygon} representing the union (an
	 *         empty geometry for empty input)
	 */
	public Geometry geometry(double maxSegmentLength) {
		return geometry(new GeometryFactory(), maxSegmentLength);
	}

	/**
	 * Computes (or reuses the cached) union topology and linearizes it.
	 *
	 * @param gf               geometry factory used to build the output
	 * @param maxSegmentLength maximum chord length for linearized arc segments;
	 *                         smaller values produce smoother output
	 * @return {@code Polygon} or {@code MultiPolygon} representing the union (an
	 *         empty geometry for empty input)
	 */
	public Geometry geometry(GeometryFactory gf, double maxSegmentLength) {
		Objects.requireNonNull(gf, "gf");
		if (!(maxSegmentLength > 0)) {
			throw new IllegalArgumentException("maxSegmentLength must be positive");
		}
		compute();
		return toGeometry(topology, gf, maxSegmentLength);
	}

	/**
	 * Computes (or reuses the cached) union topology and linearizes it using a
	 * default {@link GeometryFactory}, choosing the segment count per arc from its
	 * maximum deviation from the true circular arc (the sagitta) rather than a
	 * fixed chord length.
	 *
	 * @param maxDeviation maximum distance, in the same units as the input
	 *                     coordinates, between the linearized boundary and the true
	 *                     circular arc; smaller values produce smoother output
	 * @return {@code Polygon} or {@code MultiPolygon} representing the union (an
	 *         empty geometry for empty input)
	 */
	public Geometry geometryByDeviation(double maxDeviation) {
		return geometryByDeviation(new GeometryFactory(), maxDeviation);
	}

	/**
	 * Computes union topology and linearizes it, choosing the segment count per arc
	 * from its maximum deviation from the true circular arc (the sagitta) rather
	 * than a fixed chord length.
	 *
	 * @param gf           geometry factory used to build the output
	 * @param maxDeviation maximum distance, in the same units as the input
	 *                     coordinates, between the linearized boundary and the true
	 *                     circular arc; smaller values produce smoother output
	 * @return {@code Polygon} or {@code MultiPolygon} representing the union (an
	 *         empty geometry for empty input)
	 */
	public Geometry geometryByDeviation(GeometryFactory gf, double maxDeviation) {
		Objects.requireNonNull(gf, "gf");
		if (!(maxDeviation > 0)) {
			throw new IllegalArgumentException("maxDeviation must be positive");
		}
		compute();
		return toGeometryByDeviation(topology, gf, maxDeviation);
	}

	/**
	 * Computes the union of the given circles and linearizes the result.
	 *
	 * <p>
	 * Circles are specified as {@link Coordinate} objects where {@code x} and
	 * {@code y} define the center point and {@code z} defines the radius.
	 *
	 * @param circles          circles as Coordinates ({@code x,y} = center,
	 *                         {@code z} = radius)
	 * @param maxSegmentLength maximum chord length for linearized arc segments
	 * @return {@code Polygon} or {@code MultiPolygon} representing the union
	 * @throws IllegalArgumentException if any circle has NaN radius
	 */
	public static Geometry union(List<Coordinate> circles, double maxSegmentLength) {
		Objects.requireNonNull(circles, "circles");
		CircleUnion cu = new CircleUnion(circles.size());
		for (Coordinate c : circles) {
			double r = c.getZ();
			if (Double.isNaN(r)) {
				throw new IllegalArgumentException("circle.z (radius) must be set");
			}
			cu.add(c.x, c.y, r);
		}
		return cu.geometry(maxSegmentLength);
	}

	private void compute() {
		if (topology != null) {
			return;
		}
		final int n = pos;
		if (n == 0) {
			topology = Topology.EMPTY;
			return;
		}

		// ---- build: sort by radius descending so the larger circle of any pair
		// comes first ("owner" rule collapses to i < j) and any engulfer is
		// processed before what it engulfs.
		int[] order = new int[n];
		for (int i = 0; i < n; i++) {
			order[i] = i;
		}
		sortByRadiusDesc(order, inR);

		final double[] x = new double[n];
		final double[] y = new double[n];
		final double[] r = new double[n];
		for (int i = 0; i < n; i++) {
			int o = order[i];
			x[i] = inX[o];
			y[i] = inY[o];
			r[i] = inR[o];
		}
		final KDTree index = new KDTree(x, y);

		// ---- scan ----
		final boolean[] covered = new boolean[n];

		// Drop exact-duplicate circles up front. Must precede the pair sweep:
		// otherwise one larger circle crossing several duplicates mints
		// bit-identical intersection points under distinct IDs, desyncing the ring
		// handoff in stitch. First (largest-radius) circle at a location wins.
		// Only bit-identical centers need this; near coincidences mint distinct
		// points and are caught by the engulf test below.
		int hcap = Integer.highestOneBit(Math.max(2, n * 2 - 1)) << 1;
		int hmask = hcap - 1;
		int[] slot = new int[hcap];
		Arrays.fill(slot, -1);
		for (int i = 0; i < n; i++) {
			long hx = Double.doubleToLongBits(x[i]);
			long hy = Double.doubleToLongBits(y[i]);
			long h64 = hx * 0x9E3779B97F4A7C15L;
			h64 ^= hy + 0x9E3779B97F4A7C15L + (h64 << 6) + (h64 >>> 2);
			int h = (int) (h64 ^ (h64 >>> 32)) & hmask;
			while (slot[h] != -1 && (x[slot[h]] != x[i] || y[slot[h]] != y[i])) {
				h = (h + 1) & hmask;
			}
			if (slot[h] == -1) {
				slot[h] = i;
			} else {
				covered[i] = true;
			}
		}
		// Zero-radius circles are degenerate points; they contribute no area.
		for (int i = 0; i < n; i++) {
			if (r[i] == 0) {
				covered[i] = true;
			}
		}

		// disjoint-set forest over circles; each proper pair unions its endpoints,
		// so roots partition active circles into the union's connected components.
		final int[] parent = new int[n];
		final int[] setSize = new int[n];
		for (int i = 0; i < n; i++) {
			parent[i] = i;
			setSize[i] = 1;
		}

		// growable interleaved buffers: intersection points (x, y) and pairs
		// (i, j, baseId). A pair's two points are baseId (+h side) and baseId + 1.
		double[] points = new double[1 << 12];
		int[] pairs = new int[3 << 10];
		int np = 0;
		int pairCount = 0;
		final int[] neighbors = new int[n];

		for (int i = 0; i < n; i++) {
			if (covered[i]) {
				continue;
			}
			final double xi = x[i], yi = y[i], ri = r[i];
			// any intersecting j > i has rj <= ri, hence dist < ri + rj <= 2*ri
			final int m = index.within(xi, yi, 2 * ri, i, covered, neighbors);

			for (int t = 0; t < m; t++) {
				final int j = neighbors[t];
				final double dx = x[j] - xi;
				final double dy = y[j] - yi;
				final double d2 = dx * dx + dy * dy;
				final double rj = r[j];
				final double sumR = ri + rj;
				if (d2 >= sumR * sumR) {
					continue; // disjoint (externally tangent counts as disjoint)
				}
				final double diffR = ri - rj; // >= 0 by sort order
				if (d2 <= diffR * diffR) { // engulf: i contains j
					covered[j] = true;
					continue;
				}

				// proper intersection: standard two-circle formula, points p± on the
				// radical line at distance h either side of its foot.
				final double d = Math.sqrt(d2);
				final double a = (ri * ri - rj * rj + d2) / (2.0 * d);
				final double h2 = ri * ri - a * a;
				final double h = h2 > 0 ? Math.sqrt(h2) : 0.0; // roundoff guard
				final double invD = 1.0 / d;
				final double mx = xi + a * dx * invD;
				final double my = yi + a * dy * invD;
				final double ox = -dy * h * invD;
				final double oy = dx * h * invD;

				if (np * 2 + 4 > points.length) {
					points = Arrays.copyOf(points, points.length * 2);
				}
				if (pairCount * 3 + 3 > pairs.length) {
					pairs = Arrays.copyOf(pairs, pairs.length * 2);
				}

				final int baseId = np;
				final int p = np * 2;
				points[p] = mx + ox;
				points[p + 1] = my + oy;
				points[p + 2] = mx - ox;
				points[p + 3] = my - oy;
				np += 2;

				final int q = pairCount * 3;
				pairs[q] = i;
				pairs[q + 1] = j;
				pairs[q + 2] = baseId;
				pairCount++;

				dsuUnion(parent, setSize, i, j);
			}
		}

		// compact roots -> dense component ids over active circles
		final int[] component = new int[n];
		Arrays.fill(component, -1);
		int componentCount = 0;
		for (int i = 0; i < n; i++) {
			if (covered[i]) {
				continue;
			}
			int root = dsuFind(parent, i);
			if (component[root] == -1) {
				component[root] = componentCount++;
			}
			component[i] = component[root];
		}

		// ---- arcs: per-circle interval complement -> boundary arcs ----
		// Covered-interval events grouped per circle: theta, delta (+1 start / -1
		// end) and the intersection-point ID at that endpoint. Each proper pair
		// yields one interval per circle = two events; baseDepth[c] counts
		// intervals wrapping the atan2 seam at -pi.
		final int ne2 = pairCount * 4;
		final double[] evTheta = new double[ne2];
		final byte[] evDelta = new byte[ne2];
		final int[] evId = new int[ne2];
		final int[] evCircle = new int[ne2];
		final int[] baseDepth = new int[n];
		int ne = 0;

		for (int pi = 0; pi < pairCount; pi++) {
			final int q = pi * 3;
			final int i = pairs[q];
			final int j = pairs[q + 1];
			final int baseId = pairs[q + 2];
			final int pp = baseId * 2;
			final double ax = points[pp], ay = points[pp + 1];
			final double bx = points[pp + 2], by = points[pp + 3];
			ne = addInterval(evTheta, evDelta, evId, evCircle, baseDepth, ne, i, x[i], y[i], ax, ay, baseId, bx, by, baseId + 1, x[j] - x[i], y[j] - y[i]);
			ne = addInterval(evTheta, evDelta, evId, evCircle, baseDepth, ne, j, x[j], y[j], ax, ay, baseId, bx, by, baseId + 1, x[i] - x[j], y[i] - y[j]);
		}

		// counting sort event indices by circle -> contiguous per-circle ranges
		final int[] off = new int[n + 1];
		for (int k = 0; k < ne; k++) {
			off[evCircle[k] + 1]++;
		}
		for (int c = 0; c < n; c++) {
			off[c + 1] += off[c];
		}
		final int[] evOrder = new int[ne];
		final int[] cursor = Arrays.copyOf(off, n);
		for (int k = 0; k < ne; k++) {
			evOrder[cursor[evCircle[k]]++] = k;
		}

		// boundary arcs: circle index, CCW [thetaStart, thetaEnd], endpoint point
		// IDs. startId = -1 marks a full-circle arc. Upper bound: one arc per
		// covered interval (2 per pair) plus one full circle per circle.
		final int maxArcs = 2 * pairCount + n;
		final int[] arcCircle = new int[maxArcs];
		final double[] arcT0 = new double[maxArcs];
		final double[] arcT1 = new double[maxArcs];
		final int[] arcStartId = new int[maxArcs];
		final int[] arcEndId = new int[maxArcs];
		int arcCount = 0;

		for (int c = 0; c < n; c++) {
			if (covered[c]) {
				continue;
			}
			final int lo = off[c], hi = off[c + 1];

			if (lo == hi) { // active, no proper neighbors -> whole circle is one arc
				arcCircle[arcCount] = c;
				arcT0[arcCount] = 0;
				arcT1[arcCount] = TWO_PI;
				arcStartId[arcCount] = -1;
				arcEndId[arcCount] = -1;
				arcCount++;
				continue;
			}

			// insertion-sort this circle's event indices by theta (ties: +1 before -1)
			for (int a = lo + 1; a < hi; a++) {
				final int key = evOrder[a];
				final double kt = evTheta[key];
				final int kd = evDelta[key];
				int b = a - 1;
				while (b >= lo && (evTheta[evOrder[b]] > kt || (evTheta[evOrder[b]] == kt && evDelta[evOrder[b]] < kd))) {
					evOrder[b + 1] = evOrder[b];
					b--;
				}
				evOrder[b + 1] = key;
			}

			// sweep: boundary arcs are the runs where coverage depth is 0
			int depth = baseDepth[c];
			double gapTheta = 0, seamTheta = 0;
			int gapId = -1, seamId = -1;
			boolean haveGap = false, haveSeam = false;
			for (int a = lo; a < hi; a++) {
				final int k = evOrder[a];
				if (evDelta[k] == 1) {
					if (depth == 0) {
						if (haveGap) {
							arcCircle[arcCount] = c;
							arcT0[arcCount] = gapTheta;
							arcT1[arcCount] = evTheta[k];
							arcStartId[arcCount] = gapId;
							arcEndId[arcCount] = evId[k];
							arcCount++;
							haveGap = false;
						} else {
							seamTheta = evTheta[k];
							seamId = evId[k];
							haveSeam = true;
						}
					}
					depth++;
				} else {
					depth--;
					if (depth == 0) {
						gapTheta = evTheta[k];
						gapId = evId[k];
						haveGap = true;
					}
				}
			}
			if (haveGap && haveSeam) { // arc straddling the atan2 seam
				arcCircle[arcCount] = c;
				arcT0[arcCount] = gapTheta;
				arcT1[arcCount] = seamTheta;
				arcStartId[arcCount] = gapId;
				arcEndId[arcCount] = seamId;
				arcCount++;
			}
		}

		// planar arrangement bound: the union boundary of n disks has at most
		// 6n - 12 arcs (n >= 3). Exceeding it means the sweep produced spurious
		// arcs — an internal-consistency bug.
		int active = 0;
		for (int c = 0; c < n; c++) {
			if (!covered[c]) {
				active++;
			}
		}
		if (active >= 3 && arcCount > 6 * active - 12) {
			throw new IllegalStateException("Arc count exceeds the planar 6n-12 bound.");
		}

		// ---- stitch: walk end->start handoffs by point ID into closed rings and
		// nest them into polygons by connected component ----
		final int[] arcByStartId = new int[np];
		Arrays.fill(arcByStartId, -1);
		for (int k = 0; k < arcCount; k++) {
			if (arcStartId[k] != -1) {
				arcByStartId[arcStartId[k]] = k;
			}
		}

		final ArcRing[] shellOf = new ArcRing[componentCount];
		@SuppressWarnings("unchecked")
		final List<ArcRing>[] holesOf = new List[componentCount];

		final boolean[] visited = new boolean[arcCount];
		int consumed = 0;

		for (int k0 = 0; k0 < arcCount; k0++) {
			if (visited[k0]) {
				continue;
			}
			final int comp = component[arcCircle[k0]];

			if (arcStartId[k0] == -1) { // full circle -> standalone shell ring
				visited[k0] = true;
				consumed++;
				final int c = arcCircle[k0];
				final ArcRing ring = new ArcRing();
				ring.add(c, 0.0, TWO_PI);
				fileRing(shellOf, holesOf, ring, comp, Math.PI * r[c] * r[c]);
				continue;
			}

			final ArcRing ring = new ArcRing();
			double area = 0;
			double p0x = 0, p0y = 0;
			boolean havePoint0 = false;
			int k = k0;
			for (;;) {
				visited[k] = true;
				consumed++;

				final int sp = arcStartId[k] * 2;
				final int ep = arcEndId[k] * 2;
				final double ax = points[sp], ay = points[sp + 1];
				final double bx = points[ep], by = points[ep + 1];
				final int c = arcCircle[k];

				double dth = arcT1[k] - arcT0[k];
				if (dth < 0) {
					dth += TWO_PI;
				}
				ring.add(c, arcT0[k], arcT0[k] + dth);

				// signed ring area: circular segment between the CCW arc and its
				// chord, plus a shoelace fan through the arc endpoints
				final double rc = r[c];
				area += 0.5 * rc * rc * (dth - Math.sin(dth));
				if (!havePoint0) {
					p0x = ax;
					p0y = ay;
					havePoint0 = true;
				} else {
					area += 0.5 * ((ax - p0x) * (by - p0y) - (bx - p0x) * (ay - p0y));
				}

				final int next = arcByStartId[arcEndId[k]];
				if (next == k0) {
					break; // ring closed
				}
				if (next == -1 || visited[next]) {
					throw new IllegalStateException("Ring failed to close - arc handoff broke (degenerate input, e.g. three circles through one point).");
				}
				k = next;
			}
			fileRing(shellOf, holesOf, ring, comp, area);
		}

		if (consumed != arcCount) {
			throw new IllegalStateException("Not every arc was consumed exactly once while stitching rings.");
		}

		final ArcRing[][] polygons = new ArcRing[componentCount][];
		for (int comp = 0; comp < componentCount; comp++) {
			final ArcRing shell = shellOf[comp];
			if (shell == null) {
				throw new IllegalStateException("A connected component has no shell ring.");
			}
			final List<ArcRing> holes = holesOf[comp];
			final int hn = holes == null ? 0 : holes.size();
			final ArcRing[] poly = new ArcRing[1 + hn];
			poly[0] = shell;
			for (int h = 0; h < hn; h++) {
				poly[1 + h] = holes.get(h);
			}
			polygons[comp] = poly;
		}

		topology = new Topology(x, y, r, polygons);
	}

	/**
	 * Records circle {@code c}'s covered interval bounded by points A and B
	 * (selecting the arc toward the other circle's center) as two sweep events.
	 * Returns the updated event count.
	 */
	private static int addInterval(double[] evTheta, byte[] evDelta, int[] evId, int[] evCircle, int[] baseDepth, int ne, int c, double cx, double cy,
			double ax, double ay, int idA, double bx, double by, int idB, double dirX, double dirY) {
		final double thA = FastMath.atan2(ay - cy, ax - cx);
		final double thB = FastMath.atan2(by - cy, bx - cx);
		final double mid = FastMath.atan2(dirY, dirX);

		// CCW from thA to thB; the covered arc is the side containing mid
		double dab = thB - thA;
		if (dab < 0) {
			dab += TWO_PI;
		}
		double dam = mid - thA;
		if (dam < 0) {
			dam += TWO_PI;
		}
		final double s, e;
		final int sId, eId;
		if (dam <= dab) {
			s = thA;
			e = thB;
			sId = idA;
			eId = idB;
		} else {
			s = thB;
			e = thA;
			sId = idB;
			eId = idA;
		}

		if (s > e) {
			baseDepth[c]++; // interval wraps the seam
		}
		evTheta[ne] = s;
		evDelta[ne] = 1;
		evId[ne] = sId;
		evCircle[ne] = c;
		ne++;
		evTheta[ne] = e;
		evDelta[ne] = -1;
		evId[ne] = eId;
		evCircle[ne] = c;
		ne++;
		return ne;
	}

	private static void fileRing(ArcRing[] shellOf, List<ArcRing>[] holesOf, ArcRing ring, int comp, double area) {
		if (area >= 0) {
			if (shellOf[comp] != null) {
				throw new IllegalStateException("A connected component has more than one shell ring.");
			}
			shellOf[comp] = ring;
		} else {
			if (holesOf[comp] == null) {
				holesOf[comp] = new ArrayList<>();
			}
			holesOf[comp].add(ring);
		}
	}

	private static Geometry toGeometry(Topology topo, GeometryFactory gf, double maxSegLen) {
		final ArcRing[][] polygons = topo.polygons;
		if (polygons.length == 0) {
			return gf.createGeometryCollection();
		}
		final Polygon[] polys = new Polygon[polygons.length];
		for (int p = 0; p < polygons.length; p++) {
			final ArcRing[] rings = polygons[p];
			final LinearRing shell = toLinearRing(rings[0], topo, gf, maxSegLen);
			LinearRing[] holes = null;
			if (rings.length > 1) {
				holes = new LinearRing[rings.length - 1];
				for (int h = 1; h < rings.length; h++) {
					holes[h - 1] = toLinearRing(rings[h], topo, gf, maxSegLen);
				}
			}
			polys[p] = gf.createPolygon(shell, holes);
		}
		return polys.length == 1 ? polys[0] : gf.createMultiPolygon(polys);
	}

	private static LinearRing toLinearRing(ArcRing ring, Topology topo, GeometryFactory gf, double maxSegLen) {
		final int arcN = ring.count;
		// vertex floors keep every ring a valid LinearRing (>= 3 distinct vertices)
		final int minSegs = arcN == 1 ? 16 : (arcN == 2 ? 2 : 1);

		final int[] segs = new int[arcN];
		int total = 0;
		for (int a = 0; a < arcN; a++) {
			final double sweep = ring.t1[a] - ring.t0[a];
			int s = (int) Math.ceil((sweep * topo.r[ring.circle[a]]) / maxSegLen);
			if (s < minSegs) {
				s = minSegs;
			}
			segs[a] = s;
			total += s;
		}

		// each arc emits its start + interior samples and omits its end (the next
		// arc's shared start); the ring is closed by repeating its first vertex.
		final Coordinate[] coords = new Coordinate[total + 1];
		int idx = 0;
		for (int a = 0; a < arcN; a++) {
			final int c = ring.circle[a];
			final double cx = topo.x[c], cy = topo.y[c], rc = topo.r[c];
			final double t0 = ring.t0[a];
			final double dt = (ring.t1[a] - t0) / segs[a];
			for (int s = 0; s < segs[a]; s++) {
				final double th = t0 + s * dt;
				coords[idx++] = new Coordinate(cx + rc * Math.cos(th), cy + rc * Math.sin(th));
			}
		}
		coords[total] = new Coordinate(coords[0]);
		return gf.createLinearRing(coords);
	}

	private static Geometry toGeometryByDeviation(Topology topo, GeometryFactory gf, double maxDeviation) {
		final ArcRing[][] polygons = topo.polygons;
		if (polygons.length == 0) {
			return gf.createGeometryCollection();
		}
		final Polygon[] polys = new Polygon[polygons.length];
		for (int p = 0; p < polygons.length; p++) {
			final ArcRing[] rings = polygons[p];
			final LinearRing shell = toLinearRingByDeviation(rings[0], topo, gf, maxDeviation);
			LinearRing[] holes = null;
			if (rings.length > 1) {
				holes = new LinearRing[rings.length - 1];
				for (int h = 1; h < rings.length; h++) {
					holes[h - 1] = toLinearRingByDeviation(rings[h], topo, gf, maxDeviation);
				}
			}
			polys[p] = gf.createPolygon(shell, holes);
		}
		return polys.length == 1 ? polys[0] : gf.createMultiPolygon(polys);
	}

	/**
	 * Like {@link #toLinearRing}, but each arc's segment count is chosen from the
	 * maximum sagitta deviation (see {@link #maxAngleStep}) rather than a fixed
	 * chord length.
	 */
	private static LinearRing toLinearRingByDeviation(ArcRing ring, Topology topo, GeometryFactory gf, double maxDeviation) {
		final int arcN = ring.count;
		// vertex floors keep every ring a valid LinearRing (>= 3 distinct vertices)
		final int minSegs = arcN == 1 ? 16 : (arcN == 2 ? 2 : 1);

		final int[] segs = new int[arcN];
		int total = 0;
		for (int a = 0; a < arcN; a++) {
			final double sweep = ring.t1[a] - ring.t0[a];
			final double maxStep = maxAngleStep(topo.r[ring.circle[a]], maxDeviation);
			int s = (int) Math.ceil(sweep / maxStep);
			if (s < minSegs) {
				s = minSegs;
			}
			segs[a] = s;
			total += s;
		}

		// each arc emits its start + interior samples and omits its end (the next
		// arc's shared start); the ring is closed by repeating its first vertex.
		final Coordinate[] coords = new Coordinate[total + 1];
		int idx = 0;
		for (int a = 0; a < arcN; a++) {
			final int c = ring.circle[a];
			final double cx = topo.x[c], cy = topo.y[c], rc = topo.r[c];
			final double t0 = ring.t0[a];
			final double dt = (ring.t1[a] - t0) / segs[a];
			for (int s = 0; s < segs[a]; s++) {
				final double th = t0 + s * dt;
				coords[idx++] = new Coordinate(cx + rc * Math.cos(th), cy + rc * Math.sin(th));
			}
		}
		coords[total] = new Coordinate(coords[0]);
		return gf.createLinearRing(coords);
	}

	/**
	 * Largest angular step (radians) along a circle of radius {@code r} whose
	 * chord's sagitta does not exceed {@code maxDeviation}, i.e. the same
	 * {@code acos(1 - maxDeviation / r)} relationship used by circle linearizers
	 * such as {@code createCircle}. Degenerate ({@code r <= 0}) or oversized
	 * ({@code maxDeviation >= r}) inputs collapse to a single half-turn step since
	 * the whole arc is then within tolerance.
	 */
	private static double maxAngleStep(double r, double maxDeviation) {
		if (r <= 0 || maxDeviation >= r) {
			return Math.PI;
		}
		return Math.acos(1 - maxDeviation / r);
	}

	/** One closed boundary ring as arcs (circle index + unwrapped CCW angles). */
	private static final class ArcRing {
		int count;
		int[] circle = new int[4];
		double[] t0 = new double[4];
		double[] t1 = new double[4];

		void add(int c, double a0, double a1) {
			if (count == circle.length) {
				final int cap = count * 2;
				circle = Arrays.copyOf(circle, cap);
				t0 = Arrays.copyOf(t0, cap);
				t1 = Arrays.copyOf(t1, cap);
			}
			circle[count] = c;
			t0[count] = a0;
			t1[count] = a1;
			count++;
		}
	}

	/** Cached union topology: sorted circle data plus rings per component. */
	private static final class Topology {
		static final Topology EMPTY = new Topology(new double[0], new double[0], new double[0], new ArcRing[0][]);

		final double[] x, y, r;
		/** Per component: index 0 is the shell, the rest are holes. */
		final ArcRing[][] polygons;

		Topology(double[] x, double[] y, double[] r, ArcRing[][] polygons) {
			this.x = x;
			this.y = y;
			this.r = r;
			this.polygons = polygons;
		}
	}

	private static int dsuFind(int[] parent, int i) {
		while (parent[i] != i) {
			parent[i] = parent[parent[i]];
			i = parent[i];
		}
		return i;
	}

	private static void dsuUnion(int[] parent, int[] setSize, int a, int b) {
		int ra = dsuFind(parent, a);
		int rb = dsuFind(parent, b);
		if (ra == rb) {
			return;
		}
		if (setSize[ra] < setSize[rb]) {
			final int t = ra;
			ra = rb;
			rb = t;
		}
		parent[rb] = ra;
		setSize[ra] += setSize[rb];
	}

	/**
	 * Sorts circle indices by descending radius (ties by ascending index, for
	 * determinism). Iterative quicksort over an index array.
	 */
	private static void sortByRadiusDesc(int[] order, double[] r) {
		int[] stack = new int[64];
		int sp = 0;
		int left = 0, right = order.length - 1;

		while (true) {
			while (left < right) {
				final int pivot = order[(left + right) >>> 1];
				final double pivotR = r[pivot];
				int i = left - 1;
				int j = right + 1;

				while (true) {
					int oi, oj;
					double ri, rj;
					do {
						i++;
						oi = order[i];
						ri = r[oi];
					} while (ri > pivotR || (ri == pivotR && oi < pivot));
					do {
						j--;
						oj = order[j];
						rj = r[oj];
					} while (pivotR > rj || (pivotR == rj && pivot < oj));
					if (i >= j) {
						break;
					}
					final int tmp = order[i];
					order[i] = order[j];
					order[j] = tmp;
				}

				if (j - left < right - j - 1) {
					if (left < j) {
						if (sp + 2 > stack.length) {
							stack = Arrays.copyOf(stack, stack.length * 2);
						}
						stack[sp++] = left;
						stack[sp++] = j;
					}
					left = j + 1;
				} else {
					if (j + 1 < right) {
						if (sp + 2 > stack.length) {
							stack = Arrays.copyOf(stack, stack.length * 2);
						}
						stack[sp++] = j + 1;
						stack[sp++] = right;
					}
					right = j;
				}
			}

			if (sp == 0) {
				break;
			}
			right = stack[--sp];
			left = stack[--sp];
		}
	}

	private static final class KDTree {
		private static final int NODE_SIZE = 16;

		private final double[] xs, ys;
		private final int[] ids;
		private final int n;
		private int[] stack = new int[256];

		KDTree(double[] x, double[] y) {
			n = x.length;
			xs = x.clone();
			ys = y.clone();
			ids = new int[n];
			for (int i = 0; i < n; i++) {
				ids[i] = i;
			}
			sort(0, n - 1, 0);
		}

		/**
		 * Collects ids of points within {@code rad} of {@code (qx, qy)} that pass the
		 * sweep filter ({@code id > owner} and not covered) into {@code out}. Returns
		 * the number collected.
		 */
		int within(double qx, double qy, double rad, int owner, boolean[] covered, int[] out) {
			final double r2 = rad * rad;
			int count = 0;
			int[] st = stack;
			int sp = 0;
			int left = 0, right = n - 1, axis = 0;

			for (;;) {
				if (right - left <= NODE_SIZE) {
					for (int i = left; i <= right; i++) {
						final int id = ids[i];
						if (id > owner && !covered[id]) {
							final double dx = xs[i] - qx;
							final double dy = ys[i] - qy;
							if (dx * dx + dy * dy <= r2) {
								out[count++] = id;
							}
						}
					}
				} else {
					final int m = (left + right) >>> 1;
					final int id = ids[m];
					if (id > owner && !covered[id]) {
						final double dx = xs[m] - qx;
						final double dy = ys[m] - qy;
						if (dx * dx + dy * dy <= r2) {
							out[count++] = id;
						}
					}
					final double c = axis == 0 ? xs[m] : ys[m];
					final double q = axis == 0 ? qx : qy;
					final int nextAxis = 1 - axis;
					final boolean goLeft = q - rad <= c;
					final boolean goRight = q + rad >= c;
					if (goLeft) {
						if (goRight) {
							if (sp + 3 > st.length) {
								st = Arrays.copyOf(st, st.length * 2);
								stack = st;
							}
							st[sp++] = m + 1;
							st[sp++] = right;
							st[sp++] = nextAxis;
						}
						right = m - 1;
						axis = nextAxis;
						continue;
					} else if (goRight) {
						left = m + 1;
						axis = nextAxis;
						continue;
					}
				}
				if (sp == 0) {
					break;
				}
				axis = st[--sp];
				right = st[--sp];
				left = st[--sp];
			}
			return count;
		}

		private void sort(int left, int right, int axis) {
			if (right - left <= NODE_SIZE) {
				return;
			}
			final int m = (left + right) >>> 1;
			select(m, left, right, axis);
			sort(left, m - 1, 1 - axis);
			sort(m + 1, right, 1 - axis);
		}

		// Floyd-Rivest selection of the k-th item along the given axis
		private void select(int k, int left, int right, int axis) {
			final double[] coords = axis == 0 ? xs : ys;
			while (right > left) {
				if (right - left > 600) {
					final int nn = right - left + 1;
					final int mm = k - left + 1;
					final double z = Math.log(nn);
					final double s = 0.5 * Math.exp(2.0 * z / 3.0);
					final double sd = 0.5 * Math.sqrt(z * s * (nn - s) / nn) * (mm - nn / 2.0 < 0 ? -1 : 1);
					final int newLeft = Math.max(left, (int) Math.floor(k - (double) mm * s / nn + sd));
					final int newRight = Math.min(right, (int) Math.floor(k + (double) (nn - mm) * s / nn + sd));
					select(k, newLeft, newRight, axis);
				}
				final double t = coords[k];
				int i = left;
				int j = right;
				swap(left, k);
				if (coords[right] > t) {
					swap(left, right);
				}
				while (i < j) {
					swap(i, j);
					i++;
					j--;
					while (coords[i] < t) {
						i++;
					}
					while (coords[j] > t) {
						j--;
					}
				}
				if (coords[left] == t) {
					swap(left, j);
				} else {
					j++;
					swap(j, right);
				}
				if (j <= k) {
					left = j + 1;
				}
				if (k <= j) {
					right = j - 1;
				}
			}
		}

		private void swap(int i, int j) {
			final double tx = xs[i];
			xs[i] = xs[j];
			xs[j] = tx;
			final double ty = ys[i];
			ys[i] = ys[j];
			ys[j] = ty;
			final int ti = ids[i];
			ids[i] = ids[j];
			ids[j] = ti;
		}
	}
}
