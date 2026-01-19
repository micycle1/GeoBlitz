package com.github.micycle1.geoblitz;

import java.util.Arrays;
import java.util.Comparator;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;

// slower than JTS ConvexHull on inputs of ~1M+.
/**
 * Computes the convex hull of a set of points using a fast monotone chain
 * algorithm.
 */
public final class FastConvexHull {

	/*
	 * A Java port of Jernej Puc's "optimised version" (Python) of A. M. Andrew's
	 * monotone chain algorithm for constructing convex hulls. Unlike the Python
	 * version, this version handles duplicate coordinates (that are otherwise a
	 * degenerate case).
	 */

	private GeometryFactory geomFactory;
	private Coordinate[] inputPts;

	/**
	 * Creates a new FastConvexHull for the given Geometry's coordinates.
	 *
	 * @param geometry the geometry whose points will be used
	 */
	public FastConvexHull(Geometry geometry) {
		this(geometry.getCoordinates(), geometry.getFactory());
	}

	/**
	 * Creates a new FastConvexHull for the given points.
	 *
	 * @param pts         the input points
	 * @param geomFactory the geometry factory to use for the result
	 */
	public FastConvexHull(Coordinate[] pts, GeometryFactory geomFactory) {
		inputPts = pts;

		this.geomFactory = geomFactory;
	}

	/**
	 * Computes the convex hull.
	 *
	 * @return the convex hull geometry
	 */
	public Geometry getConvexHull() {
		if (inputPts.length == 0) {
			return geomFactory.createGeometryCollection(new Geometry[0]);
		}
		final Coordinate[] hull = computeHullCoordinates(inputPts);
		final int m = hull.length;

		if (m == 0) {
			return geomFactory.createGeometryCollection(new Geometry[0]);
		} else if (m == 1) {
			final Point p = geomFactory.createPoint(hull[0]);
			return p;
		} else if (m == 2) {
			final LineString ls = geomFactory.createLineString(hull);
			return ls;
		} else {
			// Close the ring
			final Coordinate[] ring = new Coordinate[m + 1];
			System.arraycopy(hull, 0, ring, 0, m);
			ring[m] = new Coordinate(hull[0]); // close; create new to avoid aliasing
			final LinearRing shell = geomFactory.createLinearRing(ring);
			final Polygon poly = geomFactory.createPolygon(shell, null);
			return poly;
		}
	}

	private static Coordinate[] computeHullCoordinates(final Coordinate[] pts) {
		final int n0 = pts.length;
		if (n0 == 0) {
			return new Coordinate[0];
		}
		if (n0 == 1) {
			return new Coordinate[] { pts[0] };
		}

		// Sort and deduplicate
		final int n = sortAndUnique(pts);
		if (n == 1) {
			return new Coordinate[] { pts[0] };
		}
		if (n == 2) {
			return new Coordinate[] { pts[0], pts[1] };
		}

		final Coordinate p0 = pts[0];
		final Coordinate pn = pts[n - 1];

		// Preallocate stacks (arrays) for upper and lower
		final Coordinate[] upper = new Coordinate[n];
		final Coordinate[] lower = new Coordinate[n];
		int uLen = 0;
		int lLen = 0;

		// Add left endpoint to both
		upper[uLen++] = p0;
		lower[lLen++] = p0;

		// Build middle of upper and lower sections, using the dynamic
		// separation line anchored at the last element of each stack and pointing to
		// pn.
		for (int j = 1; j < n - 1; j++) {
			final Coordinate p = pts[j];

			// Upper section: point must be "above" current sep line
			if (isAboveSep(upper[uLen - 1], pn, p)) {
				while (uLen > 1 && isLeftTurn(upper[uLen - 2], upper[uLen - 1], p)) {
					uLen--;
				}
				upper[uLen++] = p; // sep anchor implicitly updated to p
			}
			// Lower section: point must be "below" current sep line
			else if (isBelowSep(lower[lLen - 1], pn, p)) {
				while (lLen > 1 && !isLeftTurn(lower[lLen - 2], lower[lLen - 1], p)) {
					lLen--;
				}
				lower[lLen++] = p; // sep anchor implicitly updated to p
			}
			// points between the current separation lines are ignored
		}

		// Finalize with the right endpoint
		while (uLen > 1 && isLeftTurn(upper[uLen - 2], upper[uLen - 1], pn)) {
			uLen--;
		}
		while (lLen > 1 && !isLeftTurn(lower[lLen - 2], lower[lLen - 1], pn)) {
			lLen--;
		}
		upper[uLen++] = pn;

		// Merge: upper (left->right) + reversed(lower) excluding the duplicate p0
		final int hullSize = uLen + lLen - 1; // unique vertices count
		final Coordinate[] hull = new Coordinate[hullSize];
		System.arraycopy(upper, 0, hull, 0, uLen);

		int k = uLen;
		for (int i = lLen - 1; i >= 1; i--) { // skip lower[0] == p0
			hull[k++] = lower[i];
		}

		return hull;
	}

	// Cross product (b - a) x (p - a)
	// > 0 => p is left of a->b (CCW turn), < 0 => right, = 0 => collinear
	private static double cross(final Coordinate a, final Coordinate b, final Coordinate p) {
		final double bax = b.x - a.x;
		final double bay = b.y - a.y;
		final double pax = p.x - a.x;
		final double pay = p.y - a.y;
		return bax * pay - bay * pax;
	}

	private static boolean isLeftTurn(final Coordinate p, final Coordinate q, final Coordinate r) {
		return cross(p, q, r) > 0.0;
	}

	// Returns true if p is "above" the separation line from sepAnchor to pn
	// i.e. p is to the left of the directed segment sepAnchor -> pn
	private static boolean isAboveSep(final Coordinate sepAnchor, final Coordinate pn, final Coordinate p) {
		return cross(sepAnchor, pn, p) > 0.0;
	}

	// Returns true if p is "below" the separation line from sepAnchor to pn
	// i.e. p is to the right of the directed segment sepAnchor -> pn
	private static boolean isBelowSep(final Coordinate sepAnchor, final Coordinate pn, final Coordinate p) {
		return cross(sepAnchor, pn, p) < 0.0;
	}

	// In-place sort by (x, y) and deduplicate consecutive equal points.
	// Returns the number of unique points remaining in a[0..m-1].
	private static int sortAndUnique(final Coordinate[] a) {
		if (a.length <= 1) {
			return a.length;
		}

		Arrays.sort(a, XY_ORDER);
		int m = 1;
		Coordinate prev = a[0];
		for (int i = 1; i < a.length; i++) {
			final Coordinate c = a[i];
			if (c.x != prev.x || c.y != prev.y) {
				a[m++] = c;
				prev = c;
			}
		}
		return m;
	}

	// Static comparator to avoid lambda allocation each call
	private static final Comparator<Coordinate> XY_ORDER = new Comparator<>() {
		@Override
		public int compare(Coordinate a, Coordinate b) {
			if (a.x < b.x) {
				return -1;
			}
			if (a.x > b.x) {
				return 1;
			}
			// tie-breaker on y
			if (a.y < b.y) {
				return -1;
			}
			if (a.y > b.y) {
				return 1;
			}
			return 0;
		}
	};
}