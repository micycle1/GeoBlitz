package com.github.micycle1.geoblitz;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * @author Michael Carleton
 */
public class QuickConvexHull {

	private GeometryFactory geomFactory;
	private Coordinate[] inputPts;

	public QuickConvexHull(Geometry geometry) {
		this(geometry.getCoordinates(), geometry.getFactory());
	}

	public QuickConvexHull(Coordinate[] pts, GeometryFactory geomFactory) {
		inputPts = pts;
		this.geomFactory = geomFactory;
	}

	public Geometry getConvexHull() {
		var out = convexHull(inputPts);
		return lineOrPolygon(out);
	}

	/*
	 * A Java port of Jernej Puc's "optimised version" (Python) of A. M. Andrew's
	 * monotone chain algorithm for constructing convex hulls. Unlike the Python
	 * version, this version handles duplicate coordinates (that are otherwise a
	 * degenerate case).
	 */

//	public static Geometry get(Coordinate[] pts, GeometryFactory geomFactory) {
//		// -- suboptimal early uniquing - for performance testing only
//		// inputPts = UniqueCoordinateArrayFilter.filterCoordinates(pts);
//
//
//	}

	private static Coordinate[] convexHull(Coordinate[] P) {
		// Preprocess: Convert array to list for sorting
		List<Coordinate> points = new ArrayList<>();
		for (Coordinate coord : P) {
			points.add(coord);
		}

		// Sort the list of coordinates
		points.sort(new CoordinateComparator());

		// Data structures
		List<Coordinate> upper = new ArrayList<>(P.length / 2);
		List<Coordinate> lower = new ArrayList<>(P.length / 2);

		// Endpoints
		Coordinate p0 = points.get(0);
		Coordinate pn = points.get(points.size() - 1);

		// Initial line(s) of separation
		double ku = (pn.y - p0.y) / (pn.x - p0.x + 1e-12);
		double nu = p0.y - ku * p0.x;
		double kl = ku;
		double nl = nu;

		// Add left endpoint
		upper.add(p0);
		lower.add(p0);

		// Construct the middle of the upper and lower hull sections
		for (int j = 1; j < points.size() - 1; j++) {
			Coordinate p = points.get(j);

			if (p.y > ku * p.x + nu) {
				while (upper.size() > 1 && isLeftTurn(upper.get(upper.size() - 2), upper.get(upper.size() - 1), p)) {
					upper.remove(upper.size() - 1);
				}
				upper.add(p);

				// Update upper line of separation
				ku = (pn.y - p.y) / (pn.x - p.x + 1e-12);
				nu = p.y - ku * p.x;

			} else if (p.y < kl * p.x + nl) {
				while (lower.size() > 1 && !isLeftTurn(lower.get(lower.size() - 2), lower.get(lower.size() - 1), p)) {
					lower.remove(lower.size() - 1);
				}

				lower.add(p);

				// Update lower line of separation
				kl = (pn.y - p.y) / (pn.x - p.x + 1e-12);
				nl = p.y - kl * p.x;
			}
		}

		// Add right endpoint (only once due to merging that follows)
		while (upper.size() > 1 && isLeftTurn(upper.get(upper.size() - 2), upper.get(upper.size() - 1), pn)) {
			upper.remove(upper.size() - 1);
		}

		while (lower.size() > 1 && !isLeftTurn(lower.get(lower.size() - 2), lower.get(lower.size() - 1), pn)) {
			lower.remove(lower.size() - 1);
		}

		upper.add(pn);

		// Reverse lower hull section
		Collections.reverse(lower);

		// Merge hull sections
		upper.addAll(lower);

		// Convert the list of coordinates back to an array
		Coordinate[] hull = new Coordinate[upper.size()];
		return upper.toArray(hull);
	}

	/**
	 * @param vertices the vertices of a linear ring, which may or may not be
	 *                 flattened (i.e. vertices collinear)
	 * @return a 2-vertex <code>LineString</code> if the vertices are collinear;
	 *         otherwise, a <code>Polygon</code> with unnecessary (collinear)
	 *         vertices removed
	 */
	private Geometry lineOrPolygon(Coordinate[] coordinates) {

//	    coordinates = cleanRing(coordinates);
		if (coordinates.length == 3) {
			return geomFactory.createLineString(new Coordinate[] { coordinates[0], coordinates[1] });
		}
		LinearRing linearRing = geomFactory.createLinearRing(coordinates);
		return geomFactory.createPolygon(linearRing);
	}

	private static class CoordinateComparator implements Comparator<Coordinate> {
		@Override
		public int compare(final Coordinate p1, final Coordinate p2) {
			final int xComparison = Double.compare(p1.x, p2.x);
			if (xComparison == 0) {
				return Double.compare(p1.y, p2.y);
			}
			return xComparison;
		}
	}

	/**
	 * Returns True if the turn of pq -> qr is counter-clockwise and False
	 * otherwise.
	 */
	private static boolean isLeftTurn(final Coordinate p, final Coordinate q, final Coordinate r) {
		if (q.x == r.x && q.y == r.y) { // q==r
			return true; // handle duplicate point (degenerate case) -- return true to remove it
		}
		return (q.x - p.x) * (r.y - p.y) > (r.x - p.x) * (q.y - p.y);
	}

}