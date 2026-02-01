package com.github.micycle1.geoblitz;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.operation.distance.IndexedFacetDistance;

class PointDistanceIndexTest {

	private static final double TOLERANCE = 1e-9;
	private final GeometryFactory gf = new GeometryFactory();

	@Test
	void testAgainstJTSIndexedFacetDistance() {
		long seed = 12345L;
		// Create a random complex polygon
		Geometry polygon = GeomMaker.make(100, seed);

		PointDistanceIndex ild = new PointDistanceIndex(polygon);
		IndexedFacetDistance ifd = new IndexedFacetDistance(polygon);

		// Random query points
		var rnd = new Random(seed);
		for (int i = 0; i < 1000; i++) {
			double x = rnd.nextDouble() * 1200 - 100; // range slightly larger than typical 0-1000 GeomMaker
			double y = rnd.nextDouble() * 1200 - 100;
			Coordinate c = new Coordinate(x, y);
			Point p = gf.createPoint(c);

			double distJTS = ifd.distance(p);
			double distILD = ild.unsignedDistance(c);

			assertEquals(distJTS, distILD, TOLERANCE, "Distance mismatch at " + c);
		}
	}

	@Test
	void testSignedDistance() {
		Geometry poly = GeomMaker.make(50, 1337); // Generates a polygon
		PointDistanceIndex ild = new PointDistanceIndex(poly);

		// Internal point
		Point centroid = poly.getCentroid();
		// Ensure centroid is actually inside (it might not be for concave polys, but
		// GeomMaker usually makes star-ish shapes or we can checking contains)
		if (poly.contains(centroid)) {
			double d = ild.distance(centroid.getCoordinate());
			assertTrue(d >= 0, "Interior point should have positive distance");
		}

		// External point (far away)
		Coordinate outside = new Coordinate(-10000, -10000);
		double dOut = ild.distance(outside);
		assertTrue(dOut < 0, "Exterior point should have negative distance");

		// Boundary point
		Coordinate boundaryPt = poly.getCoordinates()[0];
		double dBound = ild.distance(boundaryPt);
		assertEquals(0, dBound, TOLERANCE, "Boundary point should be 0 distance");
	}
}
