package com.github.micycle1.geoblitz;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.util.GeometryCombiner;

public class EndpointSnapperTest {

	private static final GeometryFactory GF = new GeometryFactory();
	private static final double EPS = 1e-9;

	// Helpers
	private static Coordinate c(double x, double y) {
		return new Coordinate(x, y);
	}

	private static void assertCoordEquals(Coordinate expected, Coordinate actual, double eps) {
		assertNotNull(actual, "actual coordinate is null");
		assertEquals(expected.x, actual.x, eps, "x differs");
		assertEquals(expected.y, actual.y, eps, "y differs");
	}

	private static LineString ls(Coordinate... coords) {
		return GF.createLineString(coords);
	}

	private static Polygon poly(Coordinate... ring) {
		// Ensure first = last
		if (!ring[0].equals2D(ring[ring.length - 1])) {
			Coordinate[] closed = new Coordinate[ring.length + 1];
			System.arraycopy(ring, 0, closed, 0, ring.length);
			closed[closed.length - 1] = new Coordinate(ring[0]);
			ring = closed;
		}
		return GF.createPolygon(ring);
	}

	@Test
	void noChangeWhenToleranceTooSmall() {
		LineString a = ls(c(0, 0), c(1, 0));
		LineString b = ls(c(1.01, 0), c(2, 0)); // gap = 0.01
		MultiLineString in = GF.createMultiLineString(new LineString[] { a, b });

		EndpointSnapper snapper = new EndpointSnapper(0.005); // smaller than gap
		Geometry out = snapper.snapEndpoints(in, false);

		MultiLineString ml = (MultiLineString) out;
		assertCoordEquals(c(1, 0), ((LineString) ml.getGeometryN(0)).getCoordinateN(1), EPS);
		assertCoordEquals(c(1.01, 0), ((LineString) ml.getGeometryN(1)).getCoordinateN(0), EPS);
	}

	@Test
	void mutualSnappingToMeanWithoutPolygons() {
		LineString a = ls(c(0, 0), c(1, 0));
		LineString b = ls(c(1.01, 0), c(2, 0)); // gap = 0.01
		MultiLineString in = GF.createMultiLineString(new LineString[] { a, b });

		EndpointSnapper snapper = new EndpointSnapper(0.02);
		Geometry out = snapper.snapEndpoints(in, false);

		MultiLineString ml = (MultiLineString) out;
		Coordinate endA = ((LineString) ml.getGeometryN(0)).getCoordinateN(1);
		Coordinate startB = ((LineString) ml.getGeometryN(1)).getCoordinateN(0);
		// Mean of 1.00 and 1.01 = 1.005
		assertCoordEquals(c(1.005, 0), endA, 1e-12);
		assertCoordEquals(c(1.005, 0), startB, 1e-12);
	}

	@Test
	void anchoredToPolygonVertexWhenEnabled() {
		Polygon face = poly(c(1, 0), c(2, 0), c(2, 1), c(1, 1)); // vertex at (1,0)
		LineString dangling = ls(c(0, 0), c(1.008, 0.002)); // near (1,0)
		Geometry in = GeometryCombiner.combine(face, dangling);

		EndpointSnapper snapper = new EndpointSnapper(0.02);
		Geometry out = snapper.snapEndpoints(in, true);

		// Polygon should remain unchanged
		Polygon outPoly = (Polygon) ((GeometryCollection) out).getGeometryN(0);
		assertTrue(outPoly.equalsExact(face), "Polygon must be unchanged");

		LineString outLine = (LineString) ((GeometryCollection) out).getGeometryN(1);
		// End should snap exactly to polygon vertex (1,0)
		assertCoordEquals(c(1, 0), outLine.getCoordinateN(outLine.getNumPoints() - 1), 1e-12);
		// Start unchanged
		assertCoordEquals(c(0, 0), outLine.getCoordinateN(0), 1e-12);
	}

	@Test
	void polygonsIgnoredWhenSnapToValidFalse() {
		Polygon face = poly(c(1, 0), c(2, 0), c(2, 1), c(1, 1)); // vertex at (1,0)
		LineString dangling = ls(c(0, 0), c(1.008, 0.002)); // near (1,0)
		Geometry in = GeometryCombiner.combine(face, dangling);

		EndpointSnapper snapper = new EndpointSnapper(0.02);
		Geometry out = snapper.snapEndpoints(in, false);

		// Polygon unchanged
		Polygon outPoly = (Polygon) ((GeometryCollection) out).getGeometryN(0);
		assertTrue(outPoly.equalsExact(face), "Polygon must be unchanged");

		// No other endpoints nearby => no snap should occur (polygons ignored)
		LineString outLine = (LineString) ((GeometryCollection) out).getGeometryN(1);
		assertCoordEquals(c(1.008, 0.002), outLine.getCoordinateN(outLine.getNumPoints() - 1), 1e-12);
	}

	@Test
	void closedLineStringIsIgnoredAsLine() {
		// Closed LineString (ring) not provided as Polygon
		LineString ring = ls(c(0, 0), c(1, 0), c(1, 1), c(0, 0)); // closed
		LineString near = ls(c(2, 0), c(3, 0.01));

		MultiLineString in = GF.createMultiLineString(new LineString[] { ring, near });

		EndpointSnapper snapper = new EndpointSnapper(0.05);
		Geometry out = snapper.snapEndpoints(in, false);

		// Ring should be returned unchanged (not treated as a movable line)
		LineString ringOut = (LineString) ((MultiLineString) out).getGeometryN(0);
		assertTrue(ringOut.equalsExact(ring), "Closed LineString must be unchanged");
	}

	@Test
	void transitiveClusteringChainsEndpoints() {
		// E1=10.000, E2=10.015, E3=10.030 (E1-E2 within tol, E2-E3 within tol, E1-E3
		// outside)
		LineString l1 = ls(c(9, 0), c(10.000, 0));
		LineString l2 = ls(c(10.015, 0), c(11, 0));
		LineString l3 = ls(c(10.030, 0), c(12, 0));
		MultiLineString in = GF.createMultiLineString(new LineString[] { l1, l2, l3 });

		EndpointSnapper snapper = new EndpointSnapper(0.02);
		Geometry out = snapper.snapEndpoints(in, false);

		MultiLineString ml = (MultiLineString) out;
		Coordinate e1 = ((LineString) ml.getGeometryN(0)).getCoordinateN(1);
		Coordinate s2 = ((LineString) ml.getGeometryN(1)).getCoordinateN(0);
		Coordinate s3 = ((LineString) ml.getGeometryN(2)).getCoordinateN(0);

		// Mean = (10.000 + 10.015 + 10.030) / 3 = 10.015
		Coordinate mean = c(10.015, 0);
		assertCoordEquals(mean, e1, 1e-12);
		assertCoordEquals(mean, s2, 1e-12);
		assertCoordEquals(mean, s3, 1e-12);
	}

	@Test
	void endpointOnlyInternalVerticesUnchanged() {
		LineString l = ls(c(0.99, 1.02), c(0.5, 1.5), c(0.0, 2.0)); // middle vertex should remain unchanged
		Polygon anchor = poly(c(1, 1), c(2, 1), c(2, 2), c(1, 2)); // vertex at (1,1), near start

		Geometry in = GeometryCombiner.combine(anchor, l);
		EndpointSnapper snapper = new EndpointSnapper(0.05);
		Geometry out = snapper.snapEndpoints(in, true);

		LineString outLine = (LineString) ((GeometryCollection) out).getGeometryN(1);
		// Start should snap to (1,1)
		assertCoordEquals(c(1, 1), outLine.getCoordinateN(0), 1e-12);
		// Middle and end unchanged
		assertCoordEquals(c(0.5, 1.5), outLine.getCoordinateN(1), 1e-12);
		assertCoordEquals(c(0.0, 2.0), outLine.getCoordinateN(2), 1e-12);
	}

	@Test
	void anchorsFromHoleVerticesAlsoApply() {
		// Polygon with a hole; anchor should include hole ring vertices
		LinearRing shell = GF.createLinearRing(new Coordinate[] { c(0, 0), c(5, 0), c(5, 5), c(0, 5), c(0, 0) });
		LinearRing hole = GF.createLinearRing(new Coordinate[] { c(2, 2), c(3, 2), c(3, 3), c(2, 3), c(2, 2) // hole vertex at (2,2)
		});
		Polygon withHole = GF.createPolygon(shell, new LinearRing[] { hole });

		LineString nearHoleVertex = ls(c(1.99, 2.01), c(1, 1)); // near (2,2)
		Geometry in = GeometryCombiner.combine(withHole, nearHoleVertex);

		EndpointSnapper snapper = new EndpointSnapper(0.05);
		Geometry out = snapper.snapEndpoints(in, true);

		LineString outLine = (LineString) ((GeometryCollection) out).getGeometryN(1);
		assertCoordEquals(c(2, 2), outLine.getCoordinateN(0), 1e-12);
	}

	@Test
	void mixedGeometryCollectionStructurePreserved() {
		Polygon face = poly(c(0, 0), c(2, 0), c(2, 2), c(0, 2));
		LineString a = ls(c(2.01, 1), c(3, 1));
		Geometry input = GF.createGeometryCollection(new Geometry[] { face, a });

		EndpointSnapper snapper = new EndpointSnapper(0.02);
		Geometry out = snapper.snapEndpoints(input, true);

		assertTrue(out instanceof GeometryCollection);
		GeometryCollection gc = (GeometryCollection) out;

		assertTrue(gc.getGeometryN(0) instanceof Polygon);
		assertTrue(gc.getGeometryN(1) instanceof LineString);

		// Polygon unchanged
		assertTrue(gc.getGeometryN(0).equalsExact(face));

		// Line endpoint snapped to polygon vertex (2,1) does not exist; nearest vertex
		// is one of the shell corners.
		// The start point is near edge, not a vertex; unless within tolerance to a
		// vertex, it should remain if no endpoint-near-vertex.
		// Here start (2.01,1) is 0.01 from x=2 line, but vertex distances > 1.0, so no
		// snap should occur.
		LineString outLine = (LineString) gc.getGeometryN(1);
		assertCoordEquals(c(2.01, 1), outLine.getCoordinateN(0), 1e-12);
		assertCoordEquals(c(3.00, 1), outLine.getCoordinateN(1), 1e-12);
	}
}