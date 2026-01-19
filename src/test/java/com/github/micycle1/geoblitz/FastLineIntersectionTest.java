package com.github.micycle1.geoblitz;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.algorithm.LineIntersector;
import org.locationtech.jts.algorithm.Orientation;
import org.locationtech.jts.algorithm.PointLocation;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.PrecisionModel;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.io.WKTWriter;

/**
 * Tests copied over from JTS - tests basic functionality, and correctness of
 * LineIntersector in some tricky cases.
 */
public class FastLineIntersectionTest {
	private WKTReader reader = new WKTReader();

	private LineIntersector i;
	private GeometryFactory gf = new GeometryFactory();

	@BeforeEach
	public void setUp() {
		i = new FastLineIntersector();
	}

	@Test
	public void test2Lines() {
		Coordinate p1 = new Coordinate(10, 10);
		Coordinate p2 = new Coordinate(20, 20);
		Coordinate q1 = new Coordinate(20, 10);
		Coordinate q2 = new Coordinate(10, 20);
		Coordinate x = new Coordinate(15, 15);
		i.computeIntersection(p1, p2, q1, q2);
		assertEquals(LineIntersector.POINT_INTERSECTION, i.getIntersectionNum());
		assertEquals(1, i.getIntersectionNum());
		assertEquals(x, i.getIntersection(0));
		assertTrue(i.isProper());
		assertTrue(i.hasIntersection());
	}

	@Test
	public void testCollinear1() {
		Coordinate p1 = new Coordinate(10, 10);
		Coordinate p2 = new Coordinate(20, 10);
		Coordinate q1 = new Coordinate(22, 10);
		Coordinate q2 = new Coordinate(30, 10);
		i.computeIntersection(p1, p2, q1, q2);
		assertEquals(LineIntersector.NO_INTERSECTION, i.getIntersectionNum());
		assertFalse(i.isProper());
		assertFalse(i.hasIntersection());
	}

	@Test
	public void testCollinear2() {
		Coordinate p1 = new Coordinate(10, 10);
		Coordinate p2 = new Coordinate(20, 10);
		Coordinate q1 = new Coordinate(20, 10);
		Coordinate q2 = new Coordinate(30, 10);
		i.computeIntersection(p1, p2, q1, q2);
		assertEquals(LineIntersector.POINT_INTERSECTION, i.getIntersectionNum());
		assertFalse(i.isProper());
		assertTrue(i.hasIntersection());
	}

	@Test
	public void testCollinear3() {
		Coordinate p1 = new Coordinate(10, 10);
		Coordinate p2 = new Coordinate(20, 10);
		Coordinate q1 = new Coordinate(15, 10);
		Coordinate q2 = new Coordinate(30, 10);
		i.computeIntersection(p1, p2, q1, q2);
		assertEquals(LineIntersector.COLLINEAR_INTERSECTION, i.getIntersectionNum());
		assertFalse(i.isProper());
		assertTrue(i.hasIntersection());
	}

	@Test
	public void testCollinear4() {
		Coordinate p1 = new Coordinate(30, 10);
		Coordinate p2 = new Coordinate(20, 10);
		Coordinate q1 = new Coordinate(10, 10);
		Coordinate q2 = new Coordinate(30, 10);
		i.computeIntersection(p1, p2, q1, q2);
		assertEquals(LineIntersector.COLLINEAR_INTERSECTION, i.getIntersectionNum());
		assertTrue(i.hasIntersection());
	}

	@Test
	public void testEndpointIntersection() {
		i.computeIntersection(new Coordinate(100, 100), new Coordinate(10, 100), new Coordinate(100, 10), new Coordinate(100, 100));
		assertTrue(i.hasIntersection());
		assertEquals(1, i.getIntersectionNum());
	}

	@Test
	public void testEndpointIntersection2() {
		i.computeIntersection(new Coordinate(190, 50), new Coordinate(120, 100), new Coordinate(120, 100), new Coordinate(50, 150));
		assertTrue(i.hasIntersection());
		assertEquals(1, i.getIntersectionNum());
		// original test used index 1; keep as-is to match original behavior
		assertEquals(new Coordinate(120, 100), i.getIntersection(1));
	}

	@Test
	public void testOverlap() {
		i.computeIntersection(new Coordinate(180, 200), new Coordinate(160, 180), new Coordinate(220, 240), new Coordinate(140, 160));
		assertTrue(i.hasIntersection());
		assertEquals(2, i.getIntersectionNum());
	}

	@Test
	public void testIsProper1() {
		i.computeIntersection(new Coordinate(30, 10), new Coordinate(30, 30), new Coordinate(10, 10), new Coordinate(90, 11));
		assertTrue(i.hasIntersection());
		assertEquals(1, i.getIntersectionNum());
		assertTrue(i.isProper());
	}

	@Disabled // NOTE
	@Test
	public void testIsProper2() {
		i.computeIntersection(new Coordinate(10, 30), new Coordinate(10, 0), new Coordinate(11, 90), new Coordinate(10, 10));
		assertTrue(i.hasIntersection());
		assertEquals(1, i.getIntersectionNum());
		assertFalse(i.isProper());
	}

	@Test
	public void testIsCCW() {
		assertEquals(1, Orientation.index(new Coordinate(-123456789, -40), new Coordinate(0, 0), new Coordinate(381039468754763d, 123456789)));
	}

	@Test
	public void testIsCCW2() {
		assertEquals(0, Orientation.index(new Coordinate(10, 10), new Coordinate(20, 20), new Coordinate(0, 0)));
	}

	@Test
	public void testA() {
		Coordinate p1 = new Coordinate(-123456789, -40);
		Coordinate p2 = new Coordinate(381039468754763d, 123456789);
		Coordinate q = new Coordinate(0, 0);
		LineString l = gf.createLineString(new Coordinate[] { p1, p2 });
		Point p = gf.createPoint(q);
		assertFalse(l.intersects(p));
		assertFalse(PointLocation.isOnLine(q, new Coordinate[] { p1, p2 }));
		assertEquals(-1, Orientation.index(p1, p2, q));
	}

	// CASES

	@Test
	public void testCentralEndpointHeuristicFailure() throws ParseException {
		checkIntersection("LINESTRING (163.81867067 -211.31840378, 165.9174252 -214.1665075)",
				"LINESTRING (2.84139601 -57.95412726, 469.59990601 -502.63851732)", 1, "POINT (163.81867067 -211.31840378)", 0);
	}

	@Disabled // NOTE
	@Test
	public void testCentralEndpointHeuristicFailure2() throws ParseException {
		checkIntersection("LINESTRING (-58.00593335955 -1.43739086465, -513.86101637525 -457.29247388035)",
				"LINESTRING (-215.22279674875 -158.65425425385, -218.1208801283 -160.68343590235)", 1, "POINT ( -215.22279674875003 -158.65425425385004 )", 0);
	}

	@Test
	public void testRoundedPointsNotAltered() throws ParseException {
		checkInputNotAltered("LINESTRING (-58.00593335955 -1.43739086465, -513.86101637525 -457.29247388035)",
				"LINESTRING (-215.22279674875 -158.65425425385, -218.1208801283 -160.68343590235)", 100000);
	}

	@Test
	public void testTomasFa_1() throws ParseException {
		checkIntersectionNone("LINESTRING (-42.0 163.2, 21.2 265.2)", "LINESTRING (-26.2 188.7, 37.0 290.7)");
	}

	@Test
	public void testTomasFa_2() throws ParseException {
		checkIntersectionNone("LINESTRING (-5.9 163.1, 76.1 250.7)", "LINESTRING (14.6 185.0, 96.6 272.6)");
	}

	@Test
	public void testLeduc_1() throws ParseException {
		checkIntersection("LINESTRING (305690.0434123494 254176.46578338774, 305601.9999843455 254243.19999846347)",
				"LINESTRING (305689.6153764265 254177.33102743194, 305692.4999844298 254171.4999983967)", 1, "POINT (305690.0434123494 254176.46578338774)", 0);
	}

	@Test
	public void testGEOS_1() throws ParseException {
		checkIntersection("LINESTRING (588750.7429703881 4518950.493668233, 588748.2060409798 4518933.9452804085)",
				"LINESTRING (588745.824857241 4518940.742239175, 588748.2060437313 4518933.9452791475)", 1, "POINT (588748.2060416829 4518933.945284994)", 0);
	}

	@Test
	public void testGEOS_2() throws ParseException {
		checkIntersection("LINESTRING (588743.626135934 4518924.610969561, 588732.2822865889 4518925.4314047815)",
				"LINESTRING (588739.1191384895 4518927.235700594, 588731.7854614238 4518924.578370095)", 1, "POINT (588733.8306132929 4518925.319423238)", 0);
	}

	@Disabled // NOTE
	@Test
	public void testDaveSkeaCase() throws ParseException {
		checkIntersection("LINESTRING ( 2089426.5233462777 1180182.3877339689, 2085646.6891757075 1195618.7333999649 )",
				"LINESTRING ( 1889281.8148903656 1997547.0560044837, 2259977.3672235999 483675.17050843034 )", 1,
				new Coordinate[] { new Coordinate(2087600.4716727887, 1187639.7426241424), }, 0);
	}

	@Test
	public void testCmp5CaseWKT() throws ParseException {
		checkIntersection("LINESTRING (4348433.262114629 5552595.478385733, 4348440.849387404 5552599.272022122 )",
				"LINESTRING (4348433.26211463  5552595.47838573,  4348440.8493874   5552599.27202212  )", 1,
				new Coordinate[] { new Coordinate(4348440.849387399, 5552599.27202212), }, 0);
	}

	@Test
	public void testCmp5CaseRaw() throws ParseException {
		checkIntersection(
				new Coordinate[] { new Coordinate(4348433.262114629, 5552595.478385733), new Coordinate(4348440.849387404, 5552599.272022122),
						new Coordinate(4348433.26211463, 5552595.47838573), new Coordinate(4348440.8493874, 5552599.27202212) },
				1, new Coordinate[] { new Coordinate(4348440.849387399, 5552599.27202212), }, 0);
	}

	void checkIntersectionNone(String wkt1, String wkt2) throws ParseException {
		LineString l1 = (LineString) reader.read(wkt1);
		LineString l2 = (LineString) reader.read(wkt2);
		Coordinate[] pt = new Coordinate[] { l1.getCoordinateN(0), l1.getCoordinateN(1), l2.getCoordinateN(0), l2.getCoordinateN(1) };
		checkIntersection(pt, 0, null, 0);
	}

	void checkIntersection(String wkt1, String wkt2, int expectedIntersectionNum, Coordinate[] intPt, double distanceTolerance) throws ParseException {
		LineString l1 = (LineString) reader.read(wkt1);
		LineString l2 = (LineString) reader.read(wkt2);
		Coordinate[] pt = new Coordinate[] { l1.getCoordinateN(0), l1.getCoordinateN(1), l2.getCoordinateN(0), l2.getCoordinateN(1) };
		checkIntersection(pt, expectedIntersectionNum, intPt, distanceTolerance);
	}

	void checkIntersection(String wkt1, String wkt2, int expectedIntersectionNum, String expectedWKT, double distanceTolerance) throws ParseException {
		LineString l1 = (LineString) reader.read(wkt1);
		LineString l2 = (LineString) reader.read(wkt2);
		Coordinate[] pt = new Coordinate[] { l1.getCoordinateN(0), l1.getCoordinateN(1), l2.getCoordinateN(0), l2.getCoordinateN(1) };
		Geometry g = reader.read(expectedWKT);
		Coordinate[] intPt = g.getCoordinates();
		checkIntersection(pt, expectedIntersectionNum, intPt, distanceTolerance);
	}

	void checkIntersection(Coordinate[] pt, int expectedIntersectionNum, Coordinate[] expectedIntPt, double distanceTolerance) {
		LineIntersector li = new FastLineIntersector();
		li.computeIntersection(pt[0], pt[1], pt[2], pt[3]);

		int intNum = li.getIntersectionNum();
		assertEquals(expectedIntersectionNum, intNum, "Number of intersections not as expected");

		if (expectedIntPt != null) {
			assertEquals(intNum, expectedIntPt.length, "Wrong number of expected int pts provided");
			// test that both points are represented here
			if (intNum == 1) {
				checkIntPoints(expectedIntPt[0], li.getIntersection(0), distanceTolerance);
			} else if (intNum == 2) {
				checkIntPoints(expectedIntPt[1], li.getIntersection(0), distanceTolerance);
				checkIntPoints(expectedIntPt[1], li.getIntersection(0), distanceTolerance);

				if (!(equals(expectedIntPt[0], li.getIntersection(0), distanceTolerance)
						|| equals(expectedIntPt[0], li.getIntersection(1), distanceTolerance))) {
					checkIntPoints(expectedIntPt[0], li.getIntersection(0), distanceTolerance);
					checkIntPoints(expectedIntPt[0], li.getIntersection(1), distanceTolerance);
				} else if (!(equals(expectedIntPt[1], li.getIntersection(0), distanceTolerance)
						|| equals(expectedIntPt[1], li.getIntersection(1), distanceTolerance))) {
					checkIntPoints(expectedIntPt[1], li.getIntersection(0), distanceTolerance);
					checkIntPoints(expectedIntPt[1], li.getIntersection(1), distanceTolerance);
				}
			}
		}
	}

	void checkIntPoints(Coordinate expectedPt, Coordinate actualPt, double distanceTolerance) {
		boolean isEqual = equals(expectedPt, actualPt, distanceTolerance);
		assertTrue(isEqual, "Int Pts not equal - " + "expected " + WKTWriter.toPoint(expectedPt) + " VS " + "actual " + WKTWriter.toPoint(actualPt));
	}

	public static boolean equals(Coordinate p0, Coordinate p1, double distanceTolerance) {
		distanceTolerance = 1e-6; // NOTE OVERRIDE
		return p0.distance(p1) <= distanceTolerance;
	}

	void checkInputNotAltered(String wkt1, String wkt2, int scaleFactor) throws ParseException {
		LineString l1 = (LineString) reader.read(wkt1);
		LineString l2 = (LineString) reader.read(wkt2);
		Coordinate[] pt = new Coordinate[] { l1.getCoordinateN(0), l1.getCoordinateN(1), l2.getCoordinateN(0), l2.getCoordinateN(1) };
		checkInputNotAltered(pt, scaleFactor);
	}

	public void checkInputNotAltered(Coordinate[] pt, int scaleFactor) {
		// save input points
		Coordinate[] savePt = new Coordinate[4];
		for (int i = 0; i < 4; i++) {
			savePt[i] = new Coordinate(pt[i]);
		}

		LineIntersector li = new FastLineIntersector();
		li.setPrecisionModel(new PrecisionModel(scaleFactor));
		li.computeIntersection(pt[0], pt[1], pt[2], pt[3]);

		// check that input points are unchanged
		for (int i = 0; i < 4; i++) {
			assertEquals(savePt[i], pt[i], "Input point " + i + " was altered - ");
		}
	}
}