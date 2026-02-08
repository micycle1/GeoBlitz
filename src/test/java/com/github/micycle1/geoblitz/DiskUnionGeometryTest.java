package com.github.micycle1.geoblitz;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.algorithm.distance.DiscreteHausdorffDistance;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.operation.union.UnaryUnionOp;

public class DiskUnionTest {

	private static final double EPS = 1e-9;
	private static final double MAX_SEG_LEN = 0.05;
	private static final GeometryFactory GF = new GeometryFactory();

	@Test
	public void testEmpty() {
		DiskUnion.ArcBoundary boundary = DiskUnion.computeBoundaryArcs(List.of(), EPS);
		Geometry actual = DiskUnion.toJtsGeometry(boundary, GF, MAX_SEG_LEN);
		assertTrue(actual.isEmpty(), "Expected empty geometry for empty input");
	}

	@Test
	public void testDisjointCircles() {
		List<Coordinate> circles = List.of(circle(0, 0, 1), circle(5, 0, 1));
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testOverlappingPair() {
		List<Coordinate> circles = List.of(circle(0, 0, 2), circle(2, 0, 2));
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testTangentPair() {
		List<Coordinate> circles = List.of(circle(0, 0, 2), circle(4, 0, 2));
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testContainedCircle() {
		List<Coordinate> circles = List.of(circle(0, 0, 5), circle(1, 1, 1));
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testRandomCluster() {
		Random rnd = new Random(1337);
		for (int t = 0; t < 50; t++) {
			List<Coordinate> circles = new ArrayList<>();
			for (int i = 0; i < 15 + rnd.nextInt(20); i++) {
				double x = rnd.nextDouble() * 10.0;
				double y = rnd.nextDouble() * 10.0;
				double r = 0.5 + rnd.nextDouble() * 1.5;
				circles.add(circle(x, y, r));
			}
			assertMatchesJtsUnion(circles);
		}
	}

	private static void assertMatchesJtsUnion(List<Coordinate> circles) {
		Geometry expected = unionWithJts(circles, MAX_SEG_LEN);
		Geometry actual = unionWithDiskUnion(circles, EPS, MAX_SEG_LEN);

		assertEquals(expected.getNumGeometries(), actual.getNumGeometries());

		double symDiffArea = expected.symDifference(actual).getArea();
		double areaTol = Math.max(1e-4, expected.getArea() * 5e-3);

		double hausdorff = DiscreteHausdorffDistance.distance(expected, actual);
		double distTol = Math.max(1e-3, MAX_SEG_LEN * 2.0);

		boolean ok = symDiffArea <= areaTol && hausdorff <= distTol;
		assertTrue(ok, () -> "DiskUnion mismatch vs JTS union\n" + "SymDiffArea=" + symDiffArea + " (tol " + areaTol + ")\n" + "Hausdorff=" + hausdorff
				+ " (tol " + distTol + ")\n");
	}

	private static Geometry unionWithDiskUnion(List<Coordinate> circles, double eps, double maxSegLen) {
		DiskUnion.ArcBoundary boundary = DiskUnion.computeBoundaryArcs(circles, eps);
		return DiskUnion.toJtsGeometry(boundary, GF, maxSegLen).buffer(0);
	}

	private static Geometry unionWithJts(List<Coordinate> circles, double maxSegLen) {
		List<Geometry> geoms = new ArrayList<>();
		for (Coordinate c : circles) {
			geoms.add(circleToPolygon(c, maxSegLen));
		}
		return UnaryUnionOp.union(geoms);
	}

	private static Geometry circleToPolygon(Coordinate c, double maxSegLen) {
		double circumference = Math.PI * 2.0 * c.getZ();
		int segments = Math.max(16, (int) Math.ceil(circumference / maxSegLen));
		int quadSegs = Math.max(4, (int) Math.ceil(segments / 4.0));
		return GF.createPoint(new Coordinate(c.x, c.y)).buffer(c.getZ(), quadSegs);
	}

	private static Coordinate circle(double x, double y, double r) {
		return new Coordinate(x, y, r);
	}
}
