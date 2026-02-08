package com.github.micycle1.geoblitz;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.algorithm.distance.DiscreteHausdorffDistance;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.operation.valid.IsValidOp;

public class DiskUnionGeometryTest {

	private static final double EPS = 1e-10;
	private static final double MAX_SEG_LEN = 0.05;
	private static final GeometryFactory GF = new GeometryFactory();

	@Test
	public void testEmpty() {
		DiskUnion.ArcBoundary boundary = DiskUnion.computeBoundaryArcs(List.of(), EPS);
		Geometry actual = DiskUnion.toGeometry(boundary, GF, MAX_SEG_LEN);
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
	public void testDuplicatePair() {
		List<Coordinate> circles = List.of(circle(0, 0, 2), circle(0, 0, 2));
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

	@Test
	public void testDuplicateWithOverlap() {
		// Same as overlapping pair, but with a duplicate of one input.
		List<Coordinate> circles = List.of(circle(0, 0, 2), circle(0, 0, 2), // duplicate
				circle(2, 0, 2));
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testManyDuplicatesDoNotEraseUnion() {
		List<Coordinate> circles = new ArrayList<>();
		for (int i = 0; i < 20; i++) {
			circles.add(circle(0, 0, 2)); // all duplicates
		}
		circles.add(circle(6, 0, 1)); // plus one disjoint component
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testConcentricDifferentRadii() {
		// Same center; smaller is fully contained.
		List<Coordinate> circles = List.of(circle(0, 0, 5), circle(0, 0, 2));
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testNearTangentExternal() {
		// Almost tangent but disjoint by a tiny amount.
		// This stresses eps handling in circle-circle relation + angle splitting.
		double r = 2.0;
		double d = 2.0 * r + 1e-7; // slightly disjoint
		List<Coordinate> circles = List.of(circle(0, 0, r), circle(d, 0, r));
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testNearTangentOverlapping() {
		// Almost tangent but overlapping by a tiny amount.
		double r = 2.0;
		double d = 2.0 * r - 1e-7; // slightly overlapping
		List<Coordinate> circles = List.of(circle(0, 0, r), circle(d, 0, r));
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testNearlyCoincidentCentersSameRadius() {
		// Centers extremely close (<< EPS), same radius.
		// Stresses the d<=eps branch and snapping.
		double r = 3.0;
		List<Coordinate> circles = List.of(circle(0, 0, r), circle(5e-10, -5e-10, r));
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testRingOfDisksCreatesHole() {
		// Disks arranged in a ring: union should be a polygon with a hole (donut-like).
		// This stresses arc stitching + hole handling in toJtsGeometry().
		int n = 8;
		double R = 3.0; // ring radius (centers lie on this circle)
		double r = 1.4; // disk radius (overlaps neighbors, but does not reach origin)
		List<Coordinate> circles = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			double a = (2.0 * Math.PI * i) / n;
			circles.add(circle(R * Math.cos(a), R * Math.sin(a), r));
		}
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testLargeCoordinatesStress() {
		// Stresses snapping quantization/keying and numerical stability.
		double base = 1e9;
		List<Coordinate> circles = List.of(circle(base, base, 2), //
				circle(base + 3.1, base, 2), // overlap
				circle(base + 20.0, base + 1, 1) // disjoint
		);
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testPermutationInvariance() {
		List<Coordinate> circles = List.of(circle(0, 0, 2), circle(2, 0, 2), circle(1, 2, 1.5), circle(5, 0, 1));

		Geometry base = unionWithDiskUnion(circles, EPS, MAX_SEG_LEN);

		Random rnd = new Random(123);
		for (int t = 0; t < 25; t++) {
			List<Coordinate> shuffled = new ArrayList<>(circles);
			Collections.shuffle(shuffled, rnd);
			Geometry g = unionWithDiskUnion(shuffled, EPS, MAX_SEG_LEN);
			assertGeomsClose(base, g);
		}
	}

	@Test
	public void testIdempotenceUnderDuplication() {
		List<Coordinate> circles = List.of(circle(0, 0, 2), circle(2, 0, 2), circle(5, 0, 1));

		List<Coordinate> duplicated = new ArrayList<>();
		duplicated.addAll(circles);
		duplicated.addAll(circles);
		duplicated.addAll(circles);

		Geometry g1 = unionWithDiskUnion(circles, EPS, MAX_SEG_LEN);
		Geometry g2 = unionWithDiskUnion(duplicated, EPS, MAX_SEG_LEN);
		assertGeomsClose(g1, g2);
	}

	private static void assertMatchesJtsUnion(List<Coordinate> circles) {
		Geometry expected = unionWithJts(circles, MAX_SEG_LEN);
		Geometry actual = unionWithDiskUnion(circles, EPS, MAX_SEG_LEN);
		IsValidOp validator = new IsValidOp(actual);
		String m = validator.isValid() ? "" : validator.getValidationError().getMessage();
		assertTrue(validator.isValid(), "Output is invalid geometry: " + m);

		double symDiffArea = expected.symDifference(actual).getArea();
		double areaTol = Math.max(1e-4, expected.getArea() * 5e-3);

		double hausdorff = DiscreteHausdorffDistance.distance(expected, actual);
		double distTol = Math.max(1e-3, MAX_SEG_LEN * 2.0);

		boolean ok = symDiffArea <= areaTol && hausdorff <= distTol;
		assertTrue(ok, () -> "DiskUnion mismatch vs JTS union\n" + "SymDiffArea=" + symDiffArea + " (tol " + areaTol + ")\n" + "Hausdorff=" + hausdorff
				+ " (tol " + distTol + ")\n");
	}

	private static void assertGeomsClose(Geometry expected, Geometry actual) {
		double symDiffArea = expected.symDifference(actual).getArea();
		double areaTol = Math.max(1e-4, expected.getArea() * 5e-3);

		double hausdorff = DiscreteHausdorffDistance.distance(expected, actual);
		double distTol = Math.max(1e-3, MAX_SEG_LEN * 2.0);

		boolean ok = symDiffArea <= areaTol && hausdorff <= distTol;
		assertTrue(ok, () -> "Geometry mismatch\n" + "SymDiffArea=" + symDiffArea + " (tol " + areaTol + ")\n" + "Hausdorff=" + hausdorff + " (tol " + distTol
				+ ")\n");
	}

	private static Geometry unionWithDiskUnion(List<Coordinate> circles, double eps, double maxSegLen) {
		DiskUnion.ArcBoundary boundary = DiskUnion.computeBoundaryArcs(circles, eps);
		return DiskUnion.toGeometry(boundary, GF, maxSegLen);
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
