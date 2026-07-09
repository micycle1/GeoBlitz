package com.github.micycle1.geoblitz;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
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

public class CircleUnionTest {

	private static final double MAX_SEG_LEN = 0.05;
	private static final GeometryFactory GF = new GeometryFactory();

	@Test
	public void testEmpty() {
		Geometry actual = CircleUnion.union(List.of(), MAX_SEG_LEN);
		assertTrue(actual.isEmpty(), "Expected empty geometry for empty input");
	}

	@Test
	public void testSingleCircle() {
		List<Coordinate> circles = List.of(circle(1, 2, 3));
		assertMatchesJtsUnion(circles);
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
	public void testChainOfThree() {
		List<Coordinate> circles = List.of(circle(0, 0, 2), circle(3, 0, 2), circle(6, 0, 2));
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
	public void testRandomVariedRadii() {
		// mixes big engulfing circles with small ones
		Random rnd = new Random(7);
		for (int t = 0; t < 25; t++) {
			List<Coordinate> circles = new ArrayList<>();
			for (int i = 0; i < 30; i++) {
				double x = rnd.nextDouble() * 20.0;
				double y = rnd.nextDouble() * 20.0;
				double r = 0.1 + rnd.nextDouble() * 5.0;
				circles.add(circle(x, y, r));
			}
			assertMatchesJtsUnion(circles);
		}
	}

	@Test
	public void testDuplicateWithOverlap() {
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
		List<Coordinate> circles = List.of(circle(0, 0, 5), circle(0, 0, 2));
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testNearTangentExternal() {
		double r = 2.0;
		double d = 2.0 * r + 1e-7; // slightly disjoint
		List<Coordinate> circles = List.of(circle(0, 0, r), circle(d, 0, r));
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testNearTangentOverlapping() {
		double r = 2.0;
		double d = 2.0 * r - 1e-7; // slightly overlapping
		List<Coordinate> circles = List.of(circle(0, 0, r), circle(d, 0, r));
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testNearlyCoincidentCentersSameRadius() {
		double r = 3.0;
		List<Coordinate> circles = List.of(circle(0, 0, r), circle(5e-10, -5e-10, r));
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testRingOfDisksCreatesHole() {
		// Disks arranged in a ring: union should be a polygon with a hole.
		int n = 8;
		double R = 3.0;
		double r = 1.4;
		List<Coordinate> circles = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			double a = (2.0 * Math.PI * i) / n;
			circles.add(circle(R * Math.cos(a), R * Math.sin(a), r));
		}
		Geometry actual = CircleUnion.union(circles, MAX_SEG_LEN);
		assertEquals(1, actual.getNumGeometries(), "Expected a single polygon");
		assertEquals(1, ((org.locationtech.jts.geom.Polygon) actual.getGeometryN(0)).getNumInteriorRing(), "Expected one hole");
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testIslandInHoleIsSeparatePolygon() {
		// Ring of disks with an isolated disk inside the hole: two components,
		// nested by connectivity (island must be its own polygon, not a hole).
		int n = 10;
		double R = 5.0;
		double r = 1.8;
		List<Coordinate> circles = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			double a = (2.0 * Math.PI * i) / n;
			circles.add(circle(R * Math.cos(a), R * Math.sin(a), r));
		}
		circles.add(circle(0, 0, 1.0)); // island inside the hole
		Geometry actual = CircleUnion.union(circles, MAX_SEG_LEN);
		assertEquals(2, actual.getNumGeometries(), "Expected two polygons (ring + island)");
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testLargeCoordinatesStress() {
		double base = 1e9;
		List<Coordinate> circles = List.of(circle(base, base, 2), //
				circle(base + 3.1, base, 2), // overlap
				circle(base + 20.0, base + 1, 1) // disjoint
		);
		assertMatchesJtsUnion(circles);
	}

	@Test
	public void testZeroRadiusCirclesIgnored() {
		List<Coordinate> circles = List.of(circle(0, 0, 2), circle(1, 0, 0), circle(9, 9, 0));
		Geometry actual = CircleUnion.union(circles, MAX_SEG_LEN);
		Geometry expected = CircleUnion.union(List.of(circle(0, 0, 2)), MAX_SEG_LEN);
		assertGeomsClose(expected, actual);
	}

	@Test
	public void testPermutationInvariance() {
		List<Coordinate> circles = List.of(circle(0, 0, 2), circle(2, 0, 2), circle(1, 2, 1.5), circle(5, 0, 1));

		Geometry base = CircleUnion.union(circles, MAX_SEG_LEN);

		Random rnd = new Random(123);
		for (int t = 0; t < 25; t++) {
			List<Coordinate> shuffled = new ArrayList<>(circles);
			Collections.shuffle(shuffled, rnd);
			Geometry g = CircleUnion.union(shuffled, MAX_SEG_LEN);
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

		Geometry g1 = CircleUnion.union(circles, MAX_SEG_LEN);
		Geometry g2 = CircleUnion.union(duplicated, MAX_SEG_LEN);
		assertGeomsClose(g1, g2);
	}

	@Test
	public void testBuilderReuseAndInvalidation() {
		CircleUnion cu = new CircleUnion(3);
		cu.add(0, 0, 2);
		cu.add(2, 0, 2);
		Geometry g1 = cu.geometry(GF, MAX_SEG_LEN);
		// second read reuses the cached topology
		Geometry g2 = cu.geometry(GF, MAX_SEG_LEN);
		assertGeomsClose(g1, g2);

		// adding invalidates the cache
		cu.add(10, 0, 1);
		Geometry g3 = cu.geometry(GF, MAX_SEG_LEN);
		assertTrue(g3.getArea() > g1.getArea(), "Union should grow after adding a disjoint circle");
	}

	@Test
	public void testBuilderGrowsBeyondInitialCapacity() {
		// initial capacity of 1 forces repeated growth; also exercises the
		// no-arg constructor default path via a second instance.
		CircleUnion cu = new CircleUnion(1);
		CircleUnion cuNoArg = new CircleUnion();
		List<Coordinate> circles = new ArrayList<>();
		Random rnd = new Random(5);
		for (int i = 0; i < 100; i++) {
			double x = rnd.nextDouble() * 10.0;
			double y = rnd.nextDouble() * 10.0;
			double r = 0.5 + rnd.nextDouble();
			circles.add(circle(x, y, r));
			cu.add(x, y, r);
			cuNoArg.add(x, y, r);
		}
		assertEquals(100, cu.size());
		Geometry expected = unionWithJts(circles, MAX_SEG_LEN);
		assertGeomsClose(expected, cu.geometry(GF, MAX_SEG_LEN));
		assertGeomsClose(expected, cuNoArg.geometry(GF, MAX_SEG_LEN));
	}

	@Test
	public void testInvalidInputsRejected() {
		CircleUnion cu = new CircleUnion(4);
		assertThrows(IllegalArgumentException.class, () -> cu.add(0, 0, Double.NaN));
		assertThrows(IllegalArgumentException.class, () -> cu.add(0, 0, -1));
		assertThrows(IllegalArgumentException.class, () -> cu.add(Double.NaN, 0, 1));
		assertThrows(IllegalArgumentException.class, () -> new CircleUnion(-1));
		assertThrows(IllegalArgumentException.class, () -> CircleUnion.union(List.of(new Coordinate(0, 0)), MAX_SEG_LEN)); // no z
	}

	private static void assertMatchesJtsUnion(List<Coordinate> circles) {
		Geometry expected = unionWithJts(circles, MAX_SEG_LEN);
		Geometry actual = CircleUnion.union(circles, MAX_SEG_LEN);
		IsValidOp validator = new IsValidOp(actual);
		String m = validator.isValid() ? "" : validator.getValidationError().getMessage();
		assertTrue(validator.isValid(), "Output is invalid geometry: " + m);

		double symDiffArea = expected.symDifference(actual).getArea();
		double areaTol = Math.max(1e-4, expected.getArea() * 5e-3);

		double hausdorff = DiscreteHausdorffDistance.distance(expected, actual);
		double distTol = Math.max(1e-3, MAX_SEG_LEN * 2.0);

		boolean ok = symDiffArea <= areaTol && hausdorff <= distTol;
		assertTrue(ok, () -> "CircleUnion mismatch vs JTS union\n" + "SymDiffArea=" + symDiffArea + " (tol " + areaTol + ")\n" + "Hausdorff=" + hausdorff
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
