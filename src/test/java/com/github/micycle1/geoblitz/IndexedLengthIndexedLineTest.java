package com.github.micycle1.geoblitz;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.linearref.LengthIndexedLine;

/**
 * Tests equivalence between LengthIndexedLine and IndexedLengthIndexedLine.
 */
public class IndexedLengthIndexedLineTest {

	private static final double EPS = 1e-6;
	private final GeometryFactory gf = new GeometryFactory();
	private final Random rand = new Random(12345L);

	@Test
	public void testExtractPointDeterministic() {
		List<Geometry> geoms = createDeterministicGeometries();

		for (Geometry g : geoms) {
			LengthIndexedLine base = new LengthIndexedLine(g);
			IndexedLengthIndexedLine idx = new IndexedLengthIndexedLine(g);

			double total = g.getLength();

			// test a set of boundary and interior positions
			double[] tests = new double[] { 0.0, total * 0.0001, total * 0.1, total * 0.25, total * 0.5, total * 0.75, total - 1e-12, total, total + 1.0, // out
																																							// of
																																							// range
					-1.0, // negative from end
					-total * 0.5, -total - 2.0 // out of range negative
			};

			for (double t : tests) {
				Coordinate cBase = base.extractPoint(t);
				Coordinate cIdx = idx.extractPoint(t);
				assertCoordinatesEqual(g, cBase, cIdx, EPS, "extractPoint mismatch for length " + t + " on geom " + g);
			}
		}
	}

	@Test
	public void testExtractPointWithOffsetDeterministic() {
		List<Geometry> geoms = createDeterministicGeometries();
		double[] offsets = new double[] { -2.0, -0.5, 0.0, 0.5, 2.0 };

		for (Geometry g : geoms) {
			LengthIndexedLine base = new LengthIndexedLine(g);
			IndexedLengthIndexedLine idx = new IndexedLengthIndexedLine(g);

			double total = g.getLength();

			double[] tests = new double[] { 0.0, total * 0.125, total * 0.5, total - 1e-12, total, -1.0 };

			for (double t : tests) {
				for (double off : offsets) {
					Coordinate cBase = base.extractPoint(t, off);
					Coordinate cIdx = idx.extractPoint(t, off);
					assertCoordinatesEqual(g, cBase, cIdx, EPS, "extractPoint(index, offset) mismatch for index " + t + " offset " + off);
				}
			}
		}
	}

	@Test
	public void testExtractLineDeterministic() {
		List<Geometry> geoms = createDeterministicGeometries();

		for (Geometry g : geoms) {
			LengthIndexedLine base = new LengthIndexedLine(g);
			IndexedLengthIndexedLine idx = new IndexedLengthIndexedLine(g);

			double total = g.getLength();
			double[] testPoints = new double[] { 0.0, total * 0.1, total * 0.25, total * 0.5, total * 0.75, total - 1e-12, total, -0.5, -total * 0.2 };

			for (double a : testPoints) {
				for (double b : testPoints) {
					Geometry gBase = base.extractLine(a, b);
					Geometry gIdx = idx.extractLine(a, b);
					assertGeometriesEquivalent(gBase, gIdx, EPS, "extractLine mismatch for start " + a + " end " + b + " on geom " + g);
				}
			}
		}
	}

	@Test
	public void testRandomizedComparisons() {
		// Random tests over many random polylines and random indices
		for (int iter = 0; iter < 30; iter++) {
			Geometry g = createRandomPolyline(50, 10.0);
			LengthIndexedLine base = new LengthIndexedLine(g);
			IndexedLengthIndexedLine idx = new IndexedLengthIndexedLine(g);

			double total = g.getLength();

			// random indices (some negative, some out-of-range)
			for (int k = 0; k < 100; k++) {
				double a = (rand.nextDouble() * (total * 2.0)) - (total * 0.25); // roughly [-0.25T,1.75T)
				double b = (rand.nextDouble() * (total * 2.0)) - (total * 0.25);

				Coordinate ca = base.extractPoint(a);
				Coordinate cb = idx.extractPoint(a);
				assertCoordinatesEqual(g, ca, cb, EPS, "random extractPoint mismatch at index " + a);

				// offset points
				double off = (rand.nextDouble() - 0.5) * 4.0;
				Coordinate caOff = base.extractPoint(a, off);
				Coordinate cbOff = idx.extractPoint(a, off);
				assertCoordinatesEqual(g, caOff, cbOff, EPS, "random extractPoint with offset mismatch at index " + a + " offset " + off);

				Geometry gl = base.extractLine(a, b);
				Geometry gi = idx.extractLine(a, b);
				assertGeometriesEquivalent(gl, gi, EPS, "random extractLine mismatch start " + a + " end " + b);
			}
		}
	}

	// ---------- Helpers ----------

	private List<Geometry> createDeterministicGeometries() {
		List<Geometry> geoms = new ArrayList<>();

		// simple polyline
		geoms.add(gf.createLineString(new Coordinate[] { new Coordinate(0, 0), new Coordinate(10, 0), new Coordinate(10, 10) }));

		// multi-line (two components)
		LineString c0 = gf.createLineString(new Coordinate[] { new Coordinate(0, 0), new Coordinate(5, 0), new Coordinate(5, 5) });
		LineString c1 = gf.createLineString(new Coordinate[] { new Coordinate(5, 5), new Coordinate(10, 5) });
		MultiLineString mls = gf.createMultiLineString(new LineString[] { c0, c1 });
		geoms.add(mls);

		// line with zero-length component
		LineString zeroComp = gf.createLineString(new Coordinate[] { new Coordinate(100, 100), new Coordinate(100, 100) });
		LineString normal = gf.createLineString(new Coordinate[] { new Coordinate(100, 100), new Coordinate(105, 105) });
		geoms.add(gf.createMultiLineString(new LineString[] { normal }));

		// longer fixed polyline
		Coordinate[] coords = new Coordinate[7];
		coords[0] = new Coordinate(0, 0);
		coords[1] = new Coordinate(1, 2);
		coords[2] = new Coordinate(3, 3);
		coords[3] = new Coordinate(6, 2);
		coords[4] = new Coordinate(6, 6);
		coords[5] = new Coordinate(10, 6);
		coords[6] = new Coordinate(10, 10);
		geoms.add(gf.createLineString(coords));

		return geoms;
	}

	private Geometry createRandomPolyline(int nPoints, double step) {
		Coordinate[] coords = new Coordinate[nPoints];
		double x = 0.0, y = 0.0;
		for (int i = 0; i < nPoints; i++) {
			x += (rand.nextDouble() - 0.5) * step;
			y += (rand.nextDouble() - 0.5) * step;
			coords[i] = new Coordinate(x, y);
		}
		return gf.createLineString(coords);
	}

	private void assertCoordinatesEqual(Geometry context, Coordinate expected, Coordinate actual, double eps, String message) {
		assertNotNull(expected, "expected coordinate was null; context: " + context);
		assertNotNull(actual, "actual coordinate was null; context: " + context);
		double d = expected.distance(actual);
		assertTrue(d <= eps, message + " (distance = " + d + ", eps=" + eps + ") expected=" + coordToStr(expected) + " actual=" + coordToStr(actual));
	}

	private void assertGeometriesEquivalent(Geometry g1, Geometry g2, double eps, String message) {
		assertNotNull(g1, "first geometry is null");
		assertNotNull(g2, "second geometry is null");

		// Quick checks: type and length within tolerance
		assertEquals(g1.getGeometryType(), g2.getGeometryType(), message + " - geometry type differs");
		double len1 = g1.getLength();
		double len2 = g2.getLength();
		assertEquals(len1, len2, eps, message + " - geometry length differs");

		Coordinate[] c1 = g1.getCoordinates();
		Coordinate[] c2 = g2.getCoordinates();

		// Some implementations may include extra duplicate points at component
		// boundaries;
		// to be robust, compare sequences allowing possible equal leading/trailing
		// duplicates.
		// We'll require the coordinate sequences' lengths to be equal and pairwise
		// close.
		assertEquals(c1.length, c2.length, message + " - coordinate array length differs (" + c1.length + " vs " + c2.length + ")");

		for (int i = 0; i < c1.length; i++) {
			double d = c1[i].distance(c2[i]);
			assertTrue(d <= eps,
					message + " - coordinate[" + i + "] differs by " + d + " (expected=" + coordToStr(c1[i]) + ", actual=" + coordToStr(c2[i]) + ")");
		}
	}

	private String coordToStr(Coordinate c) {
		if (c == null) {
			return "null";
		}
		return String.format("(%.9f,%.9f,%.9f)", c.x, c.y, c.getZ());
	}
}