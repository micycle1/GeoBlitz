package com.github.micycle1.geoblitz;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.operation.union.CascadedPolygonUnion;

class HilbertParallelPolygonUnionTest {

	private static final GeometryFactory GF = new GeometryFactory();

	private static List<Geometry> randomBufferedPoints(int n, double radius, long seed) {
		Random rnd = new Random(seed);
		List<Geometry> polys = new ArrayList<>(n);
		for (int i = 0; i < n; i++) {
			var pt = GF.createPoint(new Coordinate(rnd.nextDouble() * 1000, rnd.nextDouble() * 1000));
			polys.add(pt.buffer(radius, 16));
		}
		return polys;
	}

	@Test
	void emptyInput() {
		Geometry result = HilbertParallelPolygonUnion.union(new ArrayList<Geometry>());
		assertTrue(result.isEmpty());
	}

	@Test
	void singlePolygon() {
		List<Geometry> input = randomBufferedPoints(1, 10, 1);
		Geometry result = HilbertParallelPolygonUnion.union(input);
		assertTrue(result.equalsTopo(input.get(0)));
	}

	@Test
	void matchesCascadedUnion_overlapping() {
		assertMatchesCascaded(randomBufferedPoints(500, 30, 42)); // dense overlap
	}

	@Test
	void matchesCascadedUnion_sparse() {
		assertMatchesCascaded(randomBufferedPoints(500, 2, 42)); // mostly disjoint
	}

	@Test
	void matchesCascadedUnion_moderate() {
		assertMatchesCascaded(randomBufferedPoints(1000, 10, 7));
	}

	private static void assertMatchesCascaded(List<Geometry> input) {
		Geometry expected = CascadedPolygonUnion.union(new ArrayList<>(input));
		Geometry actual = HilbertParallelPolygonUnion.union(new ArrayList<>(input));

		assertTrue(actual instanceof Polygonal, "result must be polygonal, was " + actual.getGeometryType());
		assertTrue(actual.isValid(), "result must be valid");
		assertEquals(expected.getNumGeometries(), actual.getNumGeometries(), "component count");
		assertEquals(expected.getArea(), actual.getArea(), expected.getArea() * 1e-9, "area");
		// symmetric difference should be (near) nothing
		double symDiff = expected.symDifference(actual).getArea();
		assertTrue(symDiff < expected.getArea() * 1e-9, "symDifference area was " + symDiff);
	}
}
