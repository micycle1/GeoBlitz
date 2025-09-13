package com.github.micycle1.geoblitz;

import static org.junit.jupiter.api.Assertions.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.LineSegment;

/**
 * Unit tests for HPRtreeX.nearestNeighbor comparing results against a
 * brute-force search.
 * <p>
 * Tests two item types: - Coordinate: items are points (envelope is the point
 * itself) - LineSegment: items are segments (envelope covers the segment
 * endpoints)
 * <p>
 * For each query the test computes the brute-force nearest item (minimum
 * distance) and asserts that the tree's nearestNeighbor returns an item whose
 * distance equals the brute-force minimum within a small epsilon. If multiple
 * items tie for the minimum, the test allows any of those tied items as the
 * valid nearest neighbor.
 */
public class HPRtreeXTest {

	private static final double EPS = 1e-9;

	@Test
	public void testNearestNeighborWithCoordinates() {
		final int N_ITEMS = 5000;
		final int N_QUERIES = 20000;
		Random rnd = new Random(42);

		HPRtreeX<Coordinate> tree = new HPRtreeX<>(16);
		List<Coordinate> items = new ArrayList<>(N_ITEMS);

		// insert random points
		for (int i = 0; i < N_ITEMS; i++) {
			double x = rnd.nextDouble() * 1000.0;
			double y = rnd.nextDouble() * 1000.0;
			Coordinate c = new Coordinate(x, y);
			items.add(c);
			tree.insert(new Envelope(x, x, y, y), c);
		}

		// building explicitly (nearestNeighbor will also trigger build if not called)
		tree.build();

		// distance function: Euclidean distance between query and point
		HPRtreeX.DistanceToItem<Coordinate> distFn = (q, item) -> q.distance(item);

		for (int qi = 0; qi < N_QUERIES; qi++) {
			Coordinate q = new Coordinate(rnd.nextDouble() * 1000.0, rnd.nextDouble() * 1000.0);

			// brute-force search
			double minDist = Double.POSITIVE_INFINITY;
			List<Coordinate> minima = new ArrayList<>();
			for (Coordinate it : items) {
				double d = q.distance(it);
				if (d + EPS < minDist) {
					minDist = d;
					minima.clear();
					minima.add(it);
				} else if (Math.abs(d - minDist) <= EPS) {
					minima.add(it);
				}
			}

			Coordinate found = tree.nearestNeighbor(q, distFn);

			// Non-empty tree should always return a nearest item
			assertNotNull(found, "Tree returned null for non-empty dataset (query " + qi + ")");

			// found item's distance must equal brute-force minimum (within EPS)
			double foundDist = distFn.distance(q, found);
			assertEquals(minDist, foundDist, EPS, "Distance mismatch for query " + qi);

			// found should be one of the brute-force minima (handles ties)
			assertTrue(minima.contains(found), "Nearest returned by tree is not one of brute-force minima for query " + qi);
		}

		// Also test empty tree returns null
		HPRtreeX<Coordinate> emptyTree = new HPRtreeX<>();
		assertNull(emptyTree.nearestNeighbor(new Coordinate(0, 0), distFn), "Empty tree must return null");
	}

	@Test
	public void testNearestNeighborWithLineSegments() {
		final int N_ITEMS = 4000;
		final int N_QUERIES = 1500;
		Random rnd = new Random(123);

		HPRtreeX<LineSegment> tree = new HPRtreeX<>(16);
		List<LineSegment> items = new ArrayList<>(N_ITEMS);

		// insert random line segments
		for (int i = 0; i < N_ITEMS; i++) {
			double x1 = rnd.nextDouble() * 1000.0;
			double y1 = rnd.nextDouble() * 1000.0;
			double x2 = rnd.nextDouble() * 1000.0;
			double y2 = rnd.nextDouble() * 1000.0;
			Coordinate p0 = new Coordinate(x1, y1);
			Coordinate p1 = new Coordinate(x2, y2);
			LineSegment seg = new LineSegment(p0, p1);
			items.add(seg);
			tree.insert(new Envelope(p0, p1), seg);
		}

		tree.build();

		// distance function: distance from query point to segment
		HPRtreeX.DistanceToItem<LineSegment> distFn = (q, seg) -> seg.distance(q);

		for (int qi = 0; qi < N_QUERIES; qi++) {
			Coordinate q = new Coordinate(rnd.nextDouble() * 1000.0, rnd.nextDouble() * 1000.0);

			// brute-force search
			double minDist = Double.POSITIVE_INFINITY;
			List<LineSegment> minima = new ArrayList<>();
			for (LineSegment seg : items) {
				double d = seg.distance(q);
				if (d + EPS < minDist) {
					minDist = d;
					minima.clear();
					minima.add(seg);
				} else if (Math.abs(d - minDist) <= EPS) {
					minima.add(seg);
				}
			}

			LineSegment found = tree.nearestNeighbor(q, distFn);

			assertNotNull(found, "Tree returned null for non-empty dataset (query " + qi + ")");

			double foundDist = distFn.distance(q, found);
			assertEquals(minDist, foundDist, EPS, "Distance mismatch for query " + qi);
			assertEquals(minima.get(0), found);

			assertTrue(minima.contains(found), "Nearest returned by tree is not one of brute-force minima for query " + qi);
		}

		// empty tree case
		HPRtreeX<LineSegment> emptyTree = new HPRtreeX<>();
		assertNull(emptyTree.nearestNeighbor(new Coordinate(0, 0), distFn), "Empty tree must return null");
	}
}