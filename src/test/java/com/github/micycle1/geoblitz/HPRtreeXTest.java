package com.github.micycle1.geoblitz;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

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

		tree.build();

		// distance functions
		HPRtreeX.DistanceToItem<Coordinate> distFn = (q, item) -> q.distance(item);
		HPRtreeX.DistanceToItem<Coordinate> distSqFn = (q, item) -> q.distanceSq(item);

		for (int qi = 0; qi < N_QUERIES; qi++) {
			Coordinate q = new Coordinate(rnd.nextDouble() * 1000.0, rnd.nextDouble() * 1000.0);

			// brute-force search (track minima to handle ties)
			double minDist = Double.POSITIVE_INFINITY;
			double minDistSq = Double.POSITIVE_INFINITY;
			List<Coordinate> minima = new ArrayList<>();

			for (Coordinate it : items) {
				double dSq = q.distanceSq(it);
				double d = Math.sqrt(dSq);

				if (d + EPS < minDist) {
					minDist = d;
					minDistSq = dSq;
					minima.clear();
					minima.add(it);
				} else if (Math.abs(d - minDist) <= EPS) {
					// tie under linear distance
					minima.add(it);
					// keep minDistSq consistent with minDist (doesn't matter for assertions below)
				}
			}

			Coordinate found = tree.nearestNeighbor(q, distFn);
			Coordinate foundSq = tree.nearestNeighborSq(q, distSqFn);

			// Non-empty tree should always return a nearest item
			assertNotNull(found, "Tree returned null for non-empty dataset (linear, query " + qi + ")");
			assertNotNull(foundSq, "Tree returned null for non-empty dataset (squared, query " + qi + ")");

			// Linear: found item's distance must equal brute-force minimum (within EPS),
			// and be one of minima
			double foundDist = distFn.distance(q, found);
			assertEquals(minDist, foundDist, EPS, "Linear distance mismatch for query " + qi);
			assertTrue(minima.contains(found), "Linear nearest is not one of brute-force minima for query " + qi);

			// Squared: found item's squared distance must equal brute-force minimum squared
			// (within tolerance)
			double foundDistSq = distSqFn.distance(q, foundSq);
			assertEquals(minDistSq, foundDistSq, 2.0 * EPS * minDist + EPS * EPS, "Squared distance mismatch for query " + qi);
			// Also validate squared result is consistent with linear minimum (tie-safe)
			assertEquals(minDist, Math.sqrt(foundDistSq), EPS, "Squared nearest not consistent with linear min for query " + qi);
		}

		// Also test empty tree returns null for both variants
		HPRtreeX<Coordinate> emptyTree = new HPRtreeX<>();
		assertNull(emptyTree.nearestNeighbor(new Coordinate(0, 0), distFn), "Empty tree (linear) must return null");
		assertNull(emptyTree.nearestNeighborSq(new Coordinate(0, 0), distSqFn), "Empty tree (squared) must return null");
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

	@Test
	public void testRangeQueryWithCoordinates() {
		final int N_ITEMS = 5000;
		final int N_QUERIES = 200;
		Random rnd = new Random(4242);

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

		tree.build();

		// distance functions
		HPRtreeX.DistanceToItem<Coordinate> distFn = (q, item) -> q.distance(item);
		HPRtreeX.DistanceToItem<Coordinate> distSqFn = (q, item) -> q.distanceSq(item);

		// Random queries with random radii
		for (int qi = 0; qi < N_QUERIES; qi++) {
			Coordinate q = new Coordinate(rnd.nextDouble() * 1000.0, rnd.nextDouble() * 1000.0);
			double r = rnd.nextDouble() * 200.0;
			double r2 = r * r;

			List<Coordinate> expected = new ArrayList<>();
			for (Coordinate it : items) {
				if (distFn.distance(q, it) <= r) {
					expected.add(it);
				}
			}

			List<Coordinate> actual = tree.rangeQuery(q, r, distFn);
			List<Coordinate> actualSq = tree.rangeQuerySq(q, r2, distSqFn);

			// Set equality (order-independent) for both variants
			assertEquals(expected.size(), actual.size(), "Linear result size mismatch for query " + qi + " with r=" + r);
			assertTrue(expected.containsAll(actual) && actual.containsAll(expected), "Linear set mismatch for query " + qi + " with r=" + r);

			assertEquals(expected.size(), actualSq.size(), "Squared result size mismatch for query " + qi + " with rSq=" + r2);
			assertTrue(expected.containsAll(actualSq) && actualSq.containsAll(expected), "Squared set mismatch for query " + qi + " with rSq=" + r2);

			// Sanity: returned items satisfy their respective thresholds
			for (Coordinate it : actual) {
				assertTrue(distFn.distance(q, it) <= r + EPS, "Linear returned item beyond radius for query " + qi);
			}
			for (Coordinate it : actualSq) {
				assertTrue(distSqFn.distance(q, it) <= r2 + (EPS * EPS), "Squared returned item beyond radius for query " + qi);
			}
		}

		// Zero radius at an existing item
		Coordinate q0 = items.get(rnd.nextInt(items.size()));
		double r0 = 0.0;
		double r0Sq = 0.0;

		List<Coordinate> expected0 = new ArrayList<>();
		for (Coordinate it : items) {
			if (distFn.distance(q0, it) <= r0) {
				expected0.add(it);
			}
		}

		List<Coordinate> actual0 = tree.rangeQuery(q0, r0, distFn);
		List<Coordinate> actual0Sq = tree.rangeQuerySq(q0, r0Sq, distSqFn);

		assertEquals(expected0.size(), actual0.size());
		assertTrue(expected0.containsAll(actual0) && actual0.containsAll(expected0));

		assertEquals(expected0.size(), actual0Sq.size());
		assertTrue(expected0.containsAll(actual0Sq) && actual0Sq.containsAll(expected0));

		// Very large radius: should return all items
		double hugeR = 1e9;
		double hugeRSq = hugeR * hugeR;

		List<Coordinate> allHits = tree.rangeQuery(new Coordinate(500.0, 500.0), hugeR, distFn);
		List<Coordinate> allHitsSq = tree.rangeQuerySq(new Coordinate(500.0, 500.0), hugeRSq, distSqFn);

		assertEquals(items.size(), allHits.size());
		assertTrue(allHits.containsAll(items) && items.containsAll(allHits));

		assertEquals(items.size(), allHitsSq.size());
		assertTrue(allHitsSq.containsAll(items) && items.containsAll(allHitsSq));
	}

	@Test
	public void testRangeQueryWithLineSegments() {
		final int N_ITEMS = 4000;
		final int N_QUERIES = 150;
		Random rnd = new Random(2121);

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

		// Random queries with random radii
		for (int qi = 0; qi < N_QUERIES; qi++) {
			Coordinate q = new Coordinate(rnd.nextDouble() * 1000.0, rnd.nextDouble() * 1000.0);
			double r = rnd.nextDouble() * 200.0;

			List<LineSegment> expected = new ArrayList<>();
			for (LineSegment seg : items) {
				if (distFn.distance(q, seg) <= r) {
					expected.add(seg);
				}
			}

			List<LineSegment> actual = tree.rangeQuery(q, r, distFn);

			assertEquals(expected.size(), actual.size(), "Result size mismatch for query " + qi + " with r=" + r);
			assertTrue(expected.containsAll(actual) && actual.containsAll(expected), "Set mismatch for query " + qi + " with r=" + r);

			for (LineSegment it : actual) {
				assertTrue(distFn.distance(q, it) <= r + EPS, "Returned item beyond radius for query " + qi);
			}
		}

		// Zero radius at a segment endpoint should include that segment
		if (!items.isEmpty()) {
			LineSegment seg = items.get(rnd.nextInt(items.size()));
			Coordinate endpoint = seg.p0;
			List<LineSegment> expected0 = new ArrayList<>();
			for (LineSegment s : items) {
				if (distFn.distance(endpoint, s) <= 0.0) {
					expected0.add(s);
				}
			}
			List<LineSegment> actual0 = tree.rangeQuery(endpoint, 0.0, distFn);
			assertEquals(expected0.size(), actual0.size());
			assertTrue(expected0.containsAll(actual0) && actual0.containsAll(expected0));
		}

		// Very large radius: should return all items
		double hugeR = 1e9;
		List<LineSegment> allHits = tree.rangeQuery(new Coordinate(250.0, 750.0), hugeR, distFn);
		assertEquals(items.size(), allHits.size());
		assertTrue(allHits.containsAll(items) && items.containsAll(allHits));
	}

	@Test
	public void testRangeQueryEdgeCases() {
		HPRtreeX<Coordinate> tree = new HPRtreeX<>(8);
		HPRtreeX.DistanceToItem<Coordinate> distFn = (q, item) -> q.distance(item);

		// Empty tree
		List<Coordinate> resEmpty = tree.rangeQuery(new Coordinate(0, 0), 10.0, distFn);
		assertTrue(resEmpty.isEmpty(), "Empty tree should return empty list");

		// Insert a couple of points
		tree = new HPRtreeX<>(8);
		Coordinate a = new Coordinate(0, 0);
		Coordinate b = new Coordinate(3, 4); // distance 5 from origin
		tree.insert(new Envelope(a, a), a);
		tree.insert(new Envelope(b, b), b);
		tree.build();

		// Negative radius -> empty
		List<Coordinate> neg = tree.rangeQuery(new Coordinate(0, 0), -1.0, distFn);
		assertTrue(neg.isEmpty(), "Negative radius should return empty list");

		// NaN radius -> empty
		List<Coordinate> nan = tree.rangeQuery(new Coordinate(0, 0), Double.NaN, distFn);
		assertTrue(nan.isEmpty(), "NaN radius should return empty list");

		// Small radius: only 'a'
		List<Coordinate> small = tree.rangeQuery(new Coordinate(0, 0), 1.0, distFn);
		assertEquals(1, small.size());
		assertTrue(small.contains(a));

		// Radius 5 should include both a and b (b is exactly 5 away)
		List<Coordinate> r5 = tree.rangeQuery(new Coordinate(0, 0), 5.0, distFn);
		assertEquals(2, r5.size());
		assertTrue(r5.contains(a) && r5.contains(b));

		// Huge radius -> all
		List<Coordinate> all = tree.rangeQuery(new Coordinate(100, 100), 1e6, distFn);
		assertEquals(2, all.size());
		assertTrue(all.contains(a) && all.contains(b));
	}
}