package com.github.micycle1.geoblitz;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertFalse;

import java.util.ArrayList;
import java.util.List;
import java.util.SplittableRandom;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.operation.relateng.RelateNG;
import org.locationtech.jts.operation.relateng.RelatePredicate;

/**
 * Verifies that FastPreparedPolygon returns identical results to JTS
 * PreparedGeometry and RelateNG across predicates and geometry types,
 * including exact boundary-touching cases.
 */
class FastPreparedPolygonTest {

	private static final GeometryFactory GF = new GeometryFactory();

	@Test
	void testAgainstJtsOnStarPolygon() {
		Polygon target = starPolygon(500, 500, 400, 64, 1);
		compareAll(target, testGeometries(target, 7));
	}

	@Test
	void testAgainstJtsOnPolygonWithHole() {
		Polygon shell = starPolygon(500, 500, 450, 48, 2);
		Polygon hole = starPolygon(500, 500, 150, 32, 3);
		Geometry target = shell.difference(hole);
		assertTrue(target instanceof Polygonal);
		List<Geometry> tests = testGeometries(target, 11);
		// point inside the hole and a small polygon inside the hole
		tests.add(GF.createPoint(new Coordinate(500, 500)));
		tests.add(GF.createPoint(new Coordinate(500, 500)).buffer(20, 8));
		compareAll(target, tests);
	}

	@Test
	void testAgainstJtsOnMultiPolygon() {
		Polygon a = starPolygon(300, 300, 250, 40, 4);
		Polygon b = starPolygon(1200, 900, 300, 40, 5);
		Geometry target = GF.createMultiPolygon(new Polygon[] { a, b });
		compareAll(target, testGeometries(target, 13));
	}

	@Test
	void testExactBoundaryCases() {
		Polygon target = starPolygon(500, 500, 400, 32, 6);
		List<Geometry> tests = new ArrayList<>();
		// the polygon itself
		tests.add(target.copy());
		// its shell as a line
		tests.add(target.getExteriorRing().copy());
		// a vertex point and an edge midpoint point
		Coordinate[] ring = target.getExteriorRing().getCoordinates();
		tests.add(GF.createPoint(ring[0]));
		tests.add(GF.createPoint(new Coordinate((ring[0].x + ring[1].x) / 2, (ring[0].y + ring[1].y) / 2)));
		// a single boundary edge as a line
		tests.add(GF.createLineString(new Coordinate[] { ring[0], ring[1] }));
		// sub-polygon sharing part of the boundary (triangle on first edge + centroid)
		Coordinate c = target.getCentroid().getCoordinate();
		tests.add(GF.createPolygon(new Coordinate[] { ring[0], ring[1], c, ring[0] }));
		compareAll(target, tests);
	}

	@Test
	void testPrepareFactoryFallback() {
		LineString line = GF.createLineString(new Coordinate[] { new Coordinate(0, 0), new Coordinate(1, 1) });
		PreparedGeometry prep = FastPreparedPolygon.prepare(line);
		assertFalse(prep instanceof FastPreparedPolygon);
		assertTrue(prep.intersects(GF.createPoint(new Coordinate(0.5, 0.5))));

		Polygon poly = starPolygon(0, 0, 10, 8, 1);
		assertTrue(FastPreparedPolygon.prepare(poly) instanceof FastPreparedPolygon);
	}

	private static void compareAll(Geometry target, List<Geometry> tests) {
		FastPreparedPolygon fast = new FastPreparedPolygon((Polygonal) target);
		PreparedGeometry jts = PreparedGeometryFactory.prepare(target);
		RelateNG ng = RelateNG.prepare(target);

		for (Geometry t : tests) {
			String wkt = t.toString();
			assertEquals(jts.intersects(t), fast.intersects(t), () -> "intersects " + wkt);
			assertEquals(jts.contains(t), fast.contains(t), () -> "contains " + wkt);
			assertEquals(jts.containsProperly(t), fast.containsProperly(t), () -> "containsProperly " + wkt);
			assertEquals(jts.covers(t), fast.covers(t), () -> "covers " + wkt);
			assertEquals(jts.disjoint(t), fast.disjoint(t), () -> "disjoint " + wkt);

			assertEquals(ng.evaluate(t, RelatePredicate.intersects()), fast.intersects(t), () -> "NG intersects " + wkt);
			assertEquals(ng.evaluate(t, RelatePredicate.contains()), fast.contains(t), () -> "NG contains " + wkt);
			assertEquals(ng.evaluate(t, RelatePredicate.covers()), fast.covers(t), () -> "NG covers " + wkt);
		}
	}

	/**
	 * Builds a mixed bag of test geometries in and around the target's envelope:
	 * points, multi-vertex linestrings and small polygons, spanning
	 * inside/outside/crossing configurations.
	 */
	private static List<Geometry> testGeometries(Geometry target, long seed) {
		SplittableRandom rnd = new SplittableRandom(seed);
		Envelope env = target.getEnvelopeInternal();
		double w = env.getWidth(), h = env.getHeight();
		double minX = env.getMinX() - 0.15 * w, maxX = env.getMaxX() + 0.15 * w;
		double minY = env.getMinY() - 0.15 * h, maxY = env.getMaxY() + 0.15 * h;

		List<Geometry> out = new ArrayList<>();
		for (int i = 0; i < 300; i++) {
			out.add(GF.createPoint(new Coordinate(rnd.nextDouble(minX, maxX), rnd.nextDouble(minY, maxY))));
		}
		for (int i = 0; i < 150; i++) {
			int nv = 2 + rnd.nextInt(4);
			Coordinate[] cs = new Coordinate[nv];
			double x = rnd.nextDouble(minX, maxX), y = rnd.nextDouble(minY, maxY);
			for (int v = 0; v < nv; v++) {
				cs[v] = new Coordinate(x, y);
				x += rnd.nextDouble(-0.1, 0.1) * w;
				y += rnd.nextDouble(-0.1, 0.1) * h;
			}
			out.add(GF.createLineString(cs));
		}
		for (int i = 0; i < 150; i++) {
			double x = rnd.nextDouble(minX, maxX), y = rnd.nextDouble(minY, maxY);
			double r = rnd.nextDouble(0.005, 0.1) * Math.min(w, h);
			out.add(GF.createPoint(new Coordinate(x, y)).buffer(r, 8));
		}
		return out;
	}

	/** Simple (non-self-intersecting) concave star polygon with random radii. */
	private static Polygon starPolygon(double cx, double cy, double maxR, int nArms, long seed) {
		SplittableRandom rnd = new SplittableRandom(seed);
		int n = nArms * 2;
		Coordinate[] cs = new Coordinate[n + 1];
		for (int i = 0; i < n; i++) {
			double angle = 2 * Math.PI * i / n;
			double r = (i % 2 == 0) ? maxR : maxR * rnd.nextDouble(0.3, 0.7);
			cs[i] = new Coordinate(cx + r * Math.cos(angle), cy + r * Math.sin(angle));
		}
		cs[n] = new Coordinate(cs[0]);
		return GF.createPolygon(cs);
	}
}
