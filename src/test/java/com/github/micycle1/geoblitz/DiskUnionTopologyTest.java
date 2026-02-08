package com.github.micycle1.geoblitz;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.valid.IsValidOp;

import com.github.micycle1.geoblitz.DiskUnion.Arc;
import com.github.micycle1.geoblitz.DiskUnion.ArcBoundary;
import com.github.micycle1.geoblitz.DiskUnion.ArcCycle;
import com.github.micycle1.geoblitz.DiskUnion.SnapVertex;

public class DiskUnionTopologyTest {

	private static final double EPS = 1e-10;
	private static final double MAX_SEG_LEN = 0.05;
	private static final GeometryFactory GF = new GeometryFactory();

	@Test
	public void testArcCyclesClosedAndVerticesBalanced_randomCluster() {
		Random rnd = new Random(1234);

		for (int t = 0; t < 50; t++) {
			List<Coordinate> circles = new ArrayList<>();
			for (int i = 0; i < 25; i++) {
				double x = rnd.nextDouble() * 10.0;
				double y = rnd.nextDouble() * 10.0;
				double r = 0.5 + rnd.nextDouble() * 1.5;
				circles.add(circle(x, y, r));
			}

			ArcBoundary b = DiskUnion.computeBoundaryArcs(circles, EPS);

			// basic sanity
			for (ArcCycle cy : b.cycles) {
				assertFalse(cy.arcs.isEmpty(), "Cycle must contain at least one arc");
				assertCycleClosed(cy);
			}

			assertVertexBalance(b);
			assertNoDuplicateArcInstances(b);
		}
	}

	@Test
	public void testBoundaryLocalInsideOutsideProperty_randomCluster() {
		Random rnd = new Random(7);

		for (int t = 0; t < 30; t++) {
			List<Coordinate> circles = new ArrayList<>();
			for (int i = 0; i < 20; i++) {
				double x = rnd.nextDouble() * 10.0;
				double y = rnd.nextDouble() * 10.0;
				double r = 0.5 + rnd.nextDouble() * 1.5;
				circles.add(circle(x, y, r));
			}

			ArcBoundary b = DiskUnion.computeBoundaryArcs(circles, EPS);

			// For a selection of arcs: slightly inside should be covered, slightly outside
			// should not.
			int checked = 0;
			for (ArcCycle cy : b.cycles) {
				for (Arc a : cy.arcs) {
					// Skip full-circle arc if start/end coincide; still valid but "outside" test
					// is fine; keep it anyway.
					assertBoundaryArcLocalProperty(a, circles);
					if (++checked >= 50) {
						break;
					}
				}
				if (checked >= 50) {
					break;
				}
			}
		}
	}

	@Test
	public void testDuplicatePair_arcLevelRegression() {
		List<Coordinate> circles = List.of(circle(0, 0, 2), circle(0, 0, 2));

		ArcBoundary b = DiskUnion.computeBoundaryArcs(circles, EPS);

		assertFalse(b.cycles.isEmpty(), "Duplicate pair must not produce empty boundary");

		double area = Math.abs(signedUnionAreaFromArcs(b));
		double expected = Math.PI * 4.0; // pi * r^2, r=2
		assertEquals(expected, area, expected * 1e-9 + 1e-8, "Area from arcs should match a single disk");
	}

	@Test
	public void testMonotonicity_area_arcLevel() {
		List<Coordinate> A = List.of(circle(0, 0, 2), circle(2, 0, 2));
		List<Coordinate> B = List.of(circle(5, 0, 1));

		List<Coordinate> AB = new ArrayList<>();
		AB.addAll(A);
		AB.addAll(B);

		double areaA = Math.abs(signedUnionAreaFromArcs(DiskUnion.computeBoundaryArcs(A, EPS)));
		double areaAB = Math.abs(signedUnionAreaFromArcs(DiskUnion.computeBoundaryArcs(AB, EPS)));

		// Union area should not decrease. Allow a tiny numeric tolerance.
		assertTrue(areaAB + 1e-10 >= areaA, () -> "Area monotonicity violated: areaA=" + areaA + " areaAB=" + areaAB);
	}

	@Test
	public void testPermutationInvariance_area_arcLevel() {
		// include a duplicate to stress invariance
		List<Coordinate> circles = new ArrayList<>(List.of(circle(0, 0, 2), circle(2, 0, 2), circle(1, 2, 1.5), circle(5, 0, 1), circle(5, 0, 1)));

		double baseArea = Math.abs(signedUnionAreaFromArcs(DiskUnion.computeBoundaryArcs(circles, EPS)));

		Random rnd = new Random(99);
		for (int t = 0; t < 50; t++) {
			Collections.shuffle(circles, rnd);
			double a = Math.abs(signedUnionAreaFromArcs(DiskUnion.computeBoundaryArcs(circles, EPS)));
			assertEquals(baseArea, a, Math.max(1e-8, baseArea * 1e-10), "Order should not change arc-derived area");
		}
	}

	@Test
	public void testRingOfDisksTopologyProducesValidPolygon() {
		// Regression: ensure arc stitching + linearization yields valid geometry
		// (no self-intersections) for a ring-shaped union.
		int n = 8;
		double R = 3.0;
		double r = 1.4;
		List<Coordinate> circles = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			double a = (2.0 * Math.PI * i) / n;
			circles.add(circle(R * Math.cos(a), R * Math.sin(a), r));
		}

		ArcBoundary boundary = DiskUnion.computeBoundaryArcs(circles, EPS);
		Geometry g = DiskUnion.toGeometry(boundary, GF, MAX_SEG_LEN);

		IsValidOp validator = new IsValidOp(g);
		String m = validator.isValid() ? "" : validator.getValidationError().getMessage();
		assertTrue(validator.isValid(), "Linearized ring should be valid: " + m);

		// Expect a single polygon with a hole for this configuration.
		assertTrue(g instanceof Polygon, "Expected a single polygon for the ring union");
		Polygon p = (Polygon) g;
		assertEquals(1, p.getNumInteriorRing(), "Ring union should have exactly one hole");
	}

	@Test
	public void testTripleIntersectionSharedVertexProducesValidGeometry() {
		// Three circles share a single intersection point (origin).
		// Stresses endpoint snapping and arc selection around a shared vertex.
		double r = 1.0;
		List<Coordinate> circles = List.of(circle(1.0, 0.0, r), circle(-0.5, Math.sqrt(3.0) * 0.5, r), circle(-0.5, -Math.sqrt(3.0) * 0.5, r));

		ArcBoundary boundary = DiskUnion.computeBoundaryArcs(circles, EPS);
		Geometry g = DiskUnion.toGeometry(boundary, GF, MAX_SEG_LEN);

		IsValidOp validator = new IsValidOp(g);
		String m = validator.isValid() ? "" : validator.getValidationError().getMessage();
		assertTrue(validator.isValid(), "Shared-vertex union should be valid: " + m);
		assertFalse(g.isEmpty(), "Shared-vertex union should be non-empty");
	}

	private static void assertCycleClosed(ArcCycle cy) {
		int n = cy.arcs.size();
		for (int i = 0; i < n; i++) {
			Arc a = cy.arcs.get(i);
			Arc b = cy.arcs.get((i + 1) % n);

			assertNotNull(a.start);
			assertNotNull(a.end);
			assertNotNull(b.start);
			assertNotNull(b.end);

			assertEquals(a.end, b.start, "Arc endpoints must stitch end->start within cycle");
		}
	}

	private static void assertVertexBalance(ArcBoundary b) {
		Map<SnapVertex, int[]> deg = new HashMap<>(); // [out,in]
		for (ArcCycle cy : b.cycles) {
			for (Arc a : cy.arcs) {
				deg.computeIfAbsent(a.start, __ -> new int[2])[0]++; // out
				deg.computeIfAbsent(a.end, __ -> new int[2])[1]++; // in
			}
		}
		for (Map.Entry<SnapVertex, int[]> e : deg.entrySet()) {
			int out = e.getValue()[0];
			int in = e.getValue()[1];
			assertEquals(in, out, "Vertex indegree must equal outdegree on closed boundary graph");
		}
	}

	private static void assertNoDuplicateArcInstances(ArcBoundary b) {
		Set<Arc> seen = Collections.newSetFromMap(new IdentityHashMap<>());
		for (ArcCycle cy : b.cycles) {
			for (Arc a : cy.arcs) {
				assertTrue(seen.add(a), "Same Arc instance appears more than once across cycles");
			}
		}
	}

	private static void assertBoundaryArcLocalProperty(Arc arc, List<Coordinate> inputCircles) {
		// Use a midpoint on the arc (safe even for long sweeps).
		double t = 0.5 * (arc.a0 + arc.a1);
		double cx = arc.circle.c.x;
		double cy = arc.circle.c.y;
		double r = arc.circle.r;

		double px = cx + r * Math.cos(t);
		double py = cy + r * Math.sin(t);

		// Radial direction (unit vector)
		double ux = (px - cx) / r;
		double uy = (py - cy) / r;

		// Choose a small step relative to scale:
		double step = Math.max(1e-7 * Math.max(1.0, r), EPS * 10.0);

		Coordinate inside = new Coordinate(px - ux * step, py - uy * step);
		Coordinate outside = new Coordinate(px + ux * step, py + uy * step);

		assertTrue(isCoveredByAnyDisk(inside, inputCircles, 0.0), "Point just inside boundary must be covered");
		assertFalse(isCoveredByAnyDisk(outside, inputCircles, 0.0), "Point just outside boundary must be uncovered");
	}

	private static boolean isCoveredByAnyDisk(Coordinate p, List<Coordinate> circles, double tol) {
		double px = p.x, py = p.y;
		for (Coordinate c : circles) {
			double dx = px - c.x;
			double dy = py - c.y;
			double rr = c.getZ() + tol;
			if (dx * dx + dy * dy <= rr * rr) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Signed area contribution of the union boundary represented as circular arcs.
	 * Assumes arcs are oriented with union interior on the left (CCW outer, CW
	 * holes).
	 */
	private static double signedUnionAreaFromArcs(ArcBoundary b) {
		double sum = 0.0;
		for (ArcCycle cy : b.cycles) {
			for (Arc a : cy.arcs) {
				sum += signedAreaOfArc(a);
			}
		}
		return sum;
	}

	/**
	 * Exact signed area contribution for a circular arc segment from angle a0 to a1
	 * (CCW). Formula: 1/2 ∫ (x dy - y dx)
	 */
	private static double signedAreaOfArc(Arc a) {
		double cx = a.circle.c.x;
		double cy = a.circle.c.y;
		double r = a.circle.r;

		double t0 = a.a0;
		double t1 = a.a1;
		double dt = t1 - t0;

		// 0.5 * (r^2*dt + cx*r*(sin t1 - sin t0) - cy*r*(cos t1 - cos t0))
		return 0.5 * (r * r * dt + cx * r * (Math.sin(t1) - Math.sin(t0)) - cy * r * (Math.cos(t1) - Math.cos(t0)));
	}

	private static Coordinate circle(double x, double y, double r) {
		return new Coordinate(x, y, r);
	}
}
