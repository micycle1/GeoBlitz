package com.github.micycle1.geoblitz;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;

public class CircleIndexTest {

    private static final double EPS = 1e-9;

    @Test
    public void testEmpty() {
        CircleIndex<String> index = new CircleIndex<>();
        assertEquals(0, index.size());
        assertFalse(index.existsInside(0, 0));
        assertNull(index.nearest(0, 0));
    }

    @Test
    public void testSingleItem() {
        CircleIndex<String> index = new CircleIndex<>();
        index.insert(10, 10, 5, "A");

        assertEquals(1, index.size());
        assertEquals(5.0, index.maxRadius(), EPS);

        // Inside
        assertTrue(index.existsInside(10, 10));
        assertTrue(index.existsInside(13, 10));

        // Outside
        assertFalse(index.existsInside(20, 20));

        // Nearest
        CircleIndex.Nearest<String> n = index.nearest(20, 20);
        assertNotNull(n);
        assertEquals("A", n.value);
        assertEquals(10, n.cx, EPS);
        assertEquals(10, n.cy, EPS);
        assertEquals(5, n.cr, EPS);
    }

    @Test
    public void testMultipleItems() {
        CircleIndex<Integer> index = new CircleIndex<>();
        index.insert(0, 0, 1, 1);
        index.insert(10, 0, 1, 2);
        index.insert(0, 10, 1, 3);

        assertEquals(3, index.size());

        // Test nearest
        assertEquals(1, index.nearest(1, 1).value);
        assertEquals(2, index.nearest(9, 1).value);
        assertEquals(3, index.nearest(1, 9).value);

        // Test existsInside
        assertTrue(index.existsInside(0.5, 0.5));
        assertTrue(index.existsInside(10.5, 0));
        assertFalse(index.existsInside(5, 5));
    }

    @Test
    public void testExistsWithinMetric() {
        CircleIndex<String> index = new CircleIndex<>();
        // Circle at (0,0) radius 2
        index.insert(0, 0, 2, "A");

        double R = index.maxRadius(); // 2

        // limit = R => existsInside
        // Point (1,0) is inside (dist 1 < 2)
        assertTrue(index.existsWithinMetric(1, 0, R, R));
        // Point (3,0) is outside (dist 3 > 2)
        assertFalse(index.existsWithinMetric(3, 0, R, R));

        // limit = R + 1 => clearance <= 1
        // Point (2.5, 0) -> dist=2.5, clearance = 0.5 <= 1 -> TRUE
        assertTrue(index.existsWithinMetric(2.5, 0, R, R + 1));

        // Point (4, 0) -> dist=4, clearance = 2 > 1 -> FALSE
        assertFalse(index.existsWithinMetric(4, 0, R, R + 1));
    }

    @Test
    public void testRandomizedVsBruteForce() {
        CircleIndex<Integer> index = new CircleIndex<>();
        List<TestCircle> circles = new ArrayList<>();
        Random rnd = new Random(42);

        int N = 1000;
        for (int i = 0; i < N; i++) {
            double x = rnd.nextDouble() * 100;
            double y = rnd.nextDouble() * 100;
            double r = 1 + rnd.nextDouble() * 4;
            index.insert(x, y, r, i);
            circles.add(new TestCircle(x, y, r, i));
        }

        assertEquals(N, index.size());

        // Query random points
        for (int i = 0; i < 200; i++) {
            double qx = rnd.nextDouble() * 100;
            double qy = rnd.nextDouble() * 100;

            CircleIndex.Nearest<Integer> nearest = index.nearest(qx, qy);
            TestCircle expected = findNearestBruteForce(circles, qx, qy);

            assertNotNull(nearest);

            // Compare clearance = hypot(q-c) - r
            double actualClearance = nearest.clearance;
            double expectedClearance = dist(qx, qy, expected) - expected.r;

            assertEquals(expectedClearance, actualClearance, EPS, "Clearance mismatch at " + i);
        }
    }

    @Test
    public void testContainsWithNested() {
        CircleIndex<String> index = new CircleIndex<>();
        index.insert(0, 0, 10, "Big");
        index.insert(0, 0, 2, "Small");

        assertTrue(index.existsInside(0, 0));
        assertTrue(index.existsInside(5, 0)); // Inside Big, outside Small
    }

    private static class TestCircle {
        double x, y, r;
        int id;

        TestCircle(double x, double y, double r, int id) {
            this.x = x;
            this.y = y;
            this.r = r;
            this.id = id;
        }
    }

    private TestCircle findNearestBruteForce(List<TestCircle> circles, double qx, double qy) {
        TestCircle best = null;
        double minClearance = Double.POSITIVE_INFINITY;

        for (TestCircle c : circles) {
            double d = dist(qx, qy, c) - c.r;
            if (d < minClearance) {
                minClearance = d;
                best = c;
            }
        }
        return best;
    }

    private double dist(double qx, double qy, TestCircle c) {
        return Math.hypot(qx - c.x, qy - c.y);
    }
}
