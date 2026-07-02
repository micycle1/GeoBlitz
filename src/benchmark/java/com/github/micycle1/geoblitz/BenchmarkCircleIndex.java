package com.github.micycle1.geoblitz;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.openjdk.jmh.annotations.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 2, time = 1)
@Measurement(iterations = 3, time = 1)
@Fork(value = 1, jvmArgsAppend = { "-Xms2g", "-Xmx4g" })
public class BenchmarkCircleIndex {

	@Param({ "1000", "10000", "50000" })
	public int nCircles;

	private CircleIndex<Object> circleIndex;
	private HPRtreeX<Object> hprTree;
	private List<Coordinate> queryPoints;
	private static final int QUERY_BATCH = 1000;
	private static final int K = 10;
	private static final Object DUMMY = new Object();

	// Made public to ensure visibility for JMH/reflection if needed
	public static class Circle {
		public double x, y, r;

		public Circle(double x, double y, double r) {
			this.x = x;
			this.y = y;
			this.r = r;
		}
	}

	@Setup(Level.Trial)
	public void setup() {
		circleIndex = new CircleIndex<>();
		hprTree = new HPRtreeX<>();

		Random rand = new Random(1337);
		double EXTENT = 1000.0;

		for (int i = 0; i < nCircles; i++) {
			double x = rand.nextDouble() * EXTENT;
			double y = rand.nextDouble() * EXTENT;
			double r = 1 + rand.nextDouble() * 20; // Radius 1..21

			// Insert into CircleIndex
			circleIndex.insert(x, y, r, DUMMY);

			// Insert into HPRtree
			Circle c = new Circle(x, y, r);
			hprTree.insert(new Envelope(x - r, x + r, y - r, y + r), c);
		}
		hprTree.build();

		// Generate query points
		queryPoints = new ArrayList<>(QUERY_BATCH);
		for (int i = 0; i < QUERY_BATCH; i++) {
			queryPoints.add(new Coordinate(rand.nextDouble() * EXTENT, rand.nextDouble() * EXTENT));
		}
	}

	@Benchmark
	public int circleIndexContains() {
		int hits = 0;
		for (int i = 0; i < QUERY_BATCH; i++) {
			Coordinate q = queryPoints.get(i);
			if (circleIndex.existsInside(q.x, q.y)) {
				hits++;
			}
		}
		return hits;
	}

	@Benchmark
	public int circleIndexKNearest() {
		int total = 0;
		for (int i = 0; i < QUERY_BATCH; i++) {
			Coordinate q = queryPoints.get(i);
			CircleIndex.Nearest<Object>[] res = circleIndex.kNearest(q.x, q.y, K);
			total += res.length;
		}

		return total;
	}

	@Benchmark
	public int hprTreeContains() {
		int hits = 0;
		for (int i = 0; i < QUERY_BATCH; i++) {
			Coordinate q = queryPoints.get(i);
			// Using fully coordinate-based constructor to avoid ambiguity
			Envelope env = new Envelope(q.x, q.x, q.y, q.y);

			// Use an array to capture the result from the lambda
			final boolean[] found = { false };

			// Query the tree for envelopes containing the point
			hprTree.query(env, (Object item) -> { // Explicit type
				Circle c = (Circle) item;
				// Exact circle check
				double dx = q.x - c.x;
				double dy = q.y - c.y;
				double r = c.r;
				if (dx * dx + dy * dy <= r * r) {
					found[0] = true;
					return false; // Stop traversal immediately if found
				}
				return true; // Continue traversal
			});

			if (found[0]) {
				hits++;
			}
		}
		return hits;
	}
}
