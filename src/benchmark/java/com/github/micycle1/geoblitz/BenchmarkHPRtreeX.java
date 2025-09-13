package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.List;
import java.util.SplittableRandom;
import java.util.concurrent.TimeUnit;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.LineSegment;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OperationsPerInvocation;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

@State(Scope.Benchmark)
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 0, time = 1)
@Measurement(iterations = 1, time = 2)
@Fork(value = 1, jvmArgsAppend = { "-Xms1g", "-Xmx1g", "-XX:+AlwaysPreTouch" })
public class BenchmarkHPRtreeX {

	// number of items in the dataset (parameterised)
	@Param({ "10000" })
	public int n;

	// range radius for rangeQuery (parameterised)
	@Param({ "20", "100" })
	public double radius;

	private static final int numQueries = 5000;
	private static final Envelope dataEnv = new Envelope(0, 1000, 0, 1000);

	// Points dataset
	private HPRtreeX<Coordinate> treePts;
	private List<Coordinate> itemsPts;
	private Coordinate[] queryPts;

	// Segments dataset
	private HPRtreeX<LineSegment> treeSegs;
	private List<LineSegment> itemsSegs;
	private Coordinate[] queryPtsForSegs;

	// Distance functions
	private HPRtreeX.DistanceToItem<Coordinate> pointDistFn;
	private HPRtreeX.DistanceToItem<LineSegment> segmentDistFn;

	@Setup(Level.Trial)
	public void setup() {
		final long seed = 42L;
		final SplittableRandom rnd = new SplittableRandom(seed);

		// Points
		itemsPts = new ArrayList<>(n);
		treePts = new HPRtreeX<>(16);
		for (int i = 0; i < n; i++) {
			double x = rnd.nextDouble(dataEnv.getMinX(), dataEnv.getMaxX());
			double y = rnd.nextDouble(dataEnv.getMinY(), dataEnv.getMaxY());
			Coordinate c = new Coordinate(x, y);
			itemsPts.add(c);
			treePts.insert(new Envelope(x, x, y, y), c);
		}
		treePts.build();

		// Segments
		itemsSegs = new ArrayList<>(n);
		treeSegs = new HPRtreeX<>(16);
		for (int i = 0; i < n; i++) {
			double x1 = rnd.nextDouble(dataEnv.getMinX(), dataEnv.getMaxX());
			double y1 = rnd.nextDouble(dataEnv.getMinY(), dataEnv.getMaxY());
			double x2 = rnd.nextDouble(dataEnv.getMinX(), dataEnv.getMaxX());
			double y2 = rnd.nextDouble(dataEnv.getMinY(), dataEnv.getMaxY());
			Coordinate p0 = new Coordinate(x1, y1);
			Coordinate p1 = new Coordinate(x2, y2);
			LineSegment seg = new LineSegment(p0, p1);
			itemsSegs.add(seg);
			treeSegs.insert(new Envelope(p0, p1), seg);
		}
		treeSegs.build();

		// Queries
		queryPts = new Coordinate[numQueries];
		queryPtsForSegs = new Coordinate[numQueries];
		for (int i = 0; i < numQueries; i++) {
			double x = rnd.nextDouble(dataEnv.getMinX(), dataEnv.getMaxX());
			double y = rnd.nextDouble(dataEnv.getMinY(), dataEnv.getMaxY());
			queryPts[i] = new Coordinate(x, y);

			double xs = rnd.nextDouble(dataEnv.getMinX(), dataEnv.getMaxX());
			double ys = rnd.nextDouble(dataEnv.getMinY(), dataEnv.getMaxY());
			queryPtsForSegs[i] = new Coordinate(xs, ys);
		}

		// Distance functions
		pointDistFn = (q, item) -> q.distance(item);
		segmentDistFn = (q, seg) -> seg.distance(q);

		// Touch once to ensure any lazy paths are exercised before measurement
		if (!itemsPts.isEmpty()) {
			treePts.nearestNeighbor(queryPts[0], pointDistFn);
			treePts.rangeQuery(queryPts[0], radius, pointDistFn);
		}
		if (!itemsSegs.isEmpty()) {
			treeSegs.nearestNeighbor(queryPtsForSegs[0], segmentDistFn);
			treeSegs.rangeQuery(queryPtsForSegs[0], radius, segmentDistFn);
		}
	}

	// ----------------------------
	// Nearest neighbor: Points
	// ----------------------------

	@Benchmark
	@OperationsPerInvocation(numQueries)
	public void b0nnTreePoints(Blackhole bh) {
		for (int i = 0; i < numQueries; i++) {
			bh.consume(treePts.nearestNeighbor(queryPts[i], pointDistFn));
		}
	}

	@Benchmark
	@OperationsPerInvocation(numQueries)
	public void b0nnBruteForcePoints(Blackhole bh) {
		for (int i = 0; i < numQueries; i++) {
			Coordinate q = queryPts[i];
			Coordinate best = null;
			double bestDist = Double.POSITIVE_INFINITY;
			for (int j = 0; j < itemsPts.size(); j++) {
				Coordinate it = itemsPts.get(j);
				double d = q.distance(it);
				if (d < bestDist) {
					bestDist = d;
					best = it;
					if (bestDist == 0.0)
						break;
				}
			}
			bh.consume(best);
		}
	}

	// ----------------------------
	// Range query: Points
	// ----------------------------

	@Benchmark
	@OperationsPerInvocation(numQueries)
	public void b1rangeTreePoints(Blackhole bh) {
		for (int i = 0; i < numQueries; i++) {
			bh.consume(treePts.rangeQuery(queryPts[i], radius, pointDistFn));
		}
	}

	@Benchmark
	@OperationsPerInvocation(numQueries)
	public void b1rangeBruteForcePoints(Blackhole bh) {
		for (int i = 0; i < numQueries; i++) {
			Coordinate q = queryPts[i];
			List<Coordinate> res = new ArrayList<>();
			for (int j = 0; j < itemsPts.size(); j++) {
				Coordinate it = itemsPts.get(j);
				if (q.distance(it) <= radius) {
					res.add(it);
				}
			}
			bh.consume(res);
		}
	}

	// ----------------------------
	// Nearest neighbor: Segments
	// ----------------------------

	@Benchmark
	@OperationsPerInvocation(numQueries)
	public void b2nnTreeSegments(Blackhole bh) {
		for (int i = 0; i < numQueries; i++) {
			bh.consume(treeSegs.nearestNeighbor(queryPtsForSegs[i], segmentDistFn));
		}
	}

	@Benchmark
	@OperationsPerInvocation(numQueries)
	public void b2nnBruteForceSegments(Blackhole bh) {
		for (int i = 0; i < numQueries; i++) {
			Coordinate q = queryPtsForSegs[i];
			LineSegment best = null;
			double bestDist = Double.POSITIVE_INFINITY;
			for (int j = 0; j < itemsSegs.size(); j++) {
				LineSegment seg = itemsSegs.get(j);
				double d = seg.distance(q);
				if (d < bestDist) {
					bestDist = d;
					best = seg;
					if (bestDist == 0.0)
						break;
				}
			}
			bh.consume(best);
		}
	}

	// ----------------------------
	// Range query: Segments
	// ----------------------------

	@Benchmark
	@OperationsPerInvocation(numQueries)
	public void b3rangeTreeSegments(Blackhole bh) {
		for (int i = 0; i < numQueries; i++) {
			bh.consume(treeSegs.rangeQuery(queryPtsForSegs[i], radius, segmentDistFn));
		}
	}

	@Benchmark
	@OperationsPerInvocation(numQueries)
	public void b3rangeBruteForceSegments(Blackhole bh) {
		for (int i = 0; i < numQueries; i++) {
			Coordinate q = queryPtsForSegs[i];
			List<LineSegment> res = new ArrayList<>();
			for (int j = 0; j < itemsSegs.size(); j++) {
				LineSegment seg = itemsSegs.get(j);
				if (seg.distance(q) <= radius) {
					res.add(seg);
				}
			}
			bh.consume(res);
		}
	}
}