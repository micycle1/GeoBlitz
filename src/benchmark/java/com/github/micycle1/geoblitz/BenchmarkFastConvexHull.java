package com.github.micycle1.geoblitz;

import java.util.SplittableRandom;
import java.util.concurrent.TimeUnit;

import org.locationtech.jts.algorithm.ConvexHull;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

@State(Scope.Benchmark)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 1, time = 1)
@Measurement(iterations = 2, time = 2)
@Fork(value = 1)
public class BenchmarkFastConvexHull {

	private ConvexHull jtsHull;
	private FastConvexHull fastHull;
	private Coordinate[] testPoints;
	private GeometryFactory gf;

	@Param({ "10000", "100000", "500000", "1000000" })
	public int n;

	@Setup(Level.Iteration) // Interation-level, in case testpoints are sorted
	public void setup() throws Exception {
//		final long seed = 42L;
		final SplittableRandom rnd = new SplittableRandom(System.currentTimeMillis());
		testPoints = new Coordinate[n];
		for (int i = 0; i < n; i++) {
			final double x = rnd.nextDouble();
			final double y = rnd.nextDouble();
			testPoints[i] = new Coordinate(x, y);
		}

		gf = new GeometryFactory();
	}

	@Benchmark
//	@OperationsPerInvocation(n)
	public void testFastConvexHull(Blackhole bh) {
		fastHull = new FastConvexHull(testPoints, gf);
		bh.consume(fastHull.getConvexHull());
	}

	@Benchmark
//	@OperationsPerInvocation(numPoints)
	public void testJTSConvexHull(Blackhole bh) {
		jtsHull = new ConvexHull(testPoints, gf);
		bh.consume(jtsHull.getConvexHull());
	}
}
