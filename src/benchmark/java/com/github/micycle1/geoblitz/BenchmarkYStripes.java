package com.github.micycle1.geoblitz;

import java.util.SplittableRandom;
import java.util.concurrent.TimeUnit;

import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
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
@Warmup(iterations = 5, time = 1)
@Measurement(iterations = 5, time = 1)
@Fork(value = 1, jvmArgsAppend = { "-Xms1g", "-Xmx1g", "-XX:+AlwaysPreTouch" })
public class BenchmarkYStripes {

	private Polygon polygon;
	private Polygon omniPolygon;
	private Coordinate[] testPoints;
	private YStripesPointInAreaLocator yStripesLocator;
	private IndexedPointInAreaLocator indexedLocator;

	private static final int numPoints = 50_000;

	// number of vertices in the polygon (parameterised)
	@Param({ "100", "1000", "10000", "100000" })
	public int n;

	@Setup(Level.Trial)
	public void setup() throws Exception {
		final long seed = 42L;

		// Build test polygon once per trial (per @Param n)
		polygon = (Polygon) GeomMaker.make(n, seed);
//		omniPolygon = ;
		// TODO test
		var q = GeometryFactory.toPolygonArray(GeomMaker.makeRandomBufferedPoints(500, 5, 8, 0, new Envelope(0, 1000, 0, 1000)));

		// Build locators here so benchmarks measure locate(), not construction
		yStripesLocator = new YStripesPointInAreaLocator(polygon);
		indexedLocator = new IndexedPointInAreaLocator(polygon);

		// Generate reproducible test points near the polygon's envelope
		testPoints = new Coordinate[numPoints];
		final SplittableRandom rnd = new SplittableRandom(seed);

		final Envelope env = polygon.getEnvelopeInternal();
		final double w = Math.max(1e-9, env.getWidth());
		final double h = Math.max(1e-9, env.getHeight());

		// Sample in a slightly expanded box to include inside/outside cases
		final double minX = env.getMinX() - 0.1 * w;
		final double maxX = env.getMaxX() + 0.1 * w;
		final double minY = env.getMinY() - 0.1 * h;
		final double maxY = env.getMaxY() + 0.1 * h;

		for (int i = 0; i < numPoints; i++) {
			final double x = rnd.nextDouble(minX, maxX);
			final double y = rnd.nextDouble(minY, maxY);
			testPoints[i] = new Coordinate(x, y);
		}

		// Ensure any lazy initialization in locators happens during warmup/setup, not
		// measurement
		yStripesLocator.locate(testPoints[0]);
		indexedLocator.locate(testPoints[0]);
	}

	@Benchmark
	@OperationsPerInvocation(numPoints)
	public void testYStripesLocator(Blackhole bh) {
		for (Coordinate c : testPoints) {
			bh.consume(yStripesLocator.locate(c));
		}
	}

	@Benchmark
	@OperationsPerInvocation(numPoints)
	public void testIndexedLocator(Blackhole bh) {
		for (Coordinate c : testPoints) {
			bh.consume(indexedLocator.locate(c));
		}
	}
}