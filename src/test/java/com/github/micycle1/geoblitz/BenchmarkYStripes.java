package com.github.micycle1.geoblitz;

import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.io.ParseException;
import org.openjdk.jmh.annotations.*;

import java.util.Random;
import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@BenchmarkMode(Mode.Throughput)
@Warmup(time = 1)
@Measurement(time = 1)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
public class BenchmarkYStripes {

	private Polygon polygon;
	private Coordinate[] testPoints;
	private YStripesPointInAreaLocator yStripesLocator;
	private IndexedPointInAreaLocator indexedLocator;

	private static final int numPoints = 50000;

	// number of vertices in the polygon (parameterised)
	@Param({ "100", "1000", "10000", "100000" })
	public int n;

	@Setup
	public void setup() throws ParseException {
		long seed = 42;
		polygon = (Polygon) TestGeomMaker.make(n, seed);

		yStripesLocator = new YStripesPointInAreaLocator(polygon);
		indexedLocator = new IndexedPointInAreaLocator(polygon);

		// Generate random test points
		testPoints = new Coordinate[numPoints];
		Random rand = new Random(seed);
		for (int i = 0; i < numPoints; i++) {
			double x = rand.nextDouble(1000);
			double y = rand.nextDouble(1000);
			testPoints[i] = new Coordinate(x, y);
		}
	}

	@Benchmark
	@OperationsPerInvocation(numPoints)
	public int testYStripesLocator() {
		int count = 0;
		for (Coordinate c : testPoints) {
			count += yStripesLocator.locate(c);
		}
		return count;
	}

	@Benchmark
	@OperationsPerInvocation(numPoints)
	public int testIndexedLocator() {
		int count = 0;
		for (Coordinate c : testPoints) {
			count += indexedLocator.locate(c);
		}
		return count;
	}
}
