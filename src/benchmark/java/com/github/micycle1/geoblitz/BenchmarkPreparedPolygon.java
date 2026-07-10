package com.github.micycle1.geoblitz;

import java.util.SplittableRandom;
import java.util.concurrent.TimeUnit;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.operation.relateng.RelateNG;
import org.locationtech.jts.operation.relateng.RelatePredicate;
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

/**
 * Benchmarks FastPreparedPolygon against JTS PreparedGeometry (PreparedPolygon)
 * and RelateNG.prepare() for repeated spatial predicates against a fixed
 * complex polygon.
 */
@State(Scope.Benchmark)
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 4, time = 1)
@Fork(value = 1, jvmArgsAppend = { "-Xms1g", "-Xmx1g", "-XX:+AlwaysPreTouch" })
public class BenchmarkPreparedPolygon {

	private static final int NUM_POINTS = 10_000;
	private static final int NUM_LINES = 1_000;
	private static final int NUM_POLYS = 1_000;

	// number of input points for the concave-hull target polygon (parameterised)
	@Param({ "1000", "10000", "100000" })
	public int n;

	private Geometry polygon;
	private Geometry[] testPoints;
	private Geometry[] testLines;
	private Geometry[] testPolys;

	private PreparedGeometry jtsPrep;
	private FastPreparedPolygon fastPrep;
	private RelateNG relateNg;

	@Setup(Level.Trial)
	public void setup() {
		final long seed = 42L;
		final GeometryFactory gf = new GeometryFactory();

		polygon = GeomMaker.make(n, seed);

		jtsPrep = PreparedGeometryFactory.prepare(polygon);
		fastPrep = new FastPreparedPolygon((Polygonal) polygon);
		relateNg = RelateNG.prepare(polygon);

		final SplittableRandom rnd = new SplittableRandom(seed);
		final Envelope env = polygon.getEnvelopeInternal();
		final double w = env.getWidth(), h = env.getHeight();
		final double minX = env.getMinX() - 0.1 * w, maxX = env.getMaxX() + 0.1 * w;
		final double minY = env.getMinY() - 0.1 * h, maxY = env.getMaxY() + 0.1 * h;

		testPoints = new Geometry[NUM_POINTS];
		for (int i = 0; i < NUM_POINTS; i++) {
			testPoints[i] = gf.createPoint(new Coordinate(rnd.nextDouble(minX, maxX), rnd.nextDouble(minY, maxY)));
		}

		testLines = new Geometry[NUM_LINES];
		for (int i = 0; i < NUM_LINES; i++) {
			int nv = 2 + rnd.nextInt(4);
			Coordinate[] cs = new Coordinate[nv];
			double x = rnd.nextDouble(minX, maxX), y = rnd.nextDouble(minY, maxY);
			for (int v = 0; v < nv; v++) {
				cs[v] = new Coordinate(x, y);
				x += rnd.nextDouble(-0.05, 0.05) * w;
				y += rnd.nextDouble(-0.05, 0.05) * h;
			}
			testLines[i] = gf.createLineString(cs);
		}

		testPolys = new Geometry[NUM_POLYS];
		for (int i = 0; i < NUM_POLYS; i++) {
			double x = rnd.nextDouble(minX, maxX), y = rnd.nextDouble(minY, maxY);
			double r = rnd.nextDouble(0.005, 0.05) * Math.min(w, h);
			testPolys[i] = gf.createPoint(new Coordinate(x, y)).buffer(r, 8);
		}

		// force lazy index construction outside the measured region
		jtsPrep.intersects(testPoints[0]);
		fastPrep.intersects(testPoints[0]);
		relateNg.evaluate(testPoints[0], RelatePredicate.intersects());
	}

	// ---------- intersects ----------

	@Benchmark
	@OperationsPerInvocation(NUM_POINTS)
	public void intersectsPointsJts(Blackhole bh) {
		for (Geometry g : testPoints) {
			bh.consume(jtsPrep.intersects(g));
		}
	}

	@Benchmark
	@OperationsPerInvocation(NUM_POINTS)
	public void intersectsPointsFast(Blackhole bh) {
		for (Geometry g : testPoints) {
			bh.consume(fastPrep.intersects(g));
		}
	}

	@Benchmark
	@OperationsPerInvocation(NUM_POINTS)
	public void intersectsPointsRelateNG(Blackhole bh) {
		for (Geometry g : testPoints) {
			bh.consume(relateNg.evaluate(g, RelatePredicate.intersects()));
		}
	}

	@Benchmark
	@OperationsPerInvocation(NUM_LINES)
	public void intersectsLinesJts(Blackhole bh) {
		for (Geometry g : testLines) {
			bh.consume(jtsPrep.intersects(g));
		}
	}

	@Benchmark
	@OperationsPerInvocation(NUM_LINES)
	public void intersectsLinesFast(Blackhole bh) {
		for (Geometry g : testLines) {
			bh.consume(fastPrep.intersects(g));
		}
	}

	@Benchmark
	@OperationsPerInvocation(NUM_LINES)
	public void intersectsLinesRelateNG(Blackhole bh) {
		for (Geometry g : testLines) {
			bh.consume(relateNg.evaluate(g, RelatePredicate.intersects()));
		}
	}

	@Benchmark
	@OperationsPerInvocation(NUM_POLYS)
	public void intersectsPolysJts(Blackhole bh) {
		for (Geometry g : testPolys) {
			bh.consume(jtsPrep.intersects(g));
		}
	}

	@Benchmark
	@OperationsPerInvocation(NUM_POLYS)
	public void intersectsPolysFast(Blackhole bh) {
		for (Geometry g : testPolys) {
			bh.consume(fastPrep.intersects(g));
		}
	}

	@Benchmark
	@OperationsPerInvocation(NUM_POLYS)
	public void intersectsPolysRelateNG(Blackhole bh) {
		for (Geometry g : testPolys) {
			bh.consume(relateNg.evaluate(g, RelatePredicate.intersects()));
		}
	}

	// ---------- contains ----------

	@Benchmark
	@OperationsPerInvocation(NUM_LINES)
	public void containsLinesJts(Blackhole bh) {
		for (Geometry g : testLines) {
			bh.consume(jtsPrep.contains(g));
		}
	}

	@Benchmark
	@OperationsPerInvocation(NUM_LINES)
	public void containsLinesFast(Blackhole bh) {
		for (Geometry g : testLines) {
			bh.consume(fastPrep.contains(g));
		}
	}

	@Benchmark
	@OperationsPerInvocation(NUM_LINES)
	public void containsLinesRelateNG(Blackhole bh) {
		for (Geometry g : testLines) {
			bh.consume(relateNg.evaluate(g, RelatePredicate.contains()));
		}
	}

	@Benchmark
	@OperationsPerInvocation(NUM_POLYS)
	public void containsPolysJts(Blackhole bh) {
		for (Geometry g : testPolys) {
			bh.consume(jtsPrep.contains(g));
		}
	}

	@Benchmark
	@OperationsPerInvocation(NUM_POLYS)
	public void containsPolysFast(Blackhole bh) {
		for (Geometry g : testPolys) {
			bh.consume(fastPrep.contains(g));
		}
	}

	@Benchmark
	@OperationsPerInvocation(NUM_POLYS)
	public void containsPolysRelateNG(Blackhole bh) {
		for (Geometry g : testPolys) {
			bh.consume(relateNg.evaluate(g, RelatePredicate.contains()));
		}
	}

	// ---------- build (prepare + first query, since all impls build lazily)
	// ----------

	@Benchmark
	@BenchmarkMode(Mode.AverageTime)
	@OutputTimeUnit(TimeUnit.MICROSECONDS)
	public boolean buildJts() {
		PreparedGeometry p = PreparedGeometryFactory.prepare(polygon);
		return p.intersects(testPoints[0]) | p.intersects(testLines[0]);
	}

	@Benchmark
	@BenchmarkMode(Mode.AverageTime)
	@OutputTimeUnit(TimeUnit.MICROSECONDS)
	public boolean buildFast() {
		FastPreparedPolygon p = new FastPreparedPolygon((Polygonal) polygon);
		return p.intersects(testPoints[0]) | p.intersects(testLines[0]);
	}

	@Benchmark
	@BenchmarkMode(Mode.AverageTime)
	@OutputTimeUnit(TimeUnit.MICROSECONDS)
	public boolean buildRelateNG() {
		RelateNG p = RelateNG.prepare(polygon);
		return p.evaluate(testPoints[0], RelatePredicate.intersects()) | p.evaluate(testLines[0], RelatePredicate.intersects());
	}
}
