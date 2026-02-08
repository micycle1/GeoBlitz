package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.List;
import java.util.SplittableRandom;
import java.util.concurrent.TimeUnit;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.union.CascadedPolygonUnion;
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
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 2, time = 2)
@Measurement(iterations = 2, time = 5)
@Fork(value = 1, jvmArgsAppend = { "-Xms1g", "-Xmx1g", "-XX:+AlwaysPreTouch" })
public class BenchmarkDiskUnion {

	private static final long SEED = 42L;
	private static final Envelope EXTENT = new Envelope(0, 1000, 0, 1000);

	private GeometryFactory gf;
	private List<Coordinate> circles;
	private List<Geometry> circlePolygons;
	private int quadSegsUsed;

	@Param({ "100", "1000", "3000" })
	public int n;

	@Param({ "5.0", "20.0" })
	public double radius;

	@Param({ "1.0" })
	public double maxSegmentLength;

	@Param({ "1e-9" })
	public double eps;

	@Setup(Level.Trial)
	public void setup() {
		gf = new GeometryFactory();
		circles = new ArrayList<>(n);
		circlePolygons = new ArrayList<>(n);
		quadSegsUsed = quadSegsForMaxSegmentLength(radius, maxSegmentLength);

		double minX = EXTENT.getMinX();
		double maxX = EXTENT.getMaxX();
		double minY = EXTENT.getMinY();
		double maxY = EXTENT.getMaxY();

		SplittableRandom rnd = new SplittableRandom(SEED);
		for (int i = 0; i < n; i++) {
			double x = rnd.nextDouble(minX, maxX);
			double y = rnd.nextDouble(minY, maxY);
			circles.add(new Coordinate(x, y, radius));
			Point p = gf.createPoint(new Coordinate(x, y));
			Polygon poly = (Polygon) p.buffer(radius, quadSegsUsed);
			circlePolygons.add(poly);
		}

		List<Coordinate> warmupCircles = circles.subList(0, Math.min(circles.size(), 16));
		var warmupBoundary = DiskUnion.computeBoundaryArcs(warmupCircles, eps);
		DiskUnion.toGeometry(warmupBoundary, gf, maxSegmentLength);

		List<Geometry> warmupPolys = new ArrayList<>(circlePolygons.subList(0, Math.min(circlePolygons.size(), 16)));
		CascadedPolygonUnion.union(warmupPolys);
	}

	@Benchmark
	public void testDiskUnion(Blackhole bh) {
		DiskUnion.ArcBoundary boundary = DiskUnion.computeBoundaryArcs(circles, eps);
		Geometry result = DiskUnion.toGeometry(boundary, gf, maxSegmentLength);
		bh.consume(result);
	}

	@Benchmark
	public void testCascadedPolygonUnion(Blackhole bh) {
		List<Geometry> input = new ArrayList<>(circlePolygons);
		Geometry result = CascadedPolygonUnion.union(input);
		bh.consume(result);
	}

	private static int quadSegsForMaxSegmentLength(double radius, double maxSegLen) {
		if (radius <= 0 || maxSegLen <= 0) {
			return 16;
		}
		double ratio = maxSegLen / (2.0 * radius);
		if (ratio >= 1.0) {
			return 1;
		}
		double angle = Math.asin(Math.max(0.0, ratio));
		int quadSegs = (int) Math.ceil(Math.PI / (4.0 * angle));
		return Math.max(1, quadSegs);
	}
}
