package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
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
public class BenchmarkHilbertParallelPolygonUnion {

	private static final long seed = 42L;
	private static final Envelope extent = new Envelope(0, 1000, 0, 1000);

	private GeometryFactory gf;
	private List<Polygon> polygons;

	// number of polygon inputs
	@Param({ "100", "1000",  })
	public int n;

	// radius of each buffered point (controls overlap with given extent)
	@Param({ "10.0" })
	public double radius;

//	@Param({ "8", "16" })
	public int quadSegs = 16;

	@Setup(Level.Trial)
	public void setup() {
		gf = new GeometryFactory();
		polygons = GeomMaker.makeRandomBufferedPoints(n, radius, quadSegs, seed, extent);

		// Trigger class loading/JIT outside measurement
		List<Geometry> warmupSmall = new ArrayList<>(polygons.subList(0, Math.min(polygons.size(), 16)));
		CascadedPolygonUnion.union(warmupSmall);
		HilbertParallelPolygonUnion.union(new ArrayList<>(warmupSmall));
	}

	// TODO node all lines, polygonise, coverage union (works on contiguous blobs only?)

	@Benchmark
	public void testCascadedPolygonUnion(Blackhole bh) {
		// copy the list to avoid any accidental mutation from algorithms
		List<Geometry> input = new ArrayList<>(polygons);
		Geometry result = CascadedPolygonUnion.union(input);
		bh.consume(result);
	}

	@Benchmark
	public void testHilbertParallelPolygonUnion(Blackhole bh) {
		// OmniUnion sorts in-place; give it a fresh list each time
		List<Geometry> input = new ArrayList<>(polygons);
		Geometry result = HilbertParallelPolygonUnion.union(input);
		bh.consume(result);
	}

//	@Benchmark
//	public void testNary(Blackhole bh) {
//		// OmniUnion sorts in-place; give it a fresh list each time
//		List<Polygon> input = new ArrayList<>(polygons);
//		Geometry result = NAryUnion.union(input, null);
//		bh.consume(result);
//	}
}