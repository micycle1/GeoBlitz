package com.github.micycle1.geoblitz;

import java.util.SplittableRandom;
import java.util.concurrent.TimeUnit;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.algorithm.RobustLineIntersector;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OperationsPerInvocation;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

/**
 * Benchmark comparing FastLineIntersector vs RobustLineIntersector.
 *
 * Measures throughput of computeIntersect(p1,p2,q1,q2) over a fixed batch of
 * randomly generated segment pairs.
 */
@State(Scope.Benchmark)
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 0, time = 1)
@Measurement(iterations = 1, time = 2)
@Fork(value = 1, jvmArgsAppend = { "-Xms1g", "-Xmx1g", "-XX:+AlwaysPreTouch" })
public class BenchmarkLineIntersector {

	// number of queries performed per benchmark invocation (kept constant)
	private static final int numQueries = 5000;

	// random seed for reproducibility
	private static final long SEED = 42L;

	// Request different problem sizes (not strictly necessary, but convenient)
	@Param({ "1000", "10000" })
	public int n; // unused as vertex count; kept for parity with your example

	// Two implementations to compare
	private FastLineIntersector fast;
	private RobustLineIntersector robust;

	// Pre-generated segment endpoints to use as queries
	private Coordinate[] p1s;
	private Coordinate[] p2s;
	private Coordinate[] q1s;
	private Coordinate[] q2s;

	@org.openjdk.jmh.annotations.Setup(org.openjdk.jmh.annotations.Level.Trial)
	public void setup() {
		SplittableRandom rnd = new SplittableRandom(SEED);

		// allocate arrays
		p1s = new Coordinate[numQueries];
		p2s = new Coordinate[numQueries];
		q1s = new Coordinate[numQueries];
		q2s = new Coordinate[numQueries];

		// Generate random segment pairs. Use a distribution that produces
		// a mixture of crossing, near-parallel, and non-intersecting cases.
		for (int i = 0; i < numQueries; i++) {
			// center points
			double cx = (rnd.nextDouble() - 0.5) * 1000.0;
			double cy = (rnd.nextDouble() - 0.5) * 1000.0;

			// random small direction and length for first segment
			double angle1 = rnd.nextDouble() * Math.PI * 2.0;
			double len1 = 0.1 + rnd.nextDouble() * 5.0;
			double dx1 = Math.cos(angle1) * len1;
			double dy1 = Math.sin(angle1) * len1;
			p1s[i] = new Coordinate(cx - dx1 * 0.5, cy - dy1 * 0.5);
			p2s[i] = new Coordinate(cx + dx1 * 0.5, cy + dy1 * 0.5);

			// random small direction and length for second segment (slightly offset)
			double angle2 = angle1 + (rnd.nextDouble() - 0.5) * 0.5; // often near-parallel
			double len2 = 0.1 + rnd.nextDouble() * 5.0;
			double dx2 = Math.cos(angle2) * len2;
			double dy2 = Math.sin(angle2) * len2;

			// offset the second segment center randomly so we get variety
			double ox = cx + (rnd.nextDouble() - 0.5) * 2.0;
			double oy = cy + (rnd.nextDouble() - 0.5) * 2.0;
			q1s[i] = new Coordinate(ox - dx2 * 0.5, oy - dy2 * 0.5);
			q2s[i] = new Coordinate(ox + dx2 * 0.5, oy + dy2 * 0.5);

			// occasionally produce an exact endpoint-touch to test endpoint logic
			if (i % 50 == 0) {
				if (rnd.nextBoolean()) {
					q1s[i] = p1s[i].copy();
				} else {
					q2s[i] = p2s[i].copy();
				}
			}
		}

		// create instances
		fast = new FastLineIntersector();
		robust = new RobustLineIntersector();

		// touch both once to avoid classloading/jit overhead during measured iterations
		fast.computeIntersection(p1s[0], p2s[0], q1s[0], q2s[0]);
		robust.computeIntersection(p1s[0], p2s[0], q1s[0], q2s[0]);
		// consume to avoid dead code elimination if JIT tries too hard
	}

	// ----------------------------
	// computeIntersect(...) benchmarks
	// ----------------------------

	@Benchmark
	@OperationsPerInvocation(numQueries)
	public void b0_computeIntersect_Fast(Blackhole bh) {
		int res = 0;
		for (int i = 0; i < numQueries; i++) {
			fast.computeIntersection(p1s[i], p2s[i], q1s[i], q2s[i]);
		}
		bh.consume(res);
	}

	@Benchmark
	@OperationsPerInvocation(numQueries)
	public void b1_computeIntersect_Robust(Blackhole bh) {
		int res = 0;
		for (int i = 0; i < numQueries; i++) {
			robust.computeIntersection(p1s[i], p2s[i], q1s[i], q2s[i]);
		}
		bh.consume(res);
	}
}