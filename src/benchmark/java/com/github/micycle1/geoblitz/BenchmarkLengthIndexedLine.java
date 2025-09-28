package com.github.micycle1.geoblitz;

import java.util.SplittableRandom;
import java.util.concurrent.TimeUnit;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.linearref.LengthIndexedLine;
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
 * Benchmark comparing LengthIndexedLine (baseline) versus
 * IndexedLengthIndexedLine (indexed) for extractPoint, extractPoint+offset, and extractLine.
 *
 * - Constructs a long LineString with 'n' vertices (n-1 segments).
 * - Generates numQueries random length indices (some negative, some out-of-range).
 * - Benchmarks throughput of performing each operation numQueries times.
 *
 * Usage (example):
 *   mvn -pl jmh:benchmark -Dbenchmark=org.locationtech.jts.linearref.BenchmarkLengthIndexedLine ...
 */
@State(Scope.Benchmark)
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 0, time = 1)
@Measurement(iterations = 1, time = 2)
@Fork(value = 1, jvmArgsAppend = { "-Xms1g", "-Xmx1g", "-XX:+AlwaysPreTouch" })
public class BenchmarkLengthIndexedLine {

    // number of vertices in the polyline (thus segments = n-1)
    @Param({ "1000", "10000" })
    public int n;

    // number of queries performed per benchmark invocation (kept constant)
    private static final int numQueries = 5000;

    // random seed for reproducibility
    private static final long SEED = 42L;

    private final GeometryFactory gf = new GeometryFactory();

    // The line geometry (single LineString)
    private Geometry line;

    // two implementations
    private LengthIndexedLine base;                 // original
    private IndexedLengthIndexedLine indexed;       // indexed variant

    // precomputed query indexes (length values), may be negative or out-of-range
    private double[] queryIndices;
    // precomputed query pairs for extractLine
    private double[] queryIndices2;

    // precomputed offsets for offseted extractPoint
    private double[] offsets;

    @org.openjdk.jmh.annotations.Setup(org.openjdk.jmh.annotations.Level.Trial)
    public void setup() {
        SplittableRandom rnd = new SplittableRandom(SEED);

        // Build a long random-ish polyline with n vertices
        Coordinate[] coords = new Coordinate[n];
        double x = 0.0, y = 0.0;
        final double step = 1.0; // small steps to produce many segments with varied lengths
        for (int i = 0; i < n; i++) {
            x += (rnd.nextDouble() - 0.5) * step;
            y += (rnd.nextDouble() - 0.5) * step;
            coords[i] = new Coordinate(x, y);
        }
        LineString ls = gf.createLineString(coords);
        this.line = ls;

        // create both implementations
        this.base = new LengthIndexedLine(line);
        this.indexed = new IndexedLengthIndexedLine(line);

        final double totalLength = line.getLength();

        // generate query indices. Some intentionally negative and some > totalLength
        queryIndices = new double[numQueries];
        queryIndices2 = new double[numQueries];
        offsets = new double[numQueries];
        for (int i = 0; i < numQueries; i++) {
            // sample index in [-0.2*T, 1.2*T]
            double t = (rnd.nextDouble() * 1.4 - 0.2) * totalLength;
            double u = (rnd.nextDouble() * 1.4 - 0.2) * totalLength;
            queryIndices[i] = t;
            queryIndices2[i] = u;
            // small offsets in [-2.0, 2.0]
            offsets[i] = (rnd.nextDouble() - 0.5) * 4.0;
        }

        // Touch both implementations once to avoid lazy-class-loading / JIT warmups interfering
        if (numQueries > 0) {
            base.extractPoint(queryIndices[0]);
            base.extractPoint(queryIndices[0], offsets[0]);
            base.extractLine(queryIndices[0], queryIndices2[0]);

            indexed.extractPoint(queryIndices[0]);
            indexed.extractPoint(queryIndices[0], offsets[0]);
            indexed.extractLine(queryIndices[0], queryIndices2[0]);
        }
    }

    // ----------------------------
    // extractPoint(index) benchmarks
    // ----------------------------

    @Benchmark
    @OperationsPerInvocation(numQueries)
    public void b0_extractPoint_Base(Blackhole bh) {
        for (int i = 0; i < numQueries; i++) {
            bh.consume(base.extractPoint(queryIndices[i]));
        }
    }

    @Benchmark
    @OperationsPerInvocation(numQueries)
    public void b1_extractPoint_Indexed(Blackhole bh) {
        for (int i = 0; i < numQueries; i++) {
            bh.consume(indexed.extractPoint(queryIndices[i]));
        }
    }

    // ----------------------------
    // extractPoint(index, offset) benchmarks
    // ----------------------------

    @Benchmark
    @OperationsPerInvocation(numQueries)
    public void b2_extractPointOffset_Base(Blackhole bh) {
        for (int i = 0; i < numQueries; i++) {
            bh.consume(base.extractPoint(queryIndices[i], offsets[i]));
        }
    }

    @Benchmark
    @OperationsPerInvocation(numQueries)
    public void b3_extractPointOffset_Indexed(Blackhole bh) {
        for (int i = 0; i < numQueries; i++) {
            bh.consume(indexed.extractPoint(queryIndices[i], offsets[i]));
        }
    }

    // ----------------------------
    // extractLine(start, end) benchmarks
    // ----------------------------

    @Benchmark
    @OperationsPerInvocation(numQueries)
    public void b4_extractLine_Base(Blackhole bh) {
        for (int i = 0; i < numQueries; i++) {
            bh.consume(base.extractLine(queryIndices[i], queryIndices2[i]));
        }
    }

    @Benchmark
    @OperationsPerInvocation(numQueries)
    public void b5_extractLine_Indexed(Blackhole bh) {
        for (int i = 0; i < numQueries; i++) {
            bh.consume(indexed.extractLine(queryIndices[i], queryIndices2[i]));
        }
    }
}