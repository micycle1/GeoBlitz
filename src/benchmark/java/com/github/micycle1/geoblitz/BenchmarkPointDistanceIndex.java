package com.github.micycle1.geoblitz;

import java.util.SplittableRandom;
import java.util.concurrent.TimeUnit;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.operation.distance.IndexedFacetDistance;
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
@Fork(value = 1, jvmArgsAppend = { "-Xms1g", "-Xmx2g" })
public class BenchmarkPointDistanceIndex {

    private Geometry geometry;
    private Coordinate[] queryPoints;
    private PointDistanceIndex ild;
    private IndexedFacetDistance ifd;

    private static final int NUM_QUERIES = 10_000;

    @Param({ "100", "1000", "5000" })
    public int nPoints;

    @Setup(Level.Trial)
    public void setup() {
        long seed = 42L;
        // Create a random geometry (likely a polygon)
        geometry = GeomMaker.make(nPoints, seed);

        // Initialize both indexes
        ild = new PointDistanceIndex(geometry);
        ifd = new IndexedFacetDistance(geometry);

        // Generate random query points
        queryPoints = new Coordinate[NUM_QUERIES];
        SplittableRandom rnd = new SplittableRandom(seed + 1);

        // Get envelope to sample points near the geometry
        var env = geometry.getEnvelopeInternal();
        double w = env.getWidth();
        double h = env.getHeight();
        double minX = env.getMinX() - w * 0.2;
        double maxX = env.getMaxX() + w * 0.2;
        double minY = env.getMinY() - h * 0.2;
        double maxY = env.getMaxY() + h * 0.2;

        for (int i = 0; i < NUM_QUERIES; i++) {
            double x = rnd.nextDouble(minX, maxX);
            double y = rnd.nextDouble(minY, maxY);
            queryPoints[i] = new Coordinate(x, y);
        }
    }

    @Benchmark
    @OperationsPerInvocation(NUM_QUERIES)
    public void benchPointDistanceIndex(Blackhole bh) {
        for (Coordinate c : queryPoints) {
            bh.consume(ild.unsignedDistance(c));
        }
    }

    @Benchmark
    @OperationsPerInvocation(NUM_QUERIES)
    public void benchIndexedFacetDistance(Blackhole bh) {
        for (Coordinate c : queryPoints) {
            Point p = geometry.getFactory().createPoint(c);
            bh.consume(ifd.distance(p));
        }
    }
}
