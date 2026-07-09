package com.github.micycle1.geoblitz;

import org.locationtech.jts.algorithm.distance.DiscreteHausdorffDistance;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.util.GeometricShapeFactory;
import org.openjdk.jmh.annotations.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;

/**
 * Benchmark comparing ProHausdorffDistance vs JTS DiscreteHausdorffDistance.
 */
@State(Scope.Benchmark)
@BenchmarkMode(Mode.Throughput)
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 2, time = 1)
@Measurement(iterations = 3, time = 1)
@Fork(value = 1, jvmArgsAppend = { "-Xms1g", "-Xmx4g" })
public class BenchmarkProHD {

	@Param({ "100", "250", "1000" })
	public int nVertex;

	private List<Geometry> geomsA;
	private List<Geometry> geomsB;
	private static final int BATCH_SIZE = 100;

	@Setup(Level.Trial)
	public void setup() {
		geomsA = new ArrayList<>(BATCH_SIZE);
		geomsB = new ArrayList<>(BATCH_SIZE);
		GeometryFactory gf = new GeometryFactory();
		GeometricShapeFactory gsf = new GeometricShapeFactory(gf);
		Random rand = new Random(1337);

		for (int i = 0; i < BATCH_SIZE; i++) {
			// Generate Geometry A
			gsf.setNumPoints(nVertex);
			gsf.setCentre(new Coordinate(rand.nextDouble() * 100, rand.nextDouble() * 100));
			gsf.setSize(rand.nextDouble() * 20 + 5);
			geomsA.add(gsf.createCircle());

			// Generate Geometry B
			gsf.setNumPoints(nVertex);
			gsf.setCentre(new Coordinate(rand.nextDouble() * 100, rand.nextDouble() * 100));
			gsf.setSize(rand.nextDouble() * 20 + 5);
			geomsB.add(gsf.createRectangle());
		}
	}

	@Benchmark
	@OperationsPerInvocation(BATCH_SIZE)
	public double benchProHD() {
		final double alpha = 0.1; // NOTE sensible alpha value
		double total = 0;
		for (int i = 0; i < BATCH_SIZE; i++) {
			var result = ProHausdorffDistance.distance(geomsA.get(i), geomsB.get(i), alpha, 0);
			total += result;
		}
		return total;
	}

	@Benchmark
	@OperationsPerInvocation(BATCH_SIZE)
	public double benchJTS() {
		double total = 0;
		for (int i = 0; i < BATCH_SIZE; i++) {
			double dist = DiscreteHausdorffDistance.distance(geomsA.get(i), geomsB.get(i));
			total += dist;
		}
		return total;
	}

	public static void main(String[] args) {
		System.out.println("Running Accuracy Check...");
		checkAccuracySummary();
	}

	private static void checkAccuracySummary() {
		int[] sizes = { 50, 200, 1000 };
		double[] alphas = { 0.01, 0.05, 0.1, 0.2, 0.5 };
		Random rand = new Random(42);
		int samples = 15;

		// Warmup to prevent JIT compilation overhead from skewing the first results
		for (int i = 0; i < 100; i++) {
			var a = GeomMaker.make(50, i);
			var b = GeomMaker.make(50, i + 1);
			DiscreteHausdorffDistance.distance(a, b);
			ProHausdorffDistance.distance(a, b, 0.1, 0);
		}

		// test on very complex shape
		System.out.println("Average Error % and Computation Time (ProHD/JTS ms) by Alpha and Vertex Count");
		System.out.println("=".repeat(115));
		System.out.printf("%-10s", "NVertex");
		for (double alpha : alphas) {
			System.out.printf(" %-20s", "α=" + String.format("%.2f", alpha));
		}
		System.out.println();
		System.out.println("-".repeat(115));

		for (int n : sizes) {
			System.out.printf("%-10d", n);

			for (double alpha : alphas) {
				double totalErrorPct = 0;
				long totalTimeJTS = 0;
				long totalTimePro = 0;

				for (int i = 0; i < samples; i++) {
					var a = GeomMaker.make(n, rand.nextLong());
					var b = GeomMaker.make(n, rand.nextLong());

					long startJTS = System.nanoTime();
					double jtsDist = DiscreteHausdorffDistance.distance(a, b);
					long endJTS = System.nanoTime();
					totalTimeJTS += (endJTS - startJTS);

					long startPro = System.nanoTime();
					double proDist = ProHausdorffDistance.distance(a, b, alpha, 0);
					long endPro = System.nanoTime();
					totalTimePro += (endPro - startPro);

					double diff = jtsDist - proDist;
					double errPct = (Math.abs(diff) / jtsDist * 100);

					totalErrorPct += errPct;
				}
				double avgError = totalErrorPct / samples;
				double avgTimeProMs = (totalTimePro / 1_000_000.0) / samples;
				double avgTimeJTSMs = (totalTimeJTS / 1_000_000.0) / samples;
				System.out.printf(" %-20s", String.format("%.2f%% (%.3f/%.3f)", avgError, avgTimeProMs, avgTimeJTSMs));
			}
			System.out.println();
		}
	}
}
