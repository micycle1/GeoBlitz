package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.List;
import java.util.SplittableRandom;

import org.locationtech.jts.algorithm.hull.ConcaveHull;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;

class GeomMaker {

	private static final GeometryFactory FACTORY = new GeometryFactory();

	static Geometry make(int points, long seed) {
		var coords = plasticJitteredLDS(0, 0, 1000, 1000, points, seed);
		var g = FACTORY.createMultiPointFromCoords(coords);
		ConcaveHull h = new ConcaveHull(g);
		var hull = h.getHull();
		return hull.buffer(0);
	}

	private static Coordinate[] plasticJitteredLDS(double xMin, double yMin, double xMax, double yMax, int n, long seed) {
	    final double w = xMax - xMin;
	    final double h = yMax - yMin;

	    final SplittableRandom random = new SplittableRandom(seed);
	    final double p = 1.32471795724474602596; // plastic constant
	    final double a1 = 1.0 / p;
	    final double a2 = 1.0 / (p * p);
	    final double c_magicNumber = 0.732;

	    final Coordinate[] points = new Coordinate[n];
	    for (int i = 0; i < n; i++) {
	        double x = (((random.nextDouble() * c_magicNumber / Math.sqrt(i + 1d) + a1 * i) % 1) * w + xMin);
	        double y = (((random.nextDouble() * c_magicNumber / Math.sqrt(i + 1d) + a2 * i) % 1) * h + yMin);
	        points[i] = new Coordinate(x, y);
	    }
	    return points;
	}

	static List<Polygon> makeRandomBufferedPoints(int count, double r, int qSegs, long seed, Envelope env) {
		SplittableRandom rnd = new SplittableRandom(seed);
		GeometryFactory gf = new GeometryFactory();

		double minX = env.getMinX();
		double maxX = env.getMaxX();
		double minY = env.getMinY();
		double maxY = env.getMaxY();

		List<Polygon> out = new ArrayList<>(count);
		for (int i = 0; i < count; i++) {
			double x = rnd.nextDouble(minX, maxX);
			double y = rnd.nextDouble(minY, maxY);
			Point p = gf.createPoint(new Coordinate(x, y));
			// buffer yields a valid polygonal disk; ensures polygonal input
			Polygon poly = (Polygon) p.buffer(r, qSegs);
			out.add(poly);
		}
		return out;
	}

}
