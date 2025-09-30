package com.github.micycle1.geoblitz;

import org.locationtech.jts.algorithm.CGAlgorithmsDD;
import org.locationtech.jts.algorithm.Distance;
import org.locationtech.jts.algorithm.LineIntersector;
import org.locationtech.jts.algorithm.RobustDeterminant;
import org.locationtech.jts.algorithm.RobustLineIntersector;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;

/**
 * Fast, pragmatic line segment intersection implementation optimized for speed.
 *
 * <p>
 * This implementation is similar to {@link RobustLineIntersector} but
 * sacrifices some of the extreme numerical robustness techniques used there
 * (for example, double-double arithmetic and determinant exactness) in favor of
 * performance.
 * </p>
 * <p>
 * Important notes and trade-offs:
 * </p>
 * <ul>
 * <li>Speed over ultimate robustness — this class is appropriate when
 * performance is the priority and occasional, extremely ill-conditioned corner
 * cases can be tolerated or handled at a higher level.</li>
 * <li>Uses a fast orientation test (with a small relative epsilon check) and
 * only falls back to {@link RobustDeterminant#signOfDet2x2} in borderline
 * cases.</li>
 * <li>No special Z interpolation is performed (the implementation returns
 * copies of input coordinates for collinear overlaps and endpoints).</li>
 * <li>A {@code precisionModel} (inherited from {@link LineIntersector}) is
 * applied to computed intersection points if present.</li>
 * <li>The instance is stateful (it records input lines, intersection points,
 * and flags such as {@code isProper}) and therefore is not thread-safe for
 * concurrent use. Create a separate instance per thread or call site if
 * concurrent access is required.</li>
 * </ul>
 * 
 * @author Michael Carleton
 */
public class FastLineIntersector extends LineIntersector {

	public FastLineIntersector() {
	}

	@Override
	public void computeIntersection(Coordinate p, Coordinate p1, Coordinate p2) {
		isProper = false;
		// do between check first, since it is faster than the orientation test
		if (Envelope.intersects(p1, p2, p)) {
			if ((orientSign(p1, p2, p) == 0) && (orientSign(p2, p1, p) == 0)) {
				isProper = true;
				if (p.equals(p1) || p.equals(p2)) {
					isProper = false;
				}
				result = POINT_INTERSECTION;
				return;
			}
		}
		result = NO_INTERSECTION;
	}

	@Override
	protected int computeIntersect(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
		isProper = false;

		// first try a fast test to see if the envelopes of the lines intersect
		if (!Envelope.intersects(p1, p2, q1, q2)) {
			return NO_INTERSECTION;
		}

		// for each endpoint, compute which side of the other segment it lies
		// if both endpoints lie on the same side of the other segment,
		// the segments do not intersect
		int Pq1 = orientSign(p1, p2, q1);
		int Pq2 = orientSign(p1, p2, q2);

		if ((Pq1 > 0 && Pq2 > 0) || (Pq1 < 0 && Pq2 < 0)) {
			return NO_INTERSECTION;
		}

		int Qp1 = orientSign(q1, q2, p1);
		int Qp2 = orientSign(q1, q2, p2);

		if ((Qp1 > 0 && Qp2 > 0) || (Qp1 < 0 && Qp2 < 0)) {
			return NO_INTERSECTION;
		}
		/**
		 * Intersection is collinear if each endpoint lies on the other line.
		 */
		boolean collinear = Pq1 == 0 && Pq2 == 0 && Qp1 == 0 && Qp2 == 0;
		if (collinear) {
			return computeCollinearIntersection(p1, p2, q1, q2);
		}

		/*
		 * At this point we know that there is a single intersection point (since the
		 * lines are not collinear).
		 */

		/**
		 * Check if the intersection is an endpoint. If it is, copy the endpoint as the
		 * intersection point. Copying the point rather than computing it ensures the
		 * point has the exact value, which is important for robustness. It is
		 * sufficient to simply check for an endpoint which is on the other line, since
		 * at this point we know that the inputLines must intersect.
		 */
		Coordinate p = null;
		if (Pq1 == 0 || Pq2 == 0 || Qp1 == 0 || Qp2 == 0) {
			isProper = false;

			if (p1.equals2D(q1)) {
				p = p1;
			} else if (p1.equals2D(q2)) {
				p = p1;
			} else if (p2.equals2D(q1)) {
				p = p2;
			} else if (p2.equals2D(q2)) {
				p = p2;
			}
			/**
			 * Now check to see if any endpoint lies on the interior of the other segment.
			 */
			else if (Pq1 == 0) {
				p = q1;
			} else if (Pq2 == 0) {
				p = q2;
			} else if (Qp1 == 0) {
				p = p1;
			} else if (Qp2 == 0) {
				p = p2;
			}
		} else {
		}
		isProper = true;
		p = intersection(p1, p2, q1, q2);
		intPt[0] = p.copy();
		return POINT_INTERSECTION;
	}

	private int computeCollinearIntersection(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
		boolean q1inP = Envelope.intersects(p1, p2, q1);
		boolean q2inP = Envelope.intersects(p1, p2, q2);
		boolean p1inQ = Envelope.intersects(q1, q2, p1);
		boolean p2inQ = Envelope.intersects(q1, q2, p2);

		if (q1inP && q2inP) {
			intPt[0] = q1.copy();
			intPt[1] = q2.copy();
			return COLLINEAR_INTERSECTION;
		}
		if (p1inQ && p2inQ) {
			intPt[0] = p1.copy();
			intPt[1] = p2.copy();
			return COLLINEAR_INTERSECTION;
		}
		if (q1inP && p1inQ) {
			intPt[0] = q1.copy();
			intPt[1] = p1.copy();
			return q1.equals(p1) && !q2inP && !p2inQ ? POINT_INTERSECTION : COLLINEAR_INTERSECTION;
		}
		if (q1inP && p2inQ) {
			intPt[0] = q1.copy();
			intPt[1] = p2.copy();
			return q1.equals(p2) && !q2inP && !p1inQ ? POINT_INTERSECTION : COLLINEAR_INTERSECTION;
		}
		if (q2inP && p1inQ) {
			intPt[0] = q2.copy();
			intPt[1] = p1.copy();
			return q2.equals(p1) && !q1inP && !p2inQ ? POINT_INTERSECTION : COLLINEAR_INTERSECTION;
		}
		if (q2inP && p2inQ) {
			intPt[0] = q2.copy();
			intPt[1] = p2.copy();
			return q2.equals(p2) && !q1inP && !p1inQ ? POINT_INTERSECTION : COLLINEAR_INTERSECTION;
		}
		return NO_INTERSECTION;
	}

	/**
	 * This method computes the actual value of the intersection point. It is
	 * rounded to the precision model if being used.
	 */
	private Coordinate intersection(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
		Coordinate intPt = intersectionSafe(p1, p2, q1, q2);

		if (!isInSegmentEnvelopes(intPt)) {
			intPt = nearestEndpoint(p1, p2, q1, q2).copy();
		}
		if (precisionModel != null) {
			precisionModel.makePrecise(intPt);
		}
		return intPt;
	}

	/**
	 * Computes a segment intersection. Round-off error can cause the raw
	 * computation to fail, (usually due to the segments being approximately
	 * parallel). If this happens, a reasonable approximation is computed instead.
	 *
	 * @param p1 a segment endpoint
	 * @param p2 a segment endpoint
	 * @param q1 a segment endpoint
	 * @param q2 a segment endpoint
	 * @return the computed intersection point
	 */
	private Coordinate intersectionSafe(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
		Coordinate intPt = intersectionFast(p1, p2, q1, q2);
		if (intPt == null) {
			intPt = nearestEndpoint(p1, p2, q1, q2);
		}
		return intPt;
	}

	/**
	 * Tests whether a point lies in the envelopes of both input segments. A
	 * correctly computed intersection point should return <code>true</code> for
	 * this test. Since this test is for debugging purposes only, no attempt is made
	 * to optimize the envelope test.
	 *
	 * @return <code>true</code> if the input point lies within both input segment
	 *         envelopes
	 */
	private boolean isInSegmentEnvelopes(Coordinate intPt) {
		Envelope env0 = new Envelope(inputLines[0][0], inputLines[0][1]);
		Envelope env1 = new Envelope(inputLines[1][0], inputLines[1][1]);
		return env0.contains(intPt) && env1.contains(intPt);
	}

	/**
	 * Finds the endpoint of the segments P and Q which is closest to the other
	 * segment. This is a reasonable surrogate for the true intersection points in
	 * ill-conditioned cases (e.g. where two segments are nearly coincident, or
	 * where the endpoint of one segment lies almost on the other segment).
	 * <p>
	 * This replaces the older CentralEndpoint heuristic, which chose the wrong
	 * endpoint in some cases where the segments had very distinct slopes and one
	 * endpoint lay almost on the other segment.
	 *
	 * @param p1 an endpoint of segment P
	 * @param p2 an endpoint of segment P
	 * @param q1 an endpoint of segment Q
	 * @param q2 an endpoint of segment Q
	 * @return the nearest endpoint to the other segment
	 */
	private static Coordinate nearestEndpoint(Coordinate p1, Coordinate p2, Coordinate q1, Coordinate q2) {
		Coordinate nearestPt = p1;
		double minDist = Distance.pointToSegmentSq(p1, q1, q2);

		double dist = Distance.pointToSegmentSq(p2, q1, q2);
		if (dist < minDist) {
			minDist = dist;
			nearestPt = p2;
		}
		dist = Distance.pointToSegmentSq(q1, p1, p2);
		if (dist < minDist) {
			minDist = dist;
			nearestPt = q1;
		}
		dist = Distance.pointToSegmentSq(q2, p1, p2);
		if (dist < minDist) {
			minDist = dist;
			nearestPt = q2;
		}
		return nearestPt;
	}

	private static Coordinate intersectionFast(final Coordinate p1, final Coordinate p2, final Coordinate q1, final Coordinate q2) {
		final double rx = p2.x - p1.x;
		final double ry = p2.y - p1.y;
		final double sx = q2.x - q1.x;
		final double sy = q2.y - q1.y;

		final double den = Math.fma(rx, sy, -ry * sx); // cross(r, s) != 0 by assumption
		final double dx = q1.x - p1.x;
		final double dy = q1.y - p1.y;
		final double t = Math.fma(dx, sy, -dy * sx) / den; // cross(q1-p1, s) / den

		return new Coordinate(Math.fma(t, rx, p1.x), Math.fma(t, ry, p1.y));
	}

	private static int orientSign(final Coordinate p1, final Coordinate p2, final Coordinate q) {
		final double dx1 = p2.x - p1.x;
		final double dy1 = p2.y - p1.y;
		final double dx2 = q.x - p2.x;
		final double dy2 = q.y - p2.y;

		final double det = Math.fma(dx1, dy2, -dy1 * dx2);

		if (det > 0.0) {
			return 1;
		}
		if (det < 0.0) {
			return -1;
		}

		// If det == 0.0 exactly, or underflow/rounding made it exactly zero,
		// do a more careful check using a tiny relative epsilon based on input scale:
		double maxAbs = Math.abs(dx1 * dy2) + Math.abs(dy1 * dx2);
		// relative tolerance; tune if needed
		double eps = 1e-12 * maxAbs;
		if (det > eps) {
			return 1;
		}
		if (det < -eps) {
			return -1;
		}

		// Otherwise fall back to the robust predicate (expensive)
		return CGAlgorithmsDD.signOfDet2x2(dx1, dy1, dx2, dy2);
	}

}
