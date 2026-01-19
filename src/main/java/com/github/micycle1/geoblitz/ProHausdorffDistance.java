package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Objects;
import java.util.concurrent.ThreadLocalRandom;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;

/**
 * ProHD: Projection-based Hausdorff Distance approximation for JTS
 * {@link Geometry}.
 *
 * <p>
 * Intuition: the Hausdorff distance is typically determined by a small set of
 * "extreme" points. ProHD finds candidate extreme points by projecting each
 * geometry's sampled points onto a few informative directions (the
 * centroid-to-centroid axis and principal-component axis/axes), keeping only
 * the top/bottom {@code alpha} tails along each projection. It then computes
 * the Hausdorff distance on the union of these selected points, yielding a fast
 * underestimate of the Hausdorff distance of the sampled point sets.
 * </p>
 *
 * <p>
 * Reference: Jiuzhou Fu, Luanzheng Guo, Nathan R. Tallent, Dongfang Zhao,
 * <i>ProHD: Projection-Based Hausdorff Distance Approximation</i>,
 * arXiv:2511.18207.
 * </p>
 */
public final class ProHausdorffDistance {

	private ProHausdorffDistance() {
	}

	/**
	 * Approximates the (undirected) Hausdorff distance using ProHD: keep only
	 * points that are extreme under a few projection directions (centroid axis +
	 * PCA axis), then compute Hausdorff on those reduced point sets.
	 *
	 * @param a                non-empty geometry A
	 * @param b                non-empty geometry B
	 * @param alpha            tail fraction kept per projection direction, in (0,
	 *                         0.5]. For each direction and each set of n sampled
	 *                         points, we keep k = max(1, floor(alpha * n)) points
	 *                         with the smallest projections AND k points with the
	 *                         largest projections (≈ 2 * alpha * n per direction,
	 *                         before unions across directions). Typical values:
	 *                         0.005–0.05 (often 0.01).
	 * @param maxSegmentLength if > 0, densifies geometry segments so consecutive
	 *                         sample points are at most this far apart (in geometry
	 *                         units). If <= 0, uses only the geometry's existing
	 *                         vertices.
	 */
	public static double distance(Geometry a, Geometry b, double alpha, double maxSegmentLength) {
		Objects.requireNonNull(a, "a");
		Objects.requireNonNull(b, "b");
		if (a.isEmpty() || b.isEmpty()) {
			throw new IllegalArgumentException("Input geometries must be non-empty.");
		}
		if (!(alpha > 0.0 && alpha <= 0.5)) {
			throw new IllegalArgumentException("alpha must be in (0, 0.5].");
		}

		Coordinate[] ptsA = sampleGeometry(a, maxSegmentLength);
		Coordinate[] ptsB = sampleGeometry(b, maxSegmentLength);

		if (ptsA.length == 0 || ptsB.length == 0) {
			throw new IllegalArgumentException("Sampling produced no points (unexpected for non-empty geometry).");
		}

		// D = 2 for JTS, so m = floor(sqrt(D)) = 1
		final int D = 2;
		final int m = (int) Math.floor(Math.sqrt(D)); // = 1
		final double alphaPca = alpha / Math.max(1, m);

		BitSet selA = new BitSet(ptsA.length);
		BitSet selB = new BitSet(ptsB.length);

		// 1) Centroid direction selection
		selectByCentroidDirection(ptsA, ptsB, alpha, selA, selB);

		// 2) PCA direction(s) selection (2D => only one)
		selectByPcaDirections2D(ptsA, ptsB, alphaPca, m, selA, selB);

		Coordinate[] aSel = extractSelected(ptsA, selA);
		Coordinate[] bSel = extractSelected(ptsB, selB);

		// Compute Hausdorff on selected sets using nearest-neighbor queries via STRtree
		double hAB = directedHausdorff(aSel, bSel);
		double hBA = directedHausdorff(bSel, aSel);
		double H = Math.max(hAB, hBA);

		return H;
	}

	private static Coordinate[] sampleGeometry(Geometry g, double maxSegmentLength) {
		// Start from raw vertex coordinates
		Coordinate[] raw = g.getCoordinates();
		if (raw == null || raw.length == 0) {
			return new Coordinate[0];
		}

		if (!(maxSegmentLength > 0)) {
			// Use Arrays.copyOf instead of manual loop
			Coordinate[] out = new Coordinate[raw.length];
			System.arraycopy(raw, 0, out, 0, raw.length);
			return out;
		}

		ArrayList<Coordinate> pts = new ArrayList<>(raw.length);
		for (int i = 0; i < raw.length; i++) {
			Coordinate c0 = raw[i];
			pts.add(new Coordinate(c0));
			if (i == raw.length - 1) {
				break;
			}

			Coordinate c1 = raw[i + 1];
			double dx = c1.x - c0.x;
			double dy = c1.y - c0.y;
			double segLenSq = dx * dx + dy * dy;
			if (segLenSq <= 0) {
				continue;
			}

			double segLen = Math.sqrt(segLenSq);
			int steps = (int) Math.ceil(segLen / maxSegmentLength);
			// Insert intermediate points (exclude endpoint; it will be added as next vertex)
			if (steps > 1) {
				double stepInv = 1.0 / steps;
				for (int s = 1; s < steps; s++) {
					double t = s * stepInv;
					pts.add(new Coordinate(c0.x + t * dx, c0.y + t * dy));
				}
			}
		}
		return pts.toArray(new Coordinate[0]);
	}

	private static void selectByCentroidDirection(Coordinate[] A, Coordinate[] B, double alpha, BitSet selA, BitSet selB) {
		double[] cA = meanXY(A);
		double[] cB = meanXY(B);

		double ux = cB[0] - cA[0];
		double uy = cB[1] - cA[1];
		double nSq = ux * ux + uy * uy;
		if (nSq < 1e-24) { // degenerate: same centroid
			ux = 1.0;
			uy = 0.0;
		} else {
			double invN = 1.0 / Math.sqrt(nSq);
			ux *= invN;
			uy *= invN;
		}

		selectExtremesAlongDirection(A, ux, uy, alpha, selA);
		selectExtremesAlongDirection(B, ux, uy, alpha, selB);
	}

	/**
	 * For JTS 2D: compute top principal component of union(A,B) and select extremes
	 * along it. m is kept for API similarity with paper; for D=2, m will be 1.
	 */
	private static void selectByPcaDirections2D(Coordinate[] A, Coordinate[] B, double alpha, int m, BitSet selA, BitSet selB) {
		if (m <= 0) {
			return;
		}

		// Compute covariance of union
		int n = A.length + B.length;
		double meanX = 0, meanY = 0;
		for (Coordinate p : A) {
			meanX += p.x;
			meanY += p.y;
		}
		for (Coordinate p : B) {
			meanX += p.x;
			meanY += p.y;
		}
		double invN = 1.0 / n;
		meanX *= invN;
		meanY *= invN;

		double sxx = 0, sxy = 0, syy = 0;
		for (Coordinate p : A) {
			double dx = p.x - meanX, dy = p.y - meanY;
			sxx += dx * dx;
			sxy += dx * dy;
			syy += dy * dy;
		}
		for (Coordinate p : B) {
			double dx = p.x - meanX, dy = p.y - meanY;
			sxx += dx * dx;
			sxy += dx * dy;
			syy += dy * dy;
		}

		// Principal axis angle: 0.5 * atan2(2*sxy, sxx - syy)
		double angle = 0.5 * Math.atan2(2.0 * sxy, (sxx - syy));
		double ux = Math.cos(angle);
		double uy = Math.sin(angle);

		// Already unit length
		selectExtremesAlongDirection(A, ux, uy, alpha, selA);
		selectExtremesAlongDirection(B, ux, uy, alpha, selB);
	}

	private static void selectExtremesAlongDirection(Coordinate[] pts, double ux, double uy, double alpha, BitSet outSel) {
		int n = pts.length;
		int k = Math.max(1, (int) Math.floor(alpha * n));

		double[] proj = new double[n];
		for (int i = 0; i < n; i++) {
			Coordinate p = pts[i];
			proj[i] = p.x * ux + p.y * uy;
		}

		// Find thresholds for bottom-k and top-k using Quickselect
		double[] tmp = proj.clone();
		double low = quickSelect(tmp, k - 1); // kth smallest (0-based)
		double high = quickSelect(tmp, n - k); // (n-k)th smallest => kth largest threshold

		for (int i = 0; i < n; i++) {
			double v = proj[i];
			if (v <= low || v >= high) {
				outSel.set(i);
			}
		}
	}

	private static Coordinate[] extractSelected(Coordinate[] pts, BitSet sel) {
		Coordinate[] out = new Coordinate[sel.cardinality()];
		int idx = 0;
		for (int i = sel.nextSetBit(0); i >= 0; i = sel.nextSetBit(i + 1)) {
			out[idx++] = pts[i];
		}
		return out;
	}

	private static double[] meanXY(Coordinate[] pts) {
		double sx = 0, sy = 0;
		for (Coordinate p : pts) {
			sx += p.x;
			sy += p.y;
		}
		double invN = 1.0 / pts.length;
		return new double[] { sx * invN, sy * invN };
	}

	private static double directedHausdorff(Coordinate[] from, Coordinate[] to) {
		if (from.length == 0 || to.length == 0) {
			return 0.0;
		}

		HPRtreeX<Coordinate> tree2 = new HPRtreeX<>();
		for (Coordinate p : to) {
			tree2.insert(new Envelope(p), p);
		}

		double maxMinSq = 0.0;
		for (Coordinate q : from) {
			Coordinate nn = tree2.nearestNeighbor(q, (a, b) -> a.distance(b));
			double dSq = q.distanceSq(nn);
			if (dSq > maxMinSq) {
				maxMinSq = dSq;
			}
		}
		return Math.sqrt(maxMinSq);
	}

	private static double quickSelect(double[] a, int k) {
		if (k < 0 || k >= a.length) {
			throw new IllegalArgumentException("k out of range");
		}

		int left = 0, right = a.length - 1;
		ThreadLocalRandom rnd = ThreadLocalRandom.current();

		while (true) {
			if (left == right) {
				return a[left];
			}

			int pivotIndex = left + rnd.nextInt(right - left + 1);
			pivotIndex = partition(a, left, right, pivotIndex);

			if (k == pivotIndex) {
				return a[k];
			} else if (k < pivotIndex) {
				right = pivotIndex - 1;
			} else {
				left = pivotIndex + 1;
			}
		}
	}

	private static int partition(double[] a, int left, int right, int pivotIndex) {
		double pivotValue = a[pivotIndex];
		swap(a, pivotIndex, right);
		int storeIndex = left;

		for (int i = left; i < right; i++) {
			if (a[i] < pivotValue) {
				swap(a, storeIndex, i);
				storeIndex++;
			}
		}
		swap(a, right, storeIndex);
		return storeIndex;
	}

	private static void swap(double[] a, int i, int j) {
		double t = a[i];
		a[i] = a[j];
		a[j] = t;
	}
}