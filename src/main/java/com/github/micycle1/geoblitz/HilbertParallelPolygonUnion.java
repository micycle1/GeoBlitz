package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.IntStream;

import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.geom.util.PolygonExtracter;
import org.locationtech.jts.index.hprtree.HilbertEncoder;
import org.locationtech.jts.shape.fractal.HilbertCode;

/**
 * High-performance polygon union built on Hilbert-ordered batching and parallel
 * reduction.
 * <p>
 * This utility provides a faster alternative to JTS's
 * {@code CascadedPolygonUnion} for large collections by:
 * <ul>
 * <li>sorting inputs in-place by a Hilbert space-filling-curve key derived from
 * their envelopes, improving spatial locality of subsequent unions; and</li>
 * <li>performing the union as a parallel reduction using Java streams.</li>
 * </ul>
 * The algorithm is intended for polygonal inputs (Polygon or MultiPolygon). If
 * intermediate union results contain non-polygonal artifacts (e.g., due to
 * precision effects), those are discarded and only polygonal components are
 * returned (as a Polygon if one component remains, otherwise as a
 * MultiPolygon).
 * <p>
 * Notes:
 * <ul>
 * <li>The input list is reordered in-place.</li>
 * <li>The list must be non-empty; the result {@code GeometryFactory} is taken
 * from the first element.</li>
 * </ul>
 *
 * @author Michael Carleton
 * @see org.locationtech.jts.operation.union.CascadedPolygonUnion
 *      CascadedPolygonUnion
 */
public class HilbertParallelPolygonUnion {

	private HilbertParallelPolygonUnion() {
	}

	/**
	 * Computes the union of the provided geometry using Hilbert-order sorting and a
	 * parallel reduction.
	 * <p>
	 * The input geometries are first sorted in-place by a Hilbert curve key
	 * computed from each geometry's envelope to cluster nearby geometries. The
	 * union is then evaluated in parallel.
	 * <p>
	 * Any non-polygonal components produced during union are discarded, and the
	 * returned geometry is guaranteed to be polygonal:
	 * <ul>
	 * <li>a {@code Polygon} if the result contains a single polygon, or</li>
	 * <li>a {@code MultiPolygon} if multiple polygons remain.</li>
	 * </ul>
	 * The result is built using the {@code GeometryFactory} of the first input
	 * geometry.
	 *
	 * @param geom a geometry (intended for Polygon collection or MultiPolygon)
	 * @return a polygonal geometry representing the union (Polygon or MultiPolygon)
	 */
	public static Geometry union(Geometry geom) {
		return union(PolygonExtracter.getPolygons(geom));
	}

	/**
	 * Computes the union of the provided geometries using Hilbert-order sorting and
	 * a parallel reduction.
	 * <p>
	 * The input list is first sorted in-place by a Hilbert curve key computed from
	 * each geometry's envelope to cluster nearby geometries. The union is then
	 * evaluated in parallel.
	 * <p>
	 * Any non-polygonal components produced during union are discarded, and the
	 * returned geometry is guaranteed to be polygonal:
	 * <ul>
	 * <li>a {@code Polygon} if the result contains a single polygon, or</li>
	 * <li>a {@code MultiPolygon} if multiple polygons remain.</li>
	 * </ul>
	 * The result is built using the {@code GeometryFactory} of the first input
	 * geometry.
	 *
	 * @param geoms a non-null, non-empty list of geometries (intended for Polygon
	 *              or MultiPolygon)
	 * @return a polygonal geometry representing the union (Polygon or MultiPolygon)
	 */
	public static Geometry union(List<Geometry> geoms) {
		if (geoms.isEmpty()) {
			return new GeometryFactory().createEmpty(2);
		}
		int n = geoms.size();
		geoms = new ArrayList<>(geoms); // copy for mutation
		sort(geoms, HilbertCode.level(n)); // sort according to center point of MBR
		var factory = geoms.get(0).getFactory();

		return geoms.parallelStream().reduce((g1, g2) -> {
			var result = g1.union(g2);
			if (!(result instanceof Polygonal)) {
				var polygons = PolygonExtracter.getPolygons(result);
				if (polygons.size() == 1) {
					result = (Polygon) polygons.get(0);
				} else {
					result = factory.buildGeometry(polygons);
				}
			}
			return result;
		}).orElse(factory.createEmpty(2));
	}

	/**
	 * Sorts a list of {@link Geometry} objects in-place by their spatial order
	 * using Hilbert curve encoding of their envelopes.
	 *
	 * @param geoms the list of geometries to sort
	 * @param level the resolution level for Hilbert curve encoding
	 */
	private static void sort(List<? extends Geometry> geoms, int level) {
		int n = geoms.size();
		if (n < 2) {
			return;
		}

		Envelope globalExtent = new Envelope();
		for (Geometry g : geoms) {
			globalExtent.expandToInclude(g.getEnvelopeInternal());
		}

		HilbertEncoder encoder = new HilbertEncoder(level, globalExtent);
		int[] keys = new int[n];
		for (int i = 0; i < n; i++) {
			Envelope e = geoms.get(i).getEnvelopeInternal();
			keys[i] = encoder.encode(e);
		}
		sortInPlaceByKeys(keys, geoms);
	}

	/**
	 * Used by sort().
	 */
	private static <T> void sortInPlaceByKeys(int[] keys, List<T> values) {
		final int n = keys.length;

		Integer[] idx = IntStream.range(0, n).boxed().toArray(Integer[]::new);
		Arrays.sort(idx, Comparator.comparingInt(i -> keys[i]));

		// rearrange keys and values in-place by following permutation cycles,
		// so that both arrays are sorted according to hilbert order key.
		boolean[] seen = new boolean[n];
		for (int i = 0; i < n; i++) {
			if (seen[i] || idx[i] == i) {
				continue;
			}

			int cycleStart = i;
			int j = i;
			int savedKey = keys[j];
			T savedVal = values.get(j);

			do {
				seen[j] = true;
				int next = idx[j];
				keys[j] = keys[next];
				values.set(j, values.get(next));

				j = next;
			} while (j != cycleStart);

			keys[j] = savedKey;
			values.set(j, savedVal);
			seen[j] = true;
		}
	}
}
