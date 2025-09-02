package com.github.micycle1.geoblitz;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.IntStream;

import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.geom.util.PolygonExtracter;
import org.locationtech.jts.index.hprtree.HilbertEncoder;
import org.locationtech.jts.shape.fractal.HilbertCode;

/**
 * Polygon union.
 * 
 * @author Michael Carleton
 */
public class HilbertParallelPolygonUnion {

	private HilbertParallelPolygonUnion() {
	}

	/**
	 * Much faster CascadedPolygonUnion (similar idea, but faster sorting and
	 * parallel union).
	 */
	public static Geometry union(List<Geometry> geoms) {
		int n = geoms.size();
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
		if (n < 2)
			return;

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
			if (seen[i] || idx[i] == i)
				continue;

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
