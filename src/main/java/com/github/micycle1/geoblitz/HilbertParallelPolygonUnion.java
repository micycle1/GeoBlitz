package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.RecursiveTask;

import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.geom.util.PolygonExtracter;
import org.locationtech.jts.index.hprtree.HilbertEncoder;
import org.locationtech.jts.operation.overlayng.OverlayNG;
import org.locationtech.jts.operation.overlayng.OverlayNGRobust;
import org.locationtech.jts.shape.fractal.HilbertCode;

/**
 * High-performance polygon union built on Hilbert-ordered batching and parallel
 * reduction.
 * <p>
 * This utility provides a faster alternative to JTS's
 * {@code CascadedPolygonUnion} for large collections by:
 * <ul>
 * <li>sorting inputs in-place by a Hilbert space-filling-curve key derived from
 * their envelopes, improving spatial locality of subsequent unions;</li>
 * <li>performing the union as a balanced binary-tree reduction on the fork/join
 * common pool, so intermediate operands stay comparable in size;</li>
 * <li>skipping overlay entirely for envelope-disjoint operands, and restricting
 * overlay to the envelope-intersecting elements otherwise (disjoint elements
 * are carried through unchanged).</li>
 * </ul>
 * The algorithm is intended for polygonal inputs (Polygon or MultiPolygon). If
 * intermediate union results contain non-polygonal artifacts (e.g., due to
 * precision effects), those are discarded and only polygonal components are
 * returned (as a Polygon if one component remains, otherwise as a
 * MultiPolygon).
 * <p>
 * Notes:
 * <ul>
 * <li>The input list is not mutated.</li>
 * <li>The result {@code GeometryFactory} is taken from the first element.</li>
 * </ul>
 *
 * @author Michael Carleton
 * @see org.locationtech.jts.operation.union.CascadedPolygonUnion
 *      CascadedPolygonUnion
 */
public class HilbertParallelPolygonUnion {

	/** Below this many inputs a task recurses sequentially instead of forking. */
	private static final int FORK_THRESHOLD = 4;

	private HilbertParallelPolygonUnion() {
	}

	/**
	 * Computes the union of the provided geometry using Hilbert-order sorting and a
	 * parallel reduction.
	 * <p>
	 * The input geometries are first sorted by a Hilbert curve key computed from
	 * each geometry's envelope to cluster nearby geometries. The union is then
	 * evaluated as a balanced binary-tree reduction in parallel.
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
		@SuppressWarnings("unchecked")
		List<Geometry> polygons = PolygonExtracter.getPolygons(geom);
		return union(polygons);
	}

	/**
	 * Computes the union of the provided geometries using Hilbert-order sorting and
	 * a parallel reduction.
	 * <p>
	 * The input list is first sorted by a Hilbert curve key computed from each
	 * geometry's envelope to cluster nearby geometries. The union is then evaluated
	 * as a balanced binary-tree reduction in parallel.
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
	 * @param geoms a non-null list of geometries (intended for Polygon or
	 *              MultiPolygon)
	 * @return a polygonal geometry representing the union (Polygon or MultiPolygon)
	 */
	public static Geometry union(Collection<? extends Geometry> geoms) {
		if (geoms.isEmpty()) {
			return new GeometryFactory().createEmpty(2);
		}
		Geometry[] items = geoms.toArray(new Geometry[0]);
		sort(items, HilbertCode.level(items.length)); // sort according to center point of MBR
		GeometryFactory factory = items[0].getFactory();

		Geometry result = new UnionTask(items, 0, items.length, factory).invoke();
		return toPolygonal(result, factory);
	}

	/**
	 * Balanced binary-tree union over a slice of the Hilbert-sorted array.
	 * Fork/join keeps operand sizes comparable at every merge, unlike a sequential
	 * left-fold where one operand grows monotonically.
	 */
	private static final class UnionTask extends RecursiveTask<Geometry> {
		private static final long serialVersionUID = 1L;

		private final Geometry[] items;
		private final int lo, hi;
		private final GeometryFactory factory;

		UnionTask(Geometry[] items, int lo, int hi, GeometryFactory factory) {
			this.items = items;
			this.lo = lo;
			this.hi = hi;
			this.factory = factory;
		}

		@Override
		protected Geometry compute() {
			int len = hi - lo;
			if (len == 1) {
				return items[lo];
			}
			int mid = (lo + hi) >>> 1;
			if (len <= FORK_THRESHOLD) {
				Geometry left = new UnionTask(items, lo, mid, factory).compute();
				Geometry right = new UnionTask(items, mid, hi, factory).compute();
				return unionPair(left, right, factory);
			}
			UnionTask rightTask = new UnionTask(items, mid, hi, factory);
			rightTask.fork();
			Geometry left = new UnionTask(items, lo, mid, factory).compute();
			Geometry right = rightTask.join();
			return unionPair(left, right, factory);
		}
	}

	/**
	 * Unions two polygonal geometries, avoiding overlay work where envelopes prove
	 * it unnecessary: fully disjoint operands are combined structurally, and for
	 * partially overlapping operands only the elements inside the mutual envelope
	 * intersection are overlaid (mirroring CascadedPolygonUnion's optimization).
	 */
	private static Geometry unionPair(Geometry a, Geometry b, GeometryFactory factory) {
		Envelope envA = a.getEnvelopeInternal();
		Envelope envB = b.getEnvelopeInternal();
		if (!envA.intersects(envB)) {
			return combine(a, b, factory);
		}

		Envelope common = envA.intersection(envB);
		List<Polygon> inA = new ArrayList<>();
		List<Polygon> inB = new ArrayList<>();
		List<Polygon> out = new ArrayList<>();
		splitByEnvelope(a, common, inA, out);
		splitByEnvelope(b, common, inB, out);

		if (inA.isEmpty() || inB.isEmpty()) {
			// no element of one operand reaches the overlap region; nothing can intersect
			return combine(a, b, factory);
		}

		Geometry overlapUnion = OverlayNGRobust.overlay(buildPolygonal(inA, factory), buildPolygonal(inB, factory), OverlayNG.UNION);
		if (out.isEmpty()) {
			return toPolygonal(overlapUnion, factory);
		}
		PolygonExtracter.getPolygons(overlapUnion, out);
		return buildPolygonal(out, factory);
	}

	/**
	 * Partitions the polygon elements of {@code g} by envelope intersection with
	 * {@code env}, appending to {@code in}/{@code out}.
	 */
	private static void splitByEnvelope(Geometry g, Envelope env, List<Polygon> in, List<Polygon> out) {
		for (int i = 0; i < g.getNumGeometries(); i++) {
			Geometry part = g.getGeometryN(i);
			if (part instanceof Polygon) {
				if (part.getEnvelopeInternal().intersects(env)) {
					in.add((Polygon) part);
				} else {
					out.add((Polygon) part);
				}
			}
		}
	}

	/** Structurally merges two polygonal geometries known not to overlap. */
	private static Geometry combine(Geometry a, Geometry b, GeometryFactory factory) {
		List<Polygon> polygons = new ArrayList<>(a.getNumGeometries() + b.getNumGeometries());
		PolygonExtracter.getPolygons(a, polygons);
		PolygonExtracter.getPolygons(b, polygons);
		return buildPolygonal(polygons, factory);
	}

	private static Geometry buildPolygonal(List<Polygon> polygons, GeometryFactory factory) {
		if (polygons.isEmpty()) {
			return factory.createPolygon();
		}
		if (polygons.size() == 1) {
			return polygons.get(0);
		}
		return factory.createMultiPolygon(GeometryFactory.toPolygonArray(polygons));
	}

	/** Discards any non-polygonal overlay artifacts. */
	private static Geometry toPolygonal(Geometry g, GeometryFactory factory) {
		if (g instanceof Polygonal) {
			return g;
		}
		@SuppressWarnings("unchecked")
		List<Polygon> polygons = PolygonExtracter.getPolygons(g);
		return buildPolygonal(polygons, factory);
	}

	/**
	 * Sorts an array of {@link Geometry} objects in-place by their spatial order
	 * using Hilbert curve encoding of their envelopes.
	 *
	 * @param geoms the array of geometries to sort
	 * @param level the resolution level for Hilbert curve encoding
	 */
	private static void sort(Geometry[] geoms, int level) {
		int n = geoms.length;
		if (n < 2) {
			return;
		}

		Envelope globalExtent = new Envelope();
		for (Geometry g : geoms) {
			globalExtent.expandToInclude(g.getEnvelopeInternal());
		}

		HilbertEncoder encoder = new HilbertEncoder(level, globalExtent);
		// pack (key, index) into a long so an unboxed primitive sort suffices
		long[] packed = new long[n];
		for (int i = 0; i < n; i++) {
			packed[i] = ((long) encoder.encode(geoms[i].getEnvelopeInternal()) << 32) | i;
		}
		Arrays.sort(packed);

		Geometry[] sorted = new Geometry[n];
		for (int i = 0; i < n; i++) {
			sorted[i] = geoms[(int) packed[i]];
		}
		System.arraycopy(sorted, 0, geoms, 0, n);
	}
}
