package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.locationtech.jts.algorithm.Distance;
import org.locationtech.jts.algorithm.locate.PointOnGeometryLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.geom.util.LinearComponentExtracter;

import com.github.micycle1.geoblitz.HPRtreeX.DistanceToItem;

/**
 * Computes distances from query points to the linear components of one or two
 * input geometries, using a spatial index for fast repeated queries.
 * <p>
 * The indexed items are not individual segments; instead, each
 * {@code LineString} is partitioned into contiguous “facets” (polyline chunks)
 * consisting of {@code N} coordinates (and therefore {@code N-1} segments).
 * Each facet is inserted into an {@code HPRtree} using its envelope, and
 * nearest-facet search is performed with an exact point-to-segment distance
 * refinement. This reduces index item count compared to segment-level indexing
 * while keeping good query speed.
 * </p>
 *
 * <h2>Signed vs unsigned distance</h2>
 * <ul>
 * <li>If a polygonal {@code boundary} is provided <em>and</em> a
 * {@link PointOnGeometryLocator} is enabled, {@link #distance(Coordinate)}
 * returns a <b>signed</b> distance: negative for points in the
 * {@link Location#EXTERIOR} of the boundary, and positive/zero for
 * {@link Location#INTERIOR}/{@link Location#BOUNDARY}.</li>
 * <li>{@link #unsignedDistance(Coordinate)} always returns a non-negative
 * distance.</li>
 * </ul>
 *
 * <h2>Inputs</h2>
 * <ul>
 * <li>{@code boundary}: must be {@link Polygonal} (Polygon/MultiPolygon). Its
 * linear components are also used as distance targets (e.g., polygon rings).
 * Additionally, it may be used to determine the sign of
 * {@link #distance(Coordinate)}.</li>
 * <li>{@code obstacles}: may be any {@link Geometry}; only its linear
 * components (LineStrings) are extracted and used as distance targets.</li>
 * </ul>
 * Either {@code boundary} or {@code obstacles} may be {@code null}, but not
 * both.
 *
 * @author Michael Carleton
 */
public class IndexedLinearDistance {

	private static final int DEFAULT_FACET_COORD_COUNT = 8;

	private final Geometry boundary;
	private final Geometry obstacles;
	private final PointOnGeometryLocator boundaryLocator; // nullable if boundary == null

	private final int facetCoordCount; // N (>= 2)
	private final List<Facet> facets = new ArrayList<>();
	private final HPRtreeX<Facet> tree = new HPRtreeX<>();

	/**
	 * Constructs an index using only a boundary geometry.
	 * <p>
	 * If {@code boundary} is non-null, signed distance is enabled using an
	 * {@link YStripesPointInAreaLocator}. Facet size defaults to
	 * {@value #DEFAULT_FACET_COORD_COUNT} coordinates.
	 * </p>
	 *
	 * @param boundary polygonal boundary used both as a distance target (via its
	 *                 rings) and to determine the sign of
	 *                 {@link #distance(Coordinate)}
	 * @throws IllegalArgumentException if {@code boundary} is non-polygonal
	 */
	public IndexedLinearDistance(Geometry boundary) {
		this(boundary, null, DEFAULT_FACET_COORD_COUNT);
	}

	/**
	 * Constructs an index from a boundary and an obstacles geometry.
	 * <p>
	 * If {@code boundary} is non-null, signed distance is enabled using an
	 * {@link YStripesPointInAreaLocator}. Facet size defaults to
	 * {@value #DEFAULT_FACET_COORD_COUNT} coordinates.
	 * </p>
	 *
	 * @param boundary  polygonal boundary used for distance targets and
	 *                  (optionally) sign
	 * @param obstacles any geometry; only linear components are used as distance
	 *                  targets
	 * @throws IllegalArgumentException if both inputs are {@code null}
	 * @throws IllegalArgumentException if {@code boundary} is non-polygonal
	 */
	public IndexedLinearDistance(Geometry boundary, Geometry obstacles) {
		this(boundary, obstacles, DEFAULT_FACET_COORD_COUNT);
	}

	/**
	 * Constructs an index from a boundary and/or obstacles geometry, using an
	 * explicit facet size.
	 * <p>
	 * If {@code boundary} is non-null, signed distance is enabled using an
	 * {@link YStripesPointInAreaLocator}.
	 * </p>
	 *
	 * @param boundary        polygonal boundary used for distance targets and
	 *                        (optionally) sign
	 * @param obstacles       any geometry; only linear components are used as
	 *                        distance targets
	 * @param facetCoordCount number of coordinates per facet (must be {@code >= 2})
	 * @throws IllegalArgumentException if both inputs are {@code null}
	 * @throws IllegalArgumentException if {@code boundary} is non-polygonal
	 * @throws IllegalArgumentException if {@code facetCoordCount < 2}
	 */
	public IndexedLinearDistance(Geometry boundary, Geometry obstacles, int facetCoordCount) {
		this(boundary, boundary == null ? null : new YStripesPointInAreaLocator(boundary), obstacles, facetCoordCount);
	}

	/**
	 * Constructs an index with an explicitly supplied boundary locator.
	 * <p>
	 * Passing a null {@code boundaryLocator} disables signed distance even if
	 * {@code boundary} is non-null, causing {@link #distance(Coordinate)} to behave
	 * like {@link #unsignedDistance(Coordinate)}.
	 * </p>
	 *
	 * @param boundary        polygonal boundary used for distance targets and
	 *                        (optionally) sign
	 * @param boundaryLocator locator used to compute sign; may be null to disable
	 *                        signed distance
	 * @param obstacles       any geometry; only linear components are used as
	 *                        distance targets
	 * @param facetCoordCount number of coordinates per facet (must be {@code >= 2})
	 * @throws IllegalArgumentException if both inputs are {@code null}
	 * @throws IllegalArgumentException if {@code boundary} is non-polygonal
	 * @throws IllegalArgumentException if {@code facetCoordCount < 2}
	 */
	public IndexedLinearDistance(Geometry boundary, PointOnGeometryLocator boundaryLocator, Geometry obstacles, int facetCoordCount) {
		if (boundary == null && obstacles == null)
			throw new IllegalArgumentException("Either boundary or obstacles must be non-null");
		if (boundary != null && !(boundary instanceof Polygonal))
			throw new IllegalArgumentException("boundary must be polygonal (Polygon/MultiPolygon)");
		if (facetCoordCount < 2)
			throw new IllegalArgumentException("facetCoordCount must be >= 2");

		this.boundary = boundary;
		this.obstacles = obstacles;
		this.boundaryLocator = boundaryLocator;
		this.facetCoordCount = facetCoordCount;

		addFacetsFrom(boundary);
		addFacetsFrom(obstacles);

		for (Facet f : facets) {
			tree.insert(f.env, f);
		}
		tree.build();
	}

	/**
	 * Returns the (optionally) signed distance from the query coordinate to the
	 * nearest indexed linear component.
	 * <ul>
	 * <li>If {@code boundaryLocator} is present: negative means the point is
	 * outside the boundary.</li>
	 * <li>If {@code boundaryLocator} is absent: this equals
	 * {@link #unsignedDistance(Coordinate)}.</li>
	 * </ul>
	 *
	 * @param p query coordinate
	 * @return signed distance, unsigned distance, or {@link Double#NaN} if no
	 *         facets exist
	 */
	public double distance(Coordinate p) {
		double d = unsignedDistance(p);
		if (boundaryLocator == null)
			return d;

		int loc = boundaryLocator.locate(p);
		return (loc == Location.EXTERIOR) ? -d : d; // INTERIOR/BOUNDARY => positive/0
	}

	/**
	 * Convenience overload of {@link #distance(Coordinate)}.
	 *
	 * @param x query x
	 * @param y query y
	 * @return signed distance, unsigned distance, or {@link Double#NaN} if no
	 *         facets exist
	 */
	public double distance(double x, double y) {
		return distance(new Coordinate(x, y));
	}

	/**
	 * Returns the non-negative distance from the query coordinate to the nearest
	 * indexed linear component.
	 * <p>
	 * This method uses the spatial index nearest-neighbour search and then computes
	 * the exact distance to segments within the chosen facet.
	 * </p>
	 *
	 * @param p query coordinate
	 * @return distance {@code >= 0}, or {@link Double#NaN} if no facets exist
	 */
	public double unsignedDistance(Coordinate p) {
		if (facets.isEmpty()) {
			return Double.NaN;
		}

		final DistanceToItem<Facet> fn = (coord, f) -> f.pointDistanceSq(coord);
		Facet nn = tree.nearestNeighborSq(p, fn);
		if (nn == null) {
			return Double.NaN;
		}

		return Math.sqrt(nn.pointDistanceSq(p));
	}

	/**
	 * Convenience overload of {@link #unsignedDistance(Coordinate)}.
	 *
	 * @param x query x
	 * @param y query y
	 * @return distance {@code >= 0}, or {@link Double#NaN} if no facets exist
	 */
	public double unsignedDistance(double x, double y) {
		return unsignedDistance(new Coordinate(x, y));
	}

	private void addFacetsFrom(Geometry g) {
		if (g == null || g.isEmpty()) {
			return;
		}

		@SuppressWarnings("unchecked")
		Collection<LineString> lines = LinearComponentExtracter.getLines(g);

		final int step = facetCoordCount - 1; // advance by segments, keep 1 coord overlap between facets

		for (LineString ls : lines) {
			Coordinate[] cs = ls.getCoordinates();
			if (cs == null || cs.length < 2) {
				continue;
			}

			// Create non-overlapping segment coverage, with shared end/start coordinate
			// between facets.
			// Facet covers coords [i .. end], and we advance i += step.
			for (int i = 0; i < cs.length - 1; i += step) {
				int end = Math.min(i + facetCoordCount - 1, cs.length - 1);
				if (end <= i)
					continue;

				// Optional: skip facets that are entirely degenerate (all coords identical)
				// Quick check: ensure at least one segment has differing endpoints.
				boolean hasNonZeroSeg = false;
				for (int k = i; k < end; k++) {
					if (!cs[k].equals2D(cs[k + 1])) {
						hasNonZeroSeg = true;
						break;
					}
				}
				if (!hasNonZeroSeg)
					continue;

				facets.add(new Facet(cs, i, end));
			}
		}
	}

	/**
	 * A polyline chunk covering coordinates {@code [start..end]} inclusive (i.e.
	 * segments {@code start..end-1}), along with its envelope for indexing.
	 */
	private static final class Facet {
		final Coordinate[] cs;
		final int start; // inclusive
		final int end; // inclusive
		final Envelope env;

		Facet(Coordinate[] cs, int start, int end) {
			this.cs = cs;
			this.start = start;
			this.end = end;
			this.env = computeEnvelope(cs, start, end);
		}

		double pointDistanceSq(Coordinate p) {
			double min = Double.POSITIVE_INFINITY;
			for (int i = start; i < end; i++) {
				Coordinate a = cs[i];
				Coordinate b = cs[i + 1];
				double d2 = Distance.pointToSegmentSq(p, a, b);
				if (d2 < min) {
					min = d2;
					if (min == 0) {
						return 0;
					}
				}
			}
			return min;
		}

		private static Envelope computeEnvelope(Coordinate[] cs, int start, int end) {
			double minX = cs[start].x, maxX = cs[start].x;
			double minY = cs[start].y, maxY = cs[start].y;
			for (int i = start + 1; i <= end; i++) {
				Coordinate c = cs[i];
				if (c.x < minX)
					minX = c.x;
				else if (c.x > maxX)
					maxX = c.x;
				if (c.y < minY)
					minY = c.y;
				else if (c.y > maxY)
					maxY = c.y;
			}
			return new Envelope(minX, maxX, minY, maxY);
		}
	}
}