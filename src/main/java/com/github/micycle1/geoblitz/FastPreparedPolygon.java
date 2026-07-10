package com.github.micycle1.geoblitz;

import java.util.Collection;
import java.util.Collections;

import org.locationtech.jts.algorithm.locate.PointOnGeometryLocator;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygonal;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.geom.prep.PreparedPolygon;
import org.locationtech.jts.noding.FastSegmentSetIntersectionFinder;
import org.locationtech.jts.noding.SegmentIntersectionDetector;
import org.locationtech.jts.noding.SegmentSetMutualIntersector;
import org.locationtech.jts.noding.SegmentString;
import org.locationtech.jts.noding.SegmentStringUtil;

/**
 * A drop-in replacement for JTS {@link PreparedPolygon} with faster spatial
 * predicates.
 * <p>
 * This class reuses the (battle-tested) prepared-predicate logic of JTS but
 * swaps in faster GeoBlitz components:
 * <ul>
 * <li>{@link YStripesPointInAreaLocator} instead of
 * {@code IndexedPointInAreaLocator} for point-in-area tests;</li>
 * <li>{@link FastMCIndexSegmentSetMutualIntersector} (HPR-tree based) instead
 * of {@code MCIndexSegmentSetMutualIntersector} (STR-tree based) for
 * boundary-segment intersection tests.</li>
 * </ul>
 * Both components are built eagerly during construction, so instances are
 * immediately safe for concurrent read-only use and queries avoid the
 * per-call {@code synchronized} lazy initialization of the JTS class.
 * <p>
 * Intended for build-once, query-many workloads: repeated
 * {@code intersects}/{@code contains}/{@code covers}/{@code containsProperly}
 * (etc.) tests of many geometries against a fixed polygonal geometry.
 *
 * @author Michael Carleton
 */
public class FastPreparedPolygon extends PreparedPolygon {

	private final FastSegmentSetIntersectionFinder segIntFinder;
	private final PointOnGeometryLocator pointLocator;

	/**
	 * Creates a prepared polygon for the given polygonal geometry. All internal
	 * indexes are built eagerly, so the instance is immediately ready for
	 * concurrent read-only queries.
	 *
	 * @param poly the polygonal geometry to prepare (must not be null)
	 */
	public FastPreparedPolygon(Polygonal poly) {
		super(poly);
		pointLocator = new YStripesPointInAreaLocator(getGeometry());
		@SuppressWarnings("unchecked")
		Collection<SegmentString> segStrings = SegmentStringUtil.extractSegmentStrings(getGeometry());
		segIntFinder = new HPRSegmentSetIntersectionFinder(segStrings);
	}

	@Override
	public FastSegmentSetIntersectionFinder getIntersectionFinder() {
		return segIntFinder;
	}

	@Override
	public PointOnGeometryLocator getPointLocator() {
		return pointLocator;
	}

	/**
	 * Prepares a geometry for repeated spatial predicate evaluation. Polygonal
	 * geometries are prepared with {@link FastPreparedPolygon}; other geometry
	 * types fall back to the standard JTS {@link PreparedGeometryFactory}.
	 *
	 * @param geom the geometry to prepare (must not be null)
	 * @return a PreparedGeometry over the input
	 */
	public static PreparedGeometry prepare(Geometry geom) {
		if (geom instanceof Polygonal) {
			return new FastPreparedPolygon((Polygonal) geom);
		}
		return PreparedGeometryFactory.prepare(geom);
	}

	/**
	 * A {@link FastSegmentSetIntersectionFinder} whose segment index is a
	 * {@link FastMCIndexSegmentSetMutualIntersector} (backed by an HPR-tree)
	 * rather than the STR-tree-based JTS default. The superclass is initialized
	 * with an empty segment set so its internal index is never populated; all
	 * queries are routed to the HPR-tree index via the overridden methods.
	 */
	private static final class HPRSegmentSetIntersectionFinder extends FastSegmentSetIntersectionFinder {

		private final FastMCIndexSegmentSetMutualIntersector segSetMutInt;

		HPRSegmentSetIntersectionFinder(Collection<SegmentString> segStrings) {
			super(Collections.emptyList());
			segSetMutInt = new FastMCIndexSegmentSetMutualIntersector(segStrings);
		}

		@Override
		public SegmentSetMutualIntersector getSegmentSetIntersector() {
			return segSetMutInt;
		}

		@Override
		public boolean intersects(@SuppressWarnings("rawtypes") Collection segStrings) {
			return intersects(segStrings, new SegmentIntersectionDetector());
		}

		@Override
		public boolean intersects(@SuppressWarnings("rawtypes") Collection segStrings, SegmentIntersectionDetector intDetector) {
			segSetMutInt.process(segStrings, intDetector);
			return intDetector.hasIntersection();
		}
	}
}
