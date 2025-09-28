package com.github.micycle1.geoblitz;

import java.util.Arrays;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.linearref.LengthIndexedLine;
import org.locationtech.jts.linearref.LinearGeometryBuilder;
import org.locationtech.jts.linearref.LinearIterator;
import org.locationtech.jts.linearref.LinearLocation;

/**
 * IndexedLengthIndexedLine provides length-based linear-referencing using a
 * prebuilt cumulative-length index so that length→location queries are answered
 * in O(log n) time.
 * <p>
 * In contrast, the vanilla {@link LengthIndexedLine} locates positions by
 * scanning segments on each query (linear time), so this class is faster for
 * repeated queries on large geometries.
 * </p>
 * <p>
 * <b>When to prefer this class</b>
 * </p>
 * <ul>
 * <li>Use {@code IndexedLengthIndexedLine} when you have a large linear
 * geometry and perform many length-based queries (extractPoint, extractLine,
 * etc.).</li>
 * <li>Use the vanilla {@link LengthIndexedLine} when queries are infrequent or
 * when minimizing memory overhead is paramount.</li>
 * </ul>
 *
 * <p>
 * <b>Caveats</b>
 * </p>
 * <ul>
 * <li>Components are treated in geometry order (as with
 * {@code LengthIndexedLine}); they need not be spatially contiguous.</li>
 * <li>The class preserves the endpoint-resolution and degenerate-segment
 * handling of {@code LengthLocationMap}; care is taken for zero-length
 * segments/components.</li>
 * </ul>
 *
 * @author Michael Carleton
 */
public final class IndexedLengthIndexedLine {
	private final Geometry linearGeom;

	private final int numComps;
	// compPrefixLen[i] = total length of components [0, i)
	private final double[] compPrefixLen;
	// segPrefixLen[c][s] = cumulative length within component c at start of segment
	// s
	// size = (numSegs(c) + 1), with segPrefixLen[c][numSegs(c)] == length of
	// component c
	private final double[][] segPrefixLen;

	private final double totalLength;

	/**
	 * Constructs an index for the supplied linear geometry.
	 * <p>
	 * The index is typically built once and used for many queries.
	 * </p>
	 *
	 * @param linearGeom the linear {@link Geometry} (LineString or MultiLineString)
	 */
	public IndexedLengthIndexedLine(Geometry linearGeom) {
		this.linearGeom = linearGeom;
		this.numComps = linearGeom.getNumGeometries();

		this.compPrefixLen = new double[numComps + 1];
		this.segPrefixLen = new double[numComps][];
		double total = 0.0;

		for (int c = 0; c < numComps; c++) {
			LineString ls = (LineString) linearGeom.getGeometryN(c);
			int nPts = ls.getNumPoints();
			int nSeg = Math.max(0, nPts - 1);

			double[] prefix = new double[nSeg + 1];
			prefix[0] = 0.0;

			double compLen = 0.0;
			for (int s = 0; s < nSeg; s++) {
				Coordinate p0 = ls.getCoordinateN(s);
				Coordinate p1 = ls.getCoordinateN(s + 1);
				double len = p0.distance(p1);
				compLen += len;
				prefix[s + 1] = compLen;
			}
			this.segPrefixLen[c] = prefix;

			compPrefixLen[c] = total;
			total += compLen;
		}
		compPrefixLen[numComps] = total;
		this.totalLength = total;
	}

	/**
	 *
	 * Computes the {@link Coordinate} for the point on the line at the given length
	 * index.
	 * <p>
	 * Negative indices are interpreted as measured from the end of the geometry.
	 * Out-of-range values are clamped to the line endpoints.
	 * </p>
	 *
	 * @param index the length index along the line
	 * @return the coordinate at the given index
	 */
	public Coordinate extractPoint(double index) {
		LinearLocation loc = getLocation(index, true);
		return loc.getCoordinate(linearGeom);
	}

	/**
	 *
	 * Computes the {@link Coordinate} for the point on the line at the given length
	 * index, offset laterally by the supplied distance.
	 * <p>
	 * A positive {@code offsetDistance} offsets the point to the left of the
	 * oriented segment, negative to the right.
	 * </p>
	 *
	 * @param index          the length index along the line
	 * @param offsetDistance lateral offset distance
	 * @return the offset coordinate
	 */
	public Coordinate extractPoint(double index, double offsetDistance) {
		LinearLocation loc = getLocation(index, true);
		LinearLocation locLow = loc.toLowest(linearGeom);
		return locLow.getSegment(linearGeom).pointAlongOffset(locLow.getSegmentFraction(), offsetDistance);
	}

	/**
	 *
	 * Extracts the subline between two length indices.
	 * <p>
	 * If {@code endIndex < startIndex} the returned geometry has reversed
	 * orientation. Indices are clamped to the valid range.
	 * </p>
	 *
	 * @param startIndex start length index
	 * @param endIndex   end length index
	 * @return a linear {@link Geometry} representing the subline
	 */
	public Geometry extractLine(double startIndex, double endIndex) {
		double startIndex2 = clampIndex(startIndex);
		double endIndex2 = clampIndex(endIndex);
		boolean resolveStartLower = startIndex2 == endIndex2;
		LinearLocation startLoc = getLocation(startIndex2, resolveStartLower);
		LinearLocation endLoc = getLocation(endIndex2, true);
		return extractSubline(linearGeom, startLoc, endLoc);
	}

	public double getStartIndex() {
		return 0.0;
	}

	public double getEndIndex() {
		return totalLength;
	}

	public boolean isValidIndex(double index) {
		double pos = positiveIndex(index);
		return pos >= 0.0 && pos <= totalLength;
	}

	public double clampIndex(double index) {
		double pos = positiveIndex(index);
		if (pos < 0.0) {
			return 0.0;
		}
		if (pos > totalLength) {
			return totalLength;
		}
		return pos;
	}

	private double positiveIndex(double index) {
		if (index >= 0.0) {
			return index;
		}
		return totalLength + index;
	}

	// Fast mapping from length -> LinearLocation
	private LinearLocation getLocation(double index, boolean resolveLower) {
		double length = index;
		if (length < 0.0) {
			length = totalLength + length;
		}
		if (length <= 0.0) {
			return new LinearLocation();
		}
		if (length >= totalLength) {
			LinearLocation end = LinearLocation.getEndLocation(linearGeom);
			if (resolveLower) {
				return end;
			}
			return resolveHigher(end);
		}

		int comp = findComponent(length);
		double compStart = compPrefixLen[comp];

		// When index equals the start of a component (not the first), resolve to
		// previous endpoint
		if (length == compStart && comp > 0) {
			LinearLocation loc = new LinearLocation(comp - 1, numSegments(comp - 1), 0.0);
			return resolveLower ? loc : resolveHigher(loc);
		}

		double local = length - compStart;
		LinearLocation loc = locateWithinComponent(comp, local);
		return resolveLower ? loc : resolveHigher(loc);
	}

	private LinearLocation locateWithinComponent(int comp, double local) {
		double[] prefix = segPrefixLen[comp];
		int nSeg = numSegments(comp);

		if (local <= 0.0) {
			return new LinearLocation(comp, 0, 0.0);
		}
		if (local >= prefix[nSeg]) {
			return new LinearLocation(comp, nSeg, 0.0);
		}

		int idx = Arrays.binarySearch(prefix, local);
		if (idx >= 0) {
			// Exact boundary within component -> endpoint at segment index idx
			return new LinearLocation(comp, idx, 0.0);
		}
		int insert = -idx - 1;
		int segIndex = insert - 1;
		double segStart = prefix[segIndex];
		double segEnd = prefix[insert];
		double segLen = segEnd - segStart;

		double frac;
		if (segLen <= 0.0) {
			frac = 0.0;
		} else {
			frac = (local - segStart) / segLen;
			if (frac <= 0.0) {
				frac = 0.0;
			} else if (frac >= 1.0) {
				segIndex += 1;
				frac = 0.0;
			}
		}
		return new LinearLocation(comp, segIndex, frac);
	}

	private LinearLocation resolveHigher(LinearLocation loc) {
		if (!loc.isEndpoint(linearGeom)) {
			return loc;
		}
		int compIndex = loc.getComponentIndex();
		if (compIndex >= numComps - 1) {
			return loc;
		}
		do {
			compIndex++;
		} while (compIndex < numComps - 1 && lengthOfComponent(compIndex) == 0.0);
		return new LinearLocation(compIndex, 0, 0.0);
	}

	private int numSegments(int comp) {
		LineString ls = (LineString) linearGeom.getGeometryN(comp);
		int nPts = ls.getNumPoints();
		return Math.max(0, nPts - 1);
	}

	private double lengthOfComponent(int comp) {
		double[] prefix = segPrefixLen[comp];
		return prefix.length == 0 ? 0.0 : prefix[prefix.length - 1];
	}

	/**
	 * Finds component index such that compPrefixLen[comp] <= length <
	 * compPrefixLen[comp+1].
	 */
	private int findComponent(double length) {
		int idx = Arrays.binarySearch(compPrefixLen, length);
		if (idx >= 0) {
			return Math.min(Math.max(0, idx), numComps - 1);
		}
		int insert = -idx - 1;
		int comp = insert - 1;
		if (comp < 0) {
			comp = 0;
		}
		if (comp >= numComps) {
			comp = numComps - 1;
		}
		return comp;
	}

	private Geometry extractSubline(Geometry line, LinearLocation start, LinearLocation end) {
		if (end.compareTo(start) < 0) {
			Geometry g = computeLinearSubline(line, end, start);
			return g.reverse();
		}
		return computeLinearSubline(line, start, end);
	}

	/**
	 * Assumes start <= end. Builds a linear geometry (LineString or
	 * MultiLineString) representing the subline between the two locations.
	 */
	private Geometry computeLinearSubline(Geometry line, LinearLocation start, LinearLocation end) {
		LinearGeometryBuilder builder = new LinearGeometryBuilder(line.getFactory());
		builder.setFixInvalidLines(true);
		if (!start.isVertex()) {
			builder.add(start.getCoordinate(line));
		}

		for (LinearIterator it = new LinearIterator(line, start); it.hasNext(); it.next()) {
			if (end.compareLocationValues(it.getComponentIndex(), it.getVertexIndex(), 0.0) < 0) {
				break;
			}
			Coordinate pt = it.getSegmentStart();
			builder.add(pt);
			if (it.isEndOfLine()) {
				builder.endLine();
			}
		}

		if (!end.isVertex()) {
			builder.add(end.getCoordinate(line));
		}

		return builder.getGeometry();
	}
}
