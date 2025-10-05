package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.util.LineStringExtracter;
import org.locationtech.jts.geom.util.PolygonExtracter;
import org.locationtech.jts.index.SpatialIndex;
import org.locationtech.jts.index.hprtree.HPRtree;

/**
 * Endpoint-only snapper with optional polygon-vertex anchoring.
 *
 * <p>
 * Closes small endpoint gaps in near-coverage linework so edges meet at nodes,
 * without moving polygon vertices or interior line vertices. Operates on
 * endpoints only (not a general vertex/edge snapper).
 * </p>
 *
 * <h3>Behavior</h3>
 * <ul>
 * <li>Clusters LineString endpoints within a tolerance (transitively) and snaps
 * them to a representative.</li>
 * <li>If snap-to-valid is enabled, polygon vertices are included as immutable
 * anchors: any cluster containing a polygon vertex snaps exactly to that
 * vertex; polygon geometries themselves are never modified.</li>
 * <li>If snap-to-valid is disabled, polygons are ignored for snapping; clusters
 * composed only of endpoints snap mutually to their mean coordinate.</li>
 * <li>Closed LineStrings (rings) are treated as polygons and are not
 * moved.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <p>
 * Returns a geometry mirroring the input structure: LineStrings are replaced by
 * snapped copies; Polygons/MultiPolygons and other geometry types are returned
 * unchanged.
 * </p>
 *
 * <h3>Notes</h3>
 * <ul>
 * <li>Intended for near-coverage datasets that should form a valid coverage but
 * don’t quite join exactly (tiny endpoint gaps, near-misses).</li>
 * <li>Tolerance is in the data’s units; proximity is evaluated in XY only (Z/M
 * are preserved by copy).</li>
 * <li>Typically followed by noding and polygonization if faces need to be
 * built.</li>
 * </ul>
 *
 * @author Michael Carleton
 */
public class EndpointSnapper {

	private final double tolerance;
	private GeometryFactory gf;

	public EndpointSnapper(double tolerance) {
		this.tolerance = tolerance;
	}

	/**
	 * Snap only the endpoints of LineStrings. If snapToValid is true, polygon
	 * 
	 * vertices act as anchors: any cluster containing a polygon vertex snaps to
	 * that
	 * 
	 * vertex; polygon geometries themselves are never modified.
	 * 
	 * Returns a Geometry preserving the input's structure, with snapped
	 * LineStrings.
	 */
	public Geometry snapEndpoints(Geometry geom, boolean snapToValid) {
		if (tolerance <= 0) {
			return geom;
		}

		gf = geom.getFactory();
		// 1) Collect LineStrings (to be modified) and polygon vertices (anchors)
		List<LineString> lines = new ArrayList<>();
		List<Coordinate> polygonVerts = new ArrayList<>();
		collectGeoms(geom, lines, polygonVerts, snapToValid);

		// Build nodes: endpoints (movable), and optionally polygon vertices (anchors)
		List<Node> nodes = new ArrayList<>();
		// For mapping line -> its start/end node indices
		int[] startNodeIdx = new int[lines.size()];
		int[] endNodeIdx = new int[lines.size()];
		Arrays.fill(startNodeIdx, -1);
		Arrays.fill(endNodeIdx, -1);

		// Add endpoints
		for (int i = 0; i < lines.size(); i++) {
			LineString ls = lines.get(i);
			int n = ls.getNumPoints();
			if (n < 2) {
				continue;
			}
			Coordinate c0 = ls.getCoordinateN(0);
			Coordinate c1 = ls.getCoordinateN(n - 1);
			int idx0 = nodes.size();
			nodes.add(Node.endpoint(c0, i, 0));
			int idx1 = nodes.size();
			nodes.add(Node.endpoint(c1, i, n - 1));
			startNodeIdx[i] = idx0;
			endNodeIdx[i] = idx1;
		}

		// Add polygon vertices as anchors only if snapToValid=true
		if (snapToValid) {
			for (Coordinate pv : polygonVerts) {
				nodes.add(Node.anchor(pv));
			}
		}

		// Nothing to do?
		if (nodes.isEmpty()) {
			return rebuildPreservingStructure(geom, Collections.emptyMap());
		}

		// 2) Build spatial index over nodes
		SpatialIndex index = new HPRtree();
		for (int i = 0; i < nodes.size(); i++) {
			Envelope env = new Envelope(nodes.get(i).coord);
			env.expandBy(tolerance);
			index.insert(env, i);
		}

		// 3) Cluster with union-find (within tolerance, transitive)
		UnionFind uf = new UnionFind(nodes.size());
		for (int i = 0; i < nodes.size(); i++) {
			Coordinate ci = nodes.get(i).coord;
			Envelope q = new Envelope(ci);
			q.expandBy(tolerance);
			@SuppressWarnings("unchecked")
			List<Integer> cand = index.query(q);
			for (int j : cand) {
				if (j <= i) {
					continue;
				}
				Coordinate cj = nodes.get(j).coord;
				if (ci.distance(cj) <= tolerance) {
					uf.union(i, j);
				}
			}
		}

		// 4) Compute representative coordinate per cluster
		Map<Integer, ClusterStats> stats = new HashMap<>();
		for (int i = 0; i < nodes.size(); i++) {
			int r = uf.find(i);
			ClusterStats s = stats.computeIfAbsent(r, k -> new ClusterStats());
			Node ni = nodes.get(i);
			if (ni.isAnchor) {
				if (!s.hasAnchor) {
					s.hasAnchor = true;
					s.anchorCoord = ni.coord; // pick first polygon vertex as anchor
				}
			} else {
				s.add(ni.coord); // endpoints contribute to mean
			}
		}

		Map<Integer, Coordinate> rep = new HashMap<>();
		for (Map.Entry<Integer, ClusterStats> e : stats.entrySet()) {
			ClusterStats s = e.getValue();
			Coordinate c = s.hasAnchor ? s.anchorCoord : (s.count > 0 ? new Coordinate(s.sx / s.count, s.sy / s.count) : null);
			// If a cluster has only anchors (no endpoints), c will be anchor; if only one
			// endpoint, mean==itself.
			rep.put(e.getKey(), c);
		}

		// 5) Build snapped LineStrings and map originals to replacements
		Map<LineString, LineString> replace = new IdentityHashMap<>();
		for (int li = 0; li < lines.size(); li++) {
			LineString ls = lines.get(li);
			int n = ls.getNumPoints();
			if (n < 2) {
				replace.put(ls, ls);
				continue;
			}
			Coordinate[] coords = ls.getCoordinates();
			Coordinate[] out = new Coordinate[coords.length];
			for (int k = 0; k < coords.length; k++) {
				out[k] = new Coordinate(coords[k]);
			}

			// Start
			int niStart = startNodeIdx[li];
			if (niStart >= 0) {
				int r = uf.find(niStart);
				Coordinate target = rep.get(r);
				if (target != null) {
					out[0] = new Coordinate(target);
				}
			}
			// End
			int niEnd = endNodeIdx[li];
			if (niEnd >= 0) {
				int r = uf.find(niEnd);
				Coordinate target = rep.get(r);
				if (target != null) {
					out[out.length - 1] = new Coordinate(target);
				}
			}

			replace.put(ls, gf.createLineString(out));
		}

		// 6) Rebuild geometry preserving original structure, replacing only LineStrings
		return rebuildPreservingStructure(geom, replace);
	}

	// Collects LineStrings to modify and polygon vertices (if requested) as
	// anchors.
	private void collectGeoms(Geometry root, List<LineString> linesOut, List<Coordinate> polyVertsOut, boolean includePolygonVerts) {
		// 1) Extract non-closed LineStrings only (ignore rings/closed lines, which
		// represent polygons)
		List<LineString> tmpLines = new ArrayList<>();
		LineStringExtracter.getLines(root, tmpLines);
		for (LineString ls : tmpLines) {
			if (!ls.isClosed() && ls.getNumPoints() >= 2) {
				linesOut.add(ls);
			}
		}

		// 2) Optionally extract polygon vertices as anchors
		if (includePolygonVerts) {
			List<Polygon> polys = new ArrayList<>();
			PolygonExtracter.getPolygons(root, polys);
			for (Polygon p : polys) {
				addRingVertices(p.getExteriorRing(), polyVertsOut);
				for (int i = 0; i < p.getNumInteriorRing(); i++) {
					addRingVertices(p.getInteriorRingN(i), polyVertsOut);
				}
			}
		}
	}

	private void addRingVertices(LineString ring, List<Coordinate> out) {
		Coordinate[] cs = ring.getCoordinates();
		for (Coordinate c : cs) {
			out.add(c);
		}
	}

	// Rebuilds geometry tree, replacing only LineStrings per 'replace' map
	private Geometry rebuildPreservingStructure(Geometry g, Map<LineString, LineString> replace) {
		if (g instanceof LineString) {
			LineString ls = (LineString) g;
			LineString repl = replace.get(ls);
			return repl != null ? repl : ls;
		}
		if (g instanceof MultiLineString) {
			MultiLineString mls = (MultiLineString) g;
			LineString[] arr = new LineString[mls.getNumGeometries()];
			for (int i = 0; i < arr.length; i++) {
				LineString ls = (LineString) mls.getGeometryN(i);
				LineString repl = replace.get(ls);
				arr[i] = repl != null ? repl : ls;
			}
			return gf.createMultiLineString(arr);
		}
		if (g instanceof Polygon || g instanceof MultiPolygon || g instanceof Point || g instanceof MultiPoint) {
			// Return unchanged
			return g;
		}
		if (g instanceof GeometryCollection) {
			GeometryCollection gc = (GeometryCollection) g;
			Geometry[] parts = new Geometry[gc.getNumGeometries()];
			for (int i = 0; i < parts.length; i++) {
				parts[i] = rebuildPreservingStructure(gc.getGeometryN(i), replace);
			}
			return gf.createGeometryCollection(parts);
		}
		return g;
	}

	// Helper classes
	private static class Node {
		final Coordinate coord;
		final boolean isAnchor; // polygon vertex

		private Node(Coordinate c, boolean anchor) {
			this.coord = c;
			this.isAnchor = anchor;
		}

		static Node endpoint(Coordinate c, int lineIdx, int coordIdx) {
			return new Node(c, false);
		}

		static Node anchor(Coordinate c) {
			return new Node(c, true);
		}
	}

	private static class ClusterStats {
		double sx = 0, sy = 0;
		int count = 0;
		boolean hasAnchor = false;
		Coordinate anchorCoord = null;

		void add(Coordinate c) {
			sx += c.x;
			sy += c.y;
			count++;
		}
	}

	private static class UnionFind {
		private final int[] p;
		private final int[] r;

		UnionFind(int n) {
			p = new int[n];
			r = new int[n];
			for (int i = 0; i < n; i++) {
				p[i] = i;
			}
		}

		int find(int x) {
			return p[x] == x ? x : (p[x] = find(p[x]));
		}

		void union(int a, int b) {
			int ra = find(a), rb = find(b);
			if (ra == rb) {
				return;
			}
			if (r[ra] < r[rb]) {
				p[ra] = rb;
			} else if (r[rb] < r[ra]) {
				p[rb] = ra;
			} else {
				p[rb] = ra;
				r[ra]++;
			}
		}
	}
}