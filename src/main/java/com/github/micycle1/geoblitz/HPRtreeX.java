package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.List;

import java.util.Comparator;
import java.util.PriorityQueue;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.index.hprtree.HilbertEncoder;
import org.locationtech.jts.util.IntArrayList;

/**
 * Improves the original JTS HPRTree implementation with eXtended features:
 * <ul>
 * <li>Generics support so the tree can store arbitrary user objects of type
 * {@code T}.</li>
 * <li>Support for an <em>early-exit</em> item visitor to terminate spatial
 * queries as soon as a condition is satisfied.</li>
 * <li>Efficient nearest-neighbor search using a best-first traversal with
 * bounding-box pruning.</li>
 * </ul>
 * <p>
 * Thread-safety notes:
 * <ul>
 * <li>{@link #build()} is idempotent and uses internal synchronization so it is
 * safe to call concurrently; the index will be prepared once.</li>
 * <li>After the index is built, read/query operations are safe (no
 * modifications are performed). However, concurrent calls that attempt to
 * insert after a build will throw an exception.</li>
 * </ul>
 *
 * @param <T> the type of the items stored in the tree
 */
public class HPRtreeX<T> {

	private static final int ENV_SIZE = 4;
	private static final int HILBERT_LEVEL = 12;
	private static final int DEFAULT_NODE_CAPACITY = 16;

	private List<Item<T>> itemsToLoad = new ArrayList<>();
	private final int nodeCapacity;
	private int numItems = 0;
	private final Envelope totalExtent = new Envelope();
	private int[] layerStartIndex;
	private double[] nodeBounds;
	private double[] itemBounds;
	private Object[] itemValues; // Java cannot create a generic array of T
	private volatile boolean isBuilt = false;

	// Used to support early termination during query
	private boolean terminated = false;

	/**
	 * Creates a new index using the default node capacity.
	 *
	 * <p>
	 * The default node capacity is {@code DEFAULT_NODE_CAPACITY}. Newly created
	 * trees are initially empty and buffer inserted items in memory until the tree
	 * is built (explicitly by calling {@link #build()} or implicitly when
	 * performing a query).
	 * </p>
	 */
	public HPRtreeX() {
		this(DEFAULT_NODE_CAPACITY);
	}

	/**
	 * Creates a new index with the specified node capacity.
	 *
	 * @param nodeCapacity the number of children stored per internal node (controls
	 *                     tree fanout and memory layout). Must be &gt; 0. Larger
	 *                     values generally produce shallower trees but increase the
	 *                     per-node scan cost.
	 */
	public HPRtreeX(int nodeCapacity) {
		this.nodeCapacity = nodeCapacity;
	}

	/**
	 * Returns the number of items that have been inserted into this tree.
	 *
	 * <p>
	 * This count includes items that have been added but not yet bulk-loaded into
	 * the internal arrays (i.e. before {@link #build()} is called).
	 * </p>
	 *
	 * @return the number of inserted items
	 */
	public int size() {
		return numItems;
	}

	/**
	 * Inserts an item with its bounding envelope into the tree's insertion buffer.
	 *
	 * <p>
	 * Items may be inserted only until the index is built. Calling this method
	 * after {@link #build()} has been invoked will throw an
	 * {@link IllegalStateException}.
	 * </p>
	 *
	 * @param itemEnv the axis-aligned bounding envelope of the item (must not be
	 *                {@code null})
	 * @param item    the item to store
	 * @throws IllegalStateException if called after the tree has been built
	 */
	public void insert(Envelope itemEnv, T item) {
		if (isBuilt) {
			throw new IllegalStateException("Cannot insert items after tree is built.");
		}
		numItems++;
		itemsToLoad.add(new Item<>(itemEnv, item));
		totalExtent.expandToInclude(itemEnv);
	}

	/**
	 * Returns a list of items whose envelopes intersect the supplied search
	 * envelope.
	 *
	 * <p>
	 * The returned list contains all items that intersect {@code searchEnv}. The
	 * order of items in the returned list follows the internal node/item ordering
	 * (Hilbert-ordered blocks), but is not otherwise guaranteed.
	 * </p>
	 *
	 * <p>
	 * Calling this method will build the index if it has not already been built.
	 * </p>
	 *
	 * @param searchEnv the query envelope
	 * @return a list of matching items; an empty list if no items intersect
	 *         {@code searchEnv}
	 */
	public List<T> query(Envelope searchEnv) {
		build();

		if (!totalExtent.intersects(searchEnv))
			return new ArrayList<>();

		List<T> result = new ArrayList<>();
		query(searchEnv, item -> {
			result.add(item);
			return true; // keep visiting
		});
		return result;
	}

	/**
	 * Visits all items whose envelopes intersect the supplied {@code searchEnv},
	 * invoking the provided {@code visitor} for each matching item.
	 *
	 * <p>
	 * Visiting order follows the internal node traversal order. The visitor may
	 * return {@code false} to request early termination of the traversal; in that
	 * case no further items will be visited.
	 * </p>
	 *
	 * <p>
	 * This method will build the index if it has not already been built.
	 * </p>
	 *
	 * @param searchEnv the envelope to search
	 * @param visitor   a callback invoked for each matching item; returning
	 *                  {@code false} signals that the traversal should stop
	 *                  immediately
	 */
	public void query(Envelope searchEnv, ItemVisitor<T> visitor) {
		build();
		if (!totalExtent.intersects(searchEnv))
			return;

		terminated = false;
		if (layerStartIndex == null) {
			queryItems(0, searchEnv, visitor);
		} else {
			queryTopLayer(searchEnv, visitor);
		}
	}

	private void queryTopLayer(Envelope searchEnv, ItemVisitor<T> visitor) {
		int layerIndex = layerStartIndex.length - 2;
		int layerSize = layerSize(layerIndex);
		// query each node in layer
		for (int i = 0; i < layerSize; i += ENV_SIZE) {
			if (terminated)
				break;
			queryNode(layerIndex, i, searchEnv, visitor);
		}
	}

	private void queryNode(int layerIndex, int nodeOffset, Envelope searchEnv, ItemVisitor<T> visitor) {
		int layerStart = layerStartIndex[layerIndex];
		int nodeIndex = layerStart + nodeOffset;
		if (!intersects(nodeBounds, nodeIndex, searchEnv))
			return;
		if (terminated)
			return;

		if (layerIndex == 0) {
			int childNodesOffset = nodeOffset / ENV_SIZE * nodeCapacity;
			queryItems(childNodesOffset, searchEnv, visitor);
		} else {
			int childNodesOffset = nodeOffset * nodeCapacity;
			queryNodeChildren(layerIndex - 1, childNodesOffset, searchEnv, visitor);
		}
	}

	private static boolean intersects(double[] bounds, int nodeIndex, Envelope env) {
		boolean isBeyond = (env.getMaxX() < bounds[nodeIndex]) || (env.getMaxY() < bounds[nodeIndex + 1]) || (env.getMinX() > bounds[nodeIndex + 2])
				|| (env.getMinY() > bounds[nodeIndex + 3]);
		return !isBeyond;
	}

	private void queryNodeChildren(int layerIndex, int blockOffset, Envelope searchEnv, ItemVisitor<T> visitor) {
		int layerStart = layerStartIndex[layerIndex];
		int layerEnd = layerStartIndex[layerIndex + 1];
		for (int i = 0; i < nodeCapacity; i++) {
			if (terminated)
				break;
			int nodeOffset = blockOffset + ENV_SIZE * i;
			// don't query past layer end
			if (layerStart + nodeOffset >= layerEnd)
				break;

			queryNode(layerIndex, nodeOffset, searchEnv, visitor);
		}
	}

	@SuppressWarnings("unchecked")
	private void queryItems(int blockStart, Envelope searchEnv, ItemVisitor<T> visitor) {
		for (int i = 0; i < nodeCapacity; i++) {
			if (terminated)
				return;

			int itemIndex = blockStart + i;
			// don't query past end of items
			if (itemIndex >= numItems)
				break;
			if (intersects(itemBounds, itemIndex * ENV_SIZE, searchEnv)) {
				T value = (T) itemValues[itemIndex];
				if (!visitor.visitItem(value)) {
					terminated = true;
					return;
				}
			}
		}
	}

	private int layerSize(int layerIndex) {
		int layerStart = layerStartIndex[layerIndex];
		int layerEnd = layerStartIndex[layerIndex + 1];
		return layerEnd - layerStart;
	}

	/**
	 * Builds the internal spatial index structures from the items that have been
	 * inserted.
	 *
	 * <p>
	 * This method is idempotent and thread-safe: if multiple threads call
	 * {@code build()} concurrently the index will be prepared only once. Building
	 * performs sorting (Hilbert encoding) and allocates internal arrays; it can be
	 * relatively expensive for large datasets. If the tree contains no more than
	 * {@code nodeCapacity} items, no internal node structures are created and
	 * queries will operate by linear scan.
	 * </p>
	 */
	public void build() {
		// skip if already built
		if (!isBuilt) {
			synchronized (this) {
				if (!isBuilt) {
					prepareIndex();
					prepareItems();
					this.isBuilt = true;
				}
			}
		}
	}

	private void prepareIndex() {
		// don't need to build an empty or very small tree
		if (itemsToLoad.size() <= nodeCapacity)
			return;

		sortItems();

		layerStartIndex = computeLayerIndices(numItems, nodeCapacity);
		// allocate storage
		int nodeCount = layerStartIndex[layerStartIndex.length - 1] / 4;
		nodeBounds = createBoundsArray(nodeCount);

		// compute tree nodes
		computeLeafNodes(layerStartIndex[1]);
		for (int i = 1; i < layerStartIndex.length - 1; i++) {
			computeLayerNodes(i);
		}
	}

	private void prepareItems() {
		// copy item contents out to arrays for querying
		int boundsIndex = 0;
		int valueIndex = 0;
		itemBounds = new double[itemsToLoad.size() * 4];
		itemValues = new Object[itemsToLoad.size()];
		for (Item<T> item : itemsToLoad) {
			Envelope envelope = item.getEnvelope();
			itemBounds[boundsIndex++] = envelope.getMinX();
			itemBounds[boundsIndex++] = envelope.getMinY();
			itemBounds[boundsIndex++] = envelope.getMaxX();
			itemBounds[boundsIndex++] = envelope.getMaxY();
			itemValues[valueIndex++] = item.getItem();
		}
		// and let GC free the original list
		itemsToLoad = null;
	}

	private static double[] createBoundsArray(int size) {
		double[] a = new double[4 * size];
		for (int i = 0; i < size; i++) {
			int index = 4 * i;
			a[index] = Double.MAX_VALUE;
			a[index + 1] = Double.MAX_VALUE;
			a[index + 2] = -Double.MAX_VALUE;
			a[index + 3] = -Double.MAX_VALUE;
		}
		return a;
	}

	private void computeLayerNodes(int layerIndex) {
		int layerStart = layerStartIndex[layerIndex];
		int childLayerStart = layerStartIndex[layerIndex - 1];
		int layerSize = layerSize(layerIndex);
		int childLayerEnd = layerStart;
		for (int i = 0; i < layerSize; i += ENV_SIZE) {
			int childStart = childLayerStart + nodeCapacity * i;
			computeNodeBounds(layerStart + i, childStart, childLayerEnd);
		}
	}

	private void computeNodeBounds(int nodeIndex, int blockStart, int nodeMaxIndex) {
		for (int i = 0; i <= nodeCapacity; i++) {
			int index = blockStart + 4 * i;
			if (index >= nodeMaxIndex)
				break;
			updateNodeBounds(nodeIndex, nodeBounds[index], nodeBounds[index + 1], nodeBounds[index + 2], nodeBounds[index + 3]);
		}
	}

	private void computeLeafNodes(int layerSize) {
		for (int i = 0; i < layerSize; i += ENV_SIZE) {
			computeLeafNodeBounds(i, nodeCapacity * i / 4);
		}
	}

	private void computeLeafNodeBounds(int nodeIndex, int blockStart) {
		for (int i = 0; i <= nodeCapacity; i++) {
			int itemIndex = blockStart + i;
			if (itemIndex >= itemsToLoad.size())
				break;
			Envelope env = itemsToLoad.get(itemIndex).getEnvelope();
			updateNodeBounds(nodeIndex, env.getMinX(), env.getMinY(), env.getMaxX(), env.getMaxY());
		}
	}

	private void updateNodeBounds(int nodeIndex, double minX, double minY, double maxX, double maxY) {
		if (minX < nodeBounds[nodeIndex])
			nodeBounds[nodeIndex] = minX;
		if (minY < nodeBounds[nodeIndex + 1])
			nodeBounds[nodeIndex + 1] = minY;
		if (maxX > nodeBounds[nodeIndex + 2])
			nodeBounds[nodeIndex + 2] = maxX;
		if (maxY > nodeBounds[nodeIndex + 3])
			nodeBounds[nodeIndex + 3] = maxY;
	}

	private static int[] computeLayerIndices(int itemSize, int nodeCapacity) {
		IntArrayList layerIndexList = new IntArrayList();
		int layerSize = itemSize;
		int index = 0;
		do {
			layerIndexList.add(index);
			layerSize = numNodesToCover(layerSize, nodeCapacity);
			index += ENV_SIZE * layerSize;
		} while (layerSize > 1);
		return layerIndexList.toArray();
	}

	/**
	 * Computes the number of blocks (nodes) required to cover a given number of
	 * children.
	 * 
	 * @param nChild
	 * @param nodeCapacity
	 * @return the number of nodes needed to cover the children
	 */
	private static int numNodesToCover(int nChild, int nodeCapacity) {
		int mult = nChild / nodeCapacity;
		int total = mult * nodeCapacity;
		if (total == nChild)
			return mult;
		return mult + 1;
	}

	/**
	 * Returns the extents (minimum bounding rectangles) of the internal index
	 * nodes.
	 *
	 * <p>
	 * The returned array contains one {@link Envelope} per internal node in the
	 * index. The ordering corresponds to the internal node ordering used by this
	 * tree (i.e. the same index offsets used by queries and nearest-neighbor
	 * traversal). Each {@link Envelope} in the returned array is a newly created
	 * object and may be freely modified by the caller.
	 * </p>
	 *
	 * @return an array of internal node envelopes, or an empty array if the index
	 *         contains no internal nodes
	 */
	public Envelope[] getBounds() {
		int numNodes = nodeBounds.length / 4;
		Envelope[] bounds = new Envelope[numNodes];
		// create from largest to smallest
		for (int i = numNodes - 1; i >= 0; i--) {
			int boundIndex = 4 * i;
			bounds[i] = new Envelope(nodeBounds[boundIndex], nodeBounds[boundIndex + 2], nodeBounds[boundIndex + 1], nodeBounds[boundIndex + 3]);
		}
		return bounds;
	}

	private void sortItems() {
		HilbertEncoder encoder = new HilbertEncoder(HILBERT_LEVEL, totalExtent);
		int[] hilbertValues = new int[itemsToLoad.size()];
		int pos = 0;
		for (Item<T> item : itemsToLoad) {
			hilbertValues[pos++] = encoder.encode(item.getEnvelope());
		}
		quickSortItemsIntoNodes(hilbertValues, 0, itemsToLoad.size() - 1);
	}

	private void quickSortItemsIntoNodes(int[] values, int lo, int hi) {
		// stop sorting when left/right pointers are within the same node
		// because queryItems just searches through them all sequentially
		if (lo / nodeCapacity < hi / nodeCapacity) {
			int pivot = hoarePartition(values, lo, hi);
			quickSortItemsIntoNodes(values, lo, pivot);
			quickSortItemsIntoNodes(values, pivot + 1, hi);
		}
	}

	private int hoarePartition(int[] values, int lo, int hi) {
		int pivot = values[(lo + hi) >> 1];
		int i = lo - 1;
		int j = hi + 1;

		while (true) {
			do
				i++;
			while (values[i] < pivot);
			do
				j--;
			while (values[j] > pivot);
			if (i >= j)
				return j;
			swapItems(values, i, j);
		}
	}

	private void swapItems(int[] values, int i, int j) {
		Item<T> tmpItemp = itemsToLoad.get(i);
		itemsToLoad.set(i, itemsToLoad.get(j));
		itemsToLoad.set(j, tmpItemp);

		int tmpValue = values[i];
		values[i] = values[j];
		values[j] = tmpValue;
	}

	/**
	 * Finds and returns the item nearest to the query point {@code q}, according to
	 * the user-supplied distance function {@code distFn}.
	 *
	 * <p>
	 * The {@link DistanceToItem} callback is used to compute the distance from the
	 * query coordinate to a stored item. The implementation performs a best-first
	 * traversal of the tree using bounding-box minimum-distance pruning; only items
	 * and nodes that can beat the current best distance are explored. If the index
	 * has no internal layers the search degenerates to a linear scan of all items.
	 * </p>
	 *
	 * <p>
	 * Notes and guarantees:
	 * <ul>
	 * <li>If the tree contains no items, {@code null} is returned.</li>
	 * <li>{@code distFn.distance} should return a non-negative scalar; correctness
	 * of pruning assumes that the function is consistent with geometric distance
	 * ordering.</li>
	 * <li>The returned object is the single item for which
	 * {@code distFn.distance(q, item)} is minimal (ties broken by traversal
	 * order).</li>
	 * </ul>
	 *
	 * @param q      the query coordinate
	 * @param distFn a callback that computes the distance between {@code q} and an
	 *               item
	 * @return the nearest item according to {@code distFn}, or {@code null} if the
	 *         tree is empty
	 */
	@SuppressWarnings("unchecked")
	public T nearestNeighbor(Coordinate q, DistanceToItem<T> distFn) {
		build();
		if (numItems == 0) {
			return null;
		}

		final double qx = q.x, qy = q.y;

		T best = null;
		double bestDist = Double.POSITIVE_INFINITY;
		double bestDistSq = Double.POSITIVE_INFINITY;

		// Degenerate case: no internal layers; scan items linearly
		if (layerStartIndex == null) {
			for (int i = 0; i < numItems; i++) {
				T value = (T) itemValues[i];
				double d = distFn.distance(q, value);
				if (d < bestDist) {
					bestDist = d;
					bestDistSq = d * d;
					best = value;
					if (bestDist == 0.0)
						return best;
				}
			}
			return best;
		}

		final int topLayer = layerStartIndex.length - 2;
		final int topLayerSize = layerSize(topLayer);
		final int nTopNodes = topLayerSize / ENV_SIZE;
		final PriorityQueue<HEntry> H = new PriorityQueue<>(Math.max(16, nTopNodes), Comparator.comparingDouble(e -> e.minDist));

		// Seed with top-layer nodes; push only if they can beat current best
		for (int off = 0; off < topLayerSize; off += ENV_SIZE) {
			double mdSq = mindistToNodeSquared(qx, qy, topLayer, off);
			if (mdSq < bestDistSq) {
				H.add(new HEntry(topLayer, off, mdSq));
			}
		}

		while (!H.isEmpty()) {
			HEntry r = H.poll();
			if (r.minDist >= bestDistSq)
				break; // cannot improve best

			if (r.layerIndex == 0) {
				// Leaf-adjacent: scan items in this node
				int blockStart = (r.nodeOffset / ENV_SIZE) * nodeCapacity;
				for (int i = 0; i < nodeCapacity; i++) {
					int itemIndex = blockStart + i;
					if (itemIndex >= numItems)
						break;
					T value = (T) itemValues[itemIndex];
					double d = distFn.distance(q, value);
					if (d < bestDist) {
						bestDist = d;
						bestDistSq = d * d;
						best = value;
						if (bestDist == 0.0)
							return best;
					}
				}
			} else {
				// Intermediate: push child MBRs that can still beat best
				int childLayer = r.layerIndex - 1;
				int childLayerStart = layerStartIndex[childLayer];
				int childLayerEnd = layerStartIndex[childLayer + 1];
				int blockStart = r.nodeOffset * nodeCapacity;
				for (int i = 0; i < nodeCapacity; i++) {
					int childNodeOffset = blockStart + ENV_SIZE * i;
					if (childLayerStart + childNodeOffset >= childLayerEnd)
						break;
					double mdSq = mindistToNodeSquared(qx, qy, childLayer, childNodeOffset);
					if (mdSq < bestDistSq) {
						H.add(new HEntry(childLayer, childNodeOffset, mdSq));
					}
				}
			}
		}
		return best;
	}

	private double mindistToNodeSquared(double qx, double qy, int layerIndex, int nodeOffset) {
		int nodeIndex = layerStartIndex[layerIndex] + nodeOffset;
		double minX = nodeBounds[nodeIndex];
		double minY = nodeBounds[nodeIndex + 1];
		double maxX = nodeBounds[nodeIndex + 2];
		double maxY = nodeBounds[nodeIndex + 3];
		return pointToRectDistanceSquared(qx, qy, minX, minY, maxX, maxY);
	}

	// Rename: squared distance, no sqrt
	private static double pointToRectDistanceSquared(double qx, double qy, double minX, double minY, double maxX, double maxY) {
		double dx = 0.0;
		if (qx < minX)
			dx = minX - qx;
		else if (qx > maxX)
			dx = qx - maxX;

		double dy = 0.0;
		if (qy < minY)
			dy = minY - qy;
		else if (qy > maxY)
			dy = qy - maxY;

		return dx * dx + dy * dy;
	}

	/**
	 * Functional interface used to compute the distance between a query coordinate
	 * and an item stored in the tree.
	 *
	 * @param <T> item type stored in the tree
	 */
	@FunctionalInterface
	public interface DistanceToItem<T> {
		/**
		 * Computes the scalar distance from the query point {@code q} to the supplied
		 * {@code item}.
		 *
		 * <p>
		 * Implementations should return a non-negative value. The tree uses these
		 * distances for ranking candidates during nearest-neighbor search; for best
		 * performance the function should be inexpensive to compute.
		 * </p>
		 *
		 * @param q    the query coordinate
		 * @param item the item whose distance to {@code q} should be computed
		 * @return the distance from {@code q} to {@code item} (must be &gt;= 0)
		 */
		double distance(Coordinate q, T item);
	}

	/**
	 * Visitor interface used by {@link #query(Envelope, ItemVisitor)} to process
	 * items that intersect a query envelope.
	 *
	 * @param <T> the type of items visited
	 */
	public interface ItemVisitor<T> {
		/**
		 * Invoked for each item whose envelope intersects the current query envelope.
		 *
		 * @param item the item being visited
		 * @return {@code true} to continue visiting additional items, {@code false} to
		 *         terminate the traversal early
		 */
		boolean visitItem(T item);
	}

	private static class Item<T> {
		private final Envelope env;
		private final T item;

		Item(Envelope env, T item) {
			this.env = env;
			this.item = item;
		}

		public Envelope getEnvelope() {
			return env;
		}

		public T getItem() {
			return item;
		}
	}

	private static final class HEntry {
		final int layerIndex;
		final int nodeOffset; // offset in doubles within the layer (multiple of 4)
		final double minDist; // actually squared distance

		HEntry(int layerIndex, int nodeOffset, double minDistSq) {
			this.layerIndex = layerIndex;
			this.nodeOffset = nodeOffset;
			this.minDist = minDistSq;
		}
	}
}
