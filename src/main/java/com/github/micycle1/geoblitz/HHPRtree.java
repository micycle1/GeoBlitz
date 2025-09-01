package com.github.micycle1.geoblitz;

import java.util.ArrayList;
import java.util.List;

import java.util.Comparator;
import java.util.PriorityQueue;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.index.hprtree.HilbertEncoder;
import org.locationtech.jts.util.IntArrayList;

public class HHPRtree<T> {

	public interface ItemVisitor<T> {
		/**
		 * @return true to continue visiting, false to terminate early
		 */
		boolean visitItem(T item);
	}

	private static final int ENV_SIZE = 4;
	private static final int HILBERT_LEVEL = 12;
	private static final int DEFAULT_NODE_CAPACITY = 16;

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

	private List<Item<T>> itemsToLoad = new ArrayList<>();

	private final int nodeCapacity;

	private int numItems = 0;

	private final Envelope totalExtent = new Envelope();

	private int[] layerStartIndex;

	private double[] nodeBounds;

	private double[] itemBounds;

	private Object[] itemValues;

	private volatile boolean isBuilt = false;

	// Used to support early termination during query
	private boolean terminated = false;

	/**
	 * Creates a new index with the default node capacity.
	 */
	public HHPRtree() {
		this(DEFAULT_NODE_CAPACITY);
	}

	/**
	 * Creates a new index with the given node capacity.
	 * 
	 * @param nodeCapacity the node capacity to use
	 */
	public HHPRtree(int nodeCapacity) {
		this.nodeCapacity = nodeCapacity;
	}

	/**
	 * Gets the number of items in the index.
	 * 
	 * @return the number of items
	 */
	public int size() {
		return numItems;
	}

	public void insert(Envelope itemEnv, T item) {
		if (isBuilt) {
			throw new IllegalStateException("Cannot insert items after tree is built.");
		}
		numItems++;
		itemsToLoad.add(new Item<>(itemEnv, item));
		totalExtent.expandToInclude(itemEnv);
	}

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
	 * Builds the index, if not already built.
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
	 * Gets the extents of the internal index nodes
	 * 
	 * @return a list of the internal node extents
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

	@FunctionalInterface
	public interface DistanceToItem<T> {
		/**
		 * Returns the distance from the query point q to the item.
		 */
		double distance(Coordinate q, T item);
	}

	@SuppressWarnings("unchecked")
	public T nearestNeighbor(Coordinate q, DistanceToItem<T> distFn) {
		build();
		if (numItems == 0)
			return null;

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

	// Update HEntry field name to reflect squared distance
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
