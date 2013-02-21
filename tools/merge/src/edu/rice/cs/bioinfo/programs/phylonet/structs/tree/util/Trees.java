/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util;

import edu.rice.cs.bioinfo.programs.phylonet.algos.lca.SchieberVishkinLCA;
import edu.rice.cs.bioinfo.programs.phylonet.structs.BitVector;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeBipartition;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import sun.java2d.SunGraphicsEnvironment;

import java.util.*;

/**
 * This class contains some basic routines that are commonly performed on trees.
 *
 * @author Derek Ruths
 */
public class Trees {

	private Trees() {}

	// static methods
	/**
	 * @param t is the tree to compute the leaf distance matrix for.
	 * @param id_lookup is the hashtable that will be populated with a lookup mapping a leaf to its index in the matrix.
	 * @param node_lookup is the array that maps a node index to the leaf itself.
	 *
	 * @return the distance matrix for all leaves in the specified tree.
	 */
	public static final double[][] getLeafDistanceMatrix(Tree t, Hashtable<TNode,Integer> id_lookup, TNode[] node_lookup) {

		// get the list of leaves
		for(TNode n : t.getNodes()) {
			if(n.isLeaf()) {
				node_lookup[id_lookup.size()] = n;
				id_lookup.put(n,id_lookup.size());
			}
		}

		// compute distances
		double[][] distances = new double[id_lookup.size()][id_lookup.size()];

		for(int i = 0; i < distances.length; i++) {
			findDistances(i, distances[i], id_lookup, node_lookup);
		}

		return distances;
	}

	protected static final void findDistances(int leaf_idx, double[] distances, Hashtable<TNode,Integer> id_lookup, TNode[] node_lookup) {
		distances[leaf_idx] = 0;
		TNode leaf = node_lookup[leaf_idx];

		recurComputeDist(leaf, leaf.getParent(), leaf.getParentDistance(), distances, id_lookup, node_lookup);
	}

	protected static final void recurComputeDist(TNode incoming, TNode curr_node, double curr_dist, double[] distances, Hashtable<TNode,Integer> id_lookup, TNode[] node_lookup) {

		// if we're a leaf, compute the distance
		if(curr_node.isLeaf()) {
			int idx = id_lookup.get(curr_node);
			distances[idx] = curr_dist;

			return;
		}

		// recur to parent if the incoming node wasn't this node's parent
		if(curr_node.getParent() != incoming && curr_node.getParent() != null) {
			recurComputeDist(curr_node, curr_node.getParent(), curr_dist + curr_node.getParentDistance(), distances, id_lookup, node_lookup);
		}

		// recur to all children
		for(TNode child : curr_node.getChildren()) {

			if(child != incoming) {
				recurComputeDist(curr_node, child, curr_dist + child.getParentDistance(), distances, id_lookup, node_lookup);
			}
		}

		return;
	}

	/**
	 * This method locates binary nodes (nodes with a parent and one child) in the tree.
	 * If the tree doesn't contain any, this method returns <code>null</code>.
	 */
	public static final TNode findBinaryNodes(Tree t) {

		for(TNode n : t.getNodes()) {
			if(n.getChildCount() == 1) {
				return n;
			}
		}

		return null;
	}

	public static final int getLeafCount(TNode n) {
		if(n.isLeaf()) {
			return 1;
		} else {
			int count = 0;

			for(TNode child : n.getChildren()) {
				count += getLeafCount(child);
			}

			return count;
		}
	}

	/**
	 * Scale the branch lengths of the tree by a specified amount.
	 */
	public static final void scaleBranchLengths(MutableTree tree, double scaling) {
		for(TMutableNode n : tree.getNodes()) {
			if(n.getParentDistance() != TMutableNode.NO_DISTANCE) {
				n.setParentDistance(n.getParentDistance() * scaling);
			}
		}

		return;
	}

	/**
	 * Apply a labeling to ensure that all nodes are labeled in the tree.  Any node
	 * in the tree specified that has either a <code>null</code> or empty string name
	 * is renamed with a unique name of the form "I#" where # is an integer chosen that
	 * makes the name unique.
	 */
	public static final void autoLabelNodes(MutableTree tree) {

		autoLabelNodes(tree, "I");

		return;
	}

	/**
	 * Apply a labeling to ensure that all nodes are labeled in the tree.  Any node
	 * in the tree specified that has either a <code>null</code> or empty string name
	 * is renamed with a unique name of the form <code>prefix</code># where # is an integer
	 * chosen that makes the name unique.
	 *
	 * @param tree is the tree to rename nodes in
	 * @param prefix is the prefix of the names that should be constructed.
	 */
	public static final void autoLabelNodes(MutableTree tree, String prefix) {
		int node_num = 0;

		for(TMutableNode n : tree.getNodes()) {
			if(n.getName() == null || n.getName().equals("")) { // TODO Make sure there aren't any other == string comparisons!
				String new_name = prefix + node_num++;

				while(tree.getNode(new_name) != null) {
					new_name = prefix + node_num++;
				}

				n.setName(new_name);
			}
		}

		return;
	}

	/**
	 * Remove all nodes with only one child in this tree.  The tree can be unrooted or rooted.
	 */
	public static final void removeBinaryNodes(MutableTree tree) {

		// walk the whole tree
		removeBinaryChildren(tree.getRoot());

		// check the root
		TMutableNode oroot = tree.getRoot();
		if(oroot.getChildCount() == 1) {
			TMutableNode n = oroot.getChildren().iterator().next();

			n.makeRoot();
			n.removeChild(oroot, false);
		}

		return;
	}


	private static final void removeBinaryChildren(TMutableNode node) {

		// copy iterator
		HashSet<TMutableNode> hs = new HashSet<TMutableNode>();
		for(TMutableNode child : node.getChildren()) {
			hs.add(child);
		}

		// remove binary children
		for(TMutableNode child : hs) {
			removeBinaryChildren(child);

			if(child.getChildCount() == 1) {
				node.removeChild(child, true);
			}
		}
	}


	/**
	 * compute if a tree is binary
	 *
	 * @param tree
	 * @return true if the tree is binary
	 * @author yy9
     */

	public static final boolean isBinary(Tree tree){
		for(TNode node: tree.getNodes()){
			if(!node.isLeaf()){
				if(node.getChildCount()>2){
					return false;
				}
			}
		}

		return true;
	}

	/**
	 * @return <code>true</code> if these trees have identical leafsets
	 * by name.
	 */
	public static final boolean leafSetsAgree(Tree t1, Tree t2) {

		// they definitely can't agree if they are different sizes
		if(t1.getLeafCount() != t2.getLeafCount()) {
			return false;
		}

		// check t1 against t2
		for(TNode node : t1.getNodes()) {
			if(node.isLeaf() && t2.getNode(node.getName()) == null) {
				return false;
			}
		}

		// check t2 against t1
		for(TNode node : t2.getNodes()) {
			if(node.isLeaf() && t1.getNode(node.getName()) == null) {
				return false;
			}
		}

		// if we got here, then they are the same.
		return true;
	}

	/**
	 * Generate a random binary tree in the object specified.  Leaves will be
	 * named '0' - num_leaves.
	 *
	 * @param t is a tree containing only a root.
	 */
	public static final void generateRandomTree(MutableTree t, int num_leaves) {

		LinkedList<TMutableNode> active_nodes = new LinkedList<TMutableNode>();

		// create all leaves
		for(int i = 0; i < num_leaves; i++) {
			active_nodes.add(t.getRoot().createChild("" + i));
		}

		// coalese the tree
		while(active_nodes.size() > 2) {
			int size = active_nodes.size();
			int pos1 = 1;
			int pos2 = 1;

			while(pos1 == pos2) {
				pos1 = (int) Math.floor(Math.random() * size);
				pos2 = (int) Math.floor(Math.random() * size);
			}

			// coalese the nodes
			TMutableNode node1 = (TMutableNode) active_nodes.get(pos1);
			TMutableNode node2 = (TMutableNode) active_nodes.get(pos2);

			active_nodes.remove(node1);
			active_nodes.remove(node2);

			TMutableNode node = t.getRoot().createChild();
			node.adoptChild(node1);
			node.adoptChild(node2);

			active_nodes.add(node);
		}

		return;
	}

	/**
	 * Compute the support for each edge in <code>support_tree</code> to be the
	 * fraction of trees in <code>trees</code> that contain that edges.
	 *
	 * @param support_tree is the tree that the supports are computed for.  The support
	 * values computed will be stored as the branch lengths of this tree.
	 *
	 * @param trees
	 */
	public static final void computeEdgeSupports(MutableTree support_tree, Iterable<Tree> trees) {

		// generate leaf assignment
		Hashtable<String,Integer> leaf_assignment = new Hashtable<String,Integer>();
		for(TNode n : support_tree.getNodes()) {
			if(n.isLeaf()) {
				leaf_assignment.put(n.getName(), leaf_assignment.size());
			}
		}

		// generate all the bipartitions
		Hashtable<BitVector,TNode> support_partitions = new Hashtable<BitVector,TNode>();
		Bipartitions.computeBipartitions(support_tree, leaf_assignment, support_partitions);

		LinkedList<Hashtable<BitVector,TNode>> tree_partitions = new LinkedList<Hashtable<BitVector,TNode>>();
		for(Tree t : trees) {
			Hashtable<BitVector,TNode> th = new Hashtable<BitVector,TNode>();
			Bipartitions.computeBipartitions(t, leaf_assignment, th);
			tree_partitions.add(th);
		}

		// compute the ratios
		for(Map.Entry<BitVector,TNode> e : support_partitions.entrySet()) {
			BitVector bvcomp = new BitVector(e.getKey());
			bvcomp.not();

			int count = 0;

			for(Hashtable<BitVector,TNode> h : tree_partitions) {
				if(h.containsKey(e.getKey()) || h.containsKey(bvcomp)) {
					count++;
				}
			}

			((TMutableNode) e.getValue()).setParentDistance(((double) count) / tree_partitions.size());
		}

		return;
	}

	/**
	 * This method computes the maximum distance across the tree specified.
	 *
	 * @param t is the tree to compute the span for.
	 *
	 * @return the maximum distance across the tree - from one taxa to another.
	 */
	public static final double computeSpan(Tree t) {
		return computeSpan(t.getRoot());
	}

	protected static final double computeSpan(TNode n) {
		if(n.isLeaf()) {
			return 0;
		} else {
			// This is the maximum span that crosses this node
			double max_crossing_span = Double.MIN_VALUE;

			// this is the maximum span of nodes BELOW this node
			double max_below_span = Double.MIN_VALUE;

			double max_hspan1 = Double.MIN_VALUE;
			double max_hspan2 = Double.MIN_VALUE;
			for(TNode child : n.getChildren()) {

				// update the crossing span
				double hspan = computeHalfSpan(child);

				if(hspan >= max_hspan1) {
					max_hspan2 = max_hspan1;
					max_hspan1 = hspan;
				} else if(hspan > max_hspan2) {
					max_hspan2 = hspan;
				}

				// update the below span
				double bspan = computeSpan(child);

				if(bspan > max_below_span) {
					max_below_span = bspan;
				}
			}

			max_crossing_span = max_hspan1 + max_hspan2;

			return Math.max(max_crossing_span, max_below_span);
		}
	}

	protected static final double computeHalfSpan(TNode n) {

		double max_child_hspan = 0;

		if(!n.isLeaf()) {
			max_child_hspan = Double.MIN_VALUE;

			for(TNode child : n.getChildren()) {
				double hspan = computeHalfSpan(child);

				if(hspan > max_child_hspan) {
					max_child_hspan = hspan;
				}
			}
		}

		return n.getParentDistance() + max_child_hspan;
	}

	/**
	 * This function prunes leaves that are in either one of the two trees, but not both. In other words, it restricts
	 * the two trees on the common set of taxa.
	 *
	 * @param tree1
	 * @param tree2
	 */
	public static void pruneLeaves(MutableTree tree1, MutableTree tree2) {
		List<String> leaves1 = new LinkedList<String>();
		List<String> leaves2 = new LinkedList<String>();

		for (TNode node : tree1.getNodes()) {
			if (node.isLeaf()) {
				leaves1.add(node.getName());
			}
		}

		for (TNode node : tree2.getNodes()) {
			if (node.isLeaf()) {
				leaves2.add(node.getName());
			}
		}

		// Find the set of leaves in common
		List<String> common = new LinkedList<String>();

		for (String name : leaves1) {
			if (leaves2.contains(name)) {
				common.add(name);
			}
		}

		// Restrict trees on those common leaves.
		tree1.constrainByLeaves(common);
		tree2.constrainByLeaves(common);
	}

	/**
	 * Build a tree from a list of clusters. All clusters are on the same set of taxa, and use the same
	 * taxon orderings.
	 *
	 * @param clusters
	 * @return
	 */
	public static Tree buildTreeFromClusters(List<STITreeCluster> clusters) {
		if (clusters == null || clusters.size() == 0) {
			System.err.println("Empty list of clusters. The function returns a null tree.");
			return null;
		}

		MutableTree tree = new STITree<Object>();

		// Create a big star tree.
		String taxa[] = clusters.get(0).getTaxa();
		for (int i = 0; i < taxa.length; i++) {
			tree.getRoot().createChild(taxa[i]);
		}

		for (STITreeCluster tc : clusters) {
			if (tc.getClusterSize() <= 1 || tc.getClusterSize() == tc.getTaxa().length) {
				continue;
			}

			// Determine which clades are to be moved.
			Set<TNode> clusterLeaves = new HashSet<TNode>();
			for (String l : tc.getClusterLeaves()) {
				TNode node = tree.getNode(l);
				clusterLeaves.add(node);
			}

			SchieberVishkinLCA lcaFinder = new SchieberVishkinLCA(tree);
			TNode lca = lcaFinder.getLCA(clusterLeaves);

			List<TNode> movedChildren = new LinkedList<TNode>();
			for (TNode child : lca.getChildren()) {
				BitSet childCluster = new BitSet(taxa.length);
				for (TNode cl : child.getLeaves()) {
					for (int i = 0; i < taxa.length; i++) {
						if (taxa[i].equals(cl.getName())) {
							childCluster.set(i);
							break;
						}
					}
				}

				BitSet temp = (BitSet) childCluster.clone();
				temp.and(tc.getCluster());
				if (temp.equals(childCluster)) {
					movedChildren.add(child);
				}
			}

			// Move those clades under a new node.
			STINode<Object> newChild = ((STINode<Object>) lca).createChild();

			while (!movedChildren.isEmpty()) {
				newChild.adoptChild((TMutableNode) movedChildren.get(0));
				movedChildren.remove(0);
			}
		}

		return tree;
	}

	/**
	 * Compute the incompatibility of two trees.
	 *
	 * @param tr1
	 * @param tr2
	 * @return [num1, num2], where num1 (num2) is the number of clusters in tr1 (tr2) that are incompatible
	 * with at least one cluster in tr2 (tr1).
	 */
	public static int[] computeNumCompatibleClusters(STITree<Object> tr1, STITree<Object> tr2) {
		String taxa[] = new String[tr1.getLeafCount()];
		int count = 0;
		for (TNode node : tr1.getNodes()) {
			if (node.isLeaf()) {
				taxa[count++] = node.getName();
			}
		}

		List<STITreeCluster> clusters1 = tr1.getClusters(taxa, false);
		List<STITreeCluster> clusters2 = tr2.getClusters(taxa, false);

		int num1 = 0;
		for (STITreeCluster c1 : clusters1) {
			boolean incomp = false;
			for (STITreeCluster c2 : clusters2) {
				if (!c1.isCompatible(c2)) {
					incomp = true;
					break;
				}
			}

			if (incomp) {
				num1++;
			}
		}

		int num2 = 0;
		for (STITreeCluster c2 : clusters2) {
			boolean incomp = false;
			for (STITreeCluster c1 : clusters1) {
				if (!c2.isCompatible(c1)) {
					incomp = true;
					break;
				}
			}

			if (incomp) {
				num2++;
			}
		}

		return new int[]{num1, num2};
	}

	/**
	 * Compute the number of incompatible bipartitions.
	 *
	 * @param tr1
	 * @param tr2
	 * @return
	 */
	public static int[] computeNumCompatibleBipartitions(STITree<Object> tr1, STITree<Object> tr2) {
		String taxa[] = new String[tr1.getLeafCount()];
		int count = 0;
		for (TNode node : tr1.getNodes()) {
			if (node.isLeaf()) {
				taxa[count++] = node.getName();
			}
		}

		List<STITreeBipartition> splits1 = tr1.getBipartitions(taxa);
		List<STITreeBipartition> splits2 = tr2.getBipartitions(taxa);

		int num1 = 0;
		for (STITreeBipartition bp1 : splits1) {
			boolean incomp = false;
			for (STITreeBipartition bp2 : splits2) {
				if (!bp1.isCompatible(bp2)) {
					incomp = true;
					break;
				}
			}

			if (incomp) {
				num1++;
			}
		}

		int num2 = 0;
		for (STITreeBipartition bp2 : splits2) {
			boolean incomp = false;
			for (STITreeBipartition bp1 : splits1) {
				if (!bp2.isCompatible(bp1)) {
					incomp = true;
					break;
				}
			}

			if (incomp) {
				num2++;
			}
		}

		return new int[]{num1, num2};
	}

	/**
	 * get all the possible clusters over given taxa
	 *
	 * @param taxa
	 * @return
	 */
	public static List<STITreeCluster> getAllClusters(String[] taxa) {
		int n = taxa.length;

		if (n <= 0) {
			System.err.println("Empty list of taxa.");
			return null;
		}
		BitSet counter = new BitSet(n);
		boolean done = false;
		List<STITreeCluster> clusters = new ArrayList<STITreeCluster>();

		while (!done) {	// Repeat until all 2^n - 1 binary strings are generated.
			int i = 0;
			while (i < n && counter.get(i)) {
				counter.clear(i);
				i++;
			}
			if (i >= n) {
				done = true;	// Already generated all binary strings.
			}
			else {
				counter.set(i, true);
				STITreeCluster tc = new STITreeCluster(taxa);
				tc.setCluster((BitSet) counter.clone());
				if(tc.getClusterSize()>1 && tc.getClusterSize()<taxa.length)
				clusters.add(tc);
			}
		}
		return clusters;
	}

	public static List<Tree> getAllBinaryResolution(Tree tr){
		List<Tree> resolvedTrees = new ArrayList<Tree>();
		//Map<TNode, List<Tree>> node2tree = new HashMap<TNode, List<Tree>>();
		resolvedTrees.add(new STITree(tr));
		for (TNode node : new PostTraversal<Object>(tr.getRoot())) {
			int childCount = node.getChildCount();
			if(childCount>2) {
				int[] childrenid = new int[childCount];
				int index = 0;
				for(TNode child: node.getChildren()){
					childrenid[index++] = child.getID();
				}
				index = 0;
				//build a basic 3 leaf unrooted tree
				STITree<Integer> subTree = new STITree<Integer>(false);
				STINode<Integer> root = subTree.getRoot();
				STINode<Integer> newnode = root.createChild();
				newnode.setData(childrenid[index++]);
				STINode<Integer> innode = root.createChild();
				newnode = innode.createChild();
				newnode.setData(childrenid[index++]);
				newnode = innode.createChild();
				newnode.setData(childrenid[index++]);
				List<Tree> resolvedsubtree = new ArrayList<Tree>();
				resolvedsubtree.add(subTree);

				for(;index<childCount; index++){
					int id = childrenid[index];
					int count = resolvedsubtree.size();
					for(int i=0; i<count; i++){
						Tree preTree = resolvedsubtree.get(i);
						for(TNode n: new PostTraversal<Object>(preTree.getRoot())){
							if(!n.isLeaf()){
								if(n.getChildCount()!=2){
									throw new RuntimeException("Not binary!");
								}
								Iterator it = n.getChildren().iterator();
								STINode<Integer> lchild = (STINode<Integer>)(it.next());
								STITree<Integer> newTree = new STITree<Integer>(preTree);
								TNode peerChild = newTree.getNode(lchild.getID());
								TNode peerParent = peerChild.getParent();
								STINode<Integer> newchild = ((STINode<Integer>)peerParent).createChild();
								newchild.adoptChild((TMutableNode)peerChild);
								newnode = newchild.createChild();
								newnode.setData(id);
								resolvedsubtree.add(newTree);

								if(!n.isRoot()){
									TMutableNode rchild = (TMutableNode)(it.next());
									newTree = new STITree(preTree);
									peerChild = newTree.getNode(rchild.getID());
									peerParent = peerChild.getParent();
									newchild = ((STINode<Integer>)peerParent).createChild();
									newchild.adoptChild((TMutableNode)peerChild);
									newnode = newchild.createChild();
									newnode.setData(id);
									resolvedsubtree.add(newTree);
								}
							}
						}
					}
					for(int i=0; i<count; i++){
						resolvedsubtree.remove(0);
					}
				}

				int count = resolvedTrees.size();
				while(count>0){
					Tree unresolvedtree = resolvedTrees.get(0);
					for(Tree subtree: resolvedsubtree){
						for(Tree rootedsubtree: subtree.getAllRootingTrees()){
							STITree<Integer> newtree = new STITree<Integer>(unresolvedtree);
							TNode changingNode = newtree.getNode(node.getID());
							List<STINode> children = ((STINode)changingNode).removeAllChildren();
							for(TNode replacingNode : new PostTraversal<Object>(rootedsubtree.getRoot())){
								if(replacingNode.isLeaf()){
									int id = ((STINode<Integer>)replacingNode).getData();
									for(TNode peerNode: children){
										if(peerNode.getID()==id){
											if(peerNode.isLeaf()){
												((STINode)replacingNode).createChild(peerNode.getName());
											}
											else{
												for(TNode child: peerNode.getChildren()){
													((STINode)replacingNode).createChild(child);
												}
											}
											((STINode)replacingNode).setName("");
											break;
										}
									}
								}
							}

							for(TNode child: rootedsubtree.getRoot().getChildren()){
								//System.out.println(child.getName());
								((STINode)changingNode).createChild(child);
							}
							Trees.removeBinaryNodes((MutableTree)newtree);
							resolvedTrees.add(newtree);
						}
					}
					resolvedTrees.remove(0);
					count--;
				}
			}
		}
		for(Tree tree: resolvedTrees){
			for(TNode node: tree.getNodes()){
				((STINode)node).setData(null);
			}
		}
		return resolvedTrees;
	}


	/**
	 * Return all possible binary trees over a given taxa list
	 */
	public static List<Tree> generateAllBinaryTrees(String[] leaves){
		List<Tree> alltrees = new ArrayList<Tree>();

		int index = 0;
		//build a basic 3 leaf unrooted tree
		STITree threeLeafTree = new STITree(false);
		STINode root = threeLeafTree.getRoot();
		root.createChild(leaves[index++]);
		STINode innode = root.createChild();
		innode.createChild(leaves[index++]);
		innode.createChild(leaves[index++]);
		alltrees.add(threeLeafTree);

		for(;index<leaves.length; index++){
			String leaf = leaves[index];
			List<Tree> temp = new ArrayList<Tree>();
			temp.addAll(alltrees);
			alltrees.clear();
			for(Tree preTree: temp){
				for(TNode n: new PostTraversal<Object>(preTree.getRoot())){
					if(!n.isLeaf()){
						Iterator it = n.getChildren().iterator();
						STINode lchild = (STINode)(it.next());
						STITree newTree = new STITree(preTree);
						TNode peerChild = newTree.getNode(lchild.getID());
						TNode peerParent = peerChild.getParent();
						STINode newchild = ((STINode<Integer>)peerParent).createChild();
						newchild.adoptChild((TMutableNode)peerChild);
						newchild.createChild(leaf);
						alltrees.add(newTree);
						if(!n.isRoot()){
							STINode rchild = (STINode)(it.next());
							STITree newTree2 = new STITree(preTree);
							TNode peerChild2 = newTree2.getNode(rchild.getID());
							TNode peerParent2 = peerChild2.getParent();
							STINode newchild2 = ((STINode<Integer>)peerParent2).createChild();
							newchild2.adoptChild((TMutableNode)peerChild2);
							newchild2.createChild(leaf);
							alltrees.add(newTree2);
						}
					}
				}
			}
			temp.clear();
		}

		List<Tree> temp = new ArrayList<Tree>();
		temp.addAll(alltrees);
		alltrees.clear();
		for(Tree unrootedTree: temp){
			alltrees.addAll(((STITree)unrootedTree).getAllRootingTrees());
		}
		temp.clear();
		return alltrees;
	}

	public static String checkMapping(List<Tree> trees, Map<String,String> taxonMap){
		String error = null;
		for(Tree tr: trees){
			for(String leaf: tr.getLeaves()){
				if(!taxonMap.containsKey(leaf)){
					error = leaf;
					return error;
				}
			}
		}
		return null;
	}

	/**
	 * Extract a tree with a given bootstrap threshold
	 */
	public static int handleBootStrapInTree(Tree tr, double threshold){
		List<TNode> nodestomodify = new ArrayList<TNode>();

		for(TNode node: new PostTraversal<Object>(tr.getRoot())){
            if(node.isLeaf() || node.isRoot()){
                continue;
            }
			Double bootstrap = ((STINode<Double>)node).getData();
			if(bootstrap==null){
				return -1;
			}
			if(bootstrap < threshold){
				nodestomodify.add(node);
			}
		}

		for(TNode node: nodestomodify){
			TNode parent = node.getParent();
			double distance = node.getParentDistance();
			for(TNode child: node.getChildren()){
				child.setParentDistance(child.getParentDistance()+distance);
			}
			((STINode)parent).removeChild((TMutableNode)node, true);
		}
		return 0;
	}


    public static boolean haveSameRootedTopology(Tree t1, Tree t2){
        if(!leafSetsAgree(t1,t2)){
            return false;
        }
        String[] taxa = t1.getLeaves();
        List<STITreeCluster> clusters1 = t1.getClusters(taxa, false);
        List<STITreeCluster> clusters2 = t2.getClusters(taxa, false);
        int fp= 0;
        for(STITreeCluster cl1: clusters1){
            if(!clusters2.contains(cl1)){
                fp++;
            }
        }
        int fn = 0;
        for(STITreeCluster cl2: clusters2){
            if(!clusters1.contains(cl2)){
                fn++;
            }
        }
        return fp==0 && fn==0;
    }

}

