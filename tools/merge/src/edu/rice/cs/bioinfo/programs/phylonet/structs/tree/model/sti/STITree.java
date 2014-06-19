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

package edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti;

import java.io.CharArrayWriter;
import java.io.IOException;
import java.io.StringReader;
import java.util.*;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickWriter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.PostTraversal;


/**
 * This class provides an implementation of a mutable tree whose nodes contain a specific type of data.
 * 
 * @author Derek Ruths
 *
 * @param <D> is the type of data each node contains
 * 
 * TODO: Fix package scope sharing of fields between tree and node
 */
public class STITree<D extends Object> implements MutableTree {

	// static fields
	/**
	 * This object is used for writing the tree.
	 */
	protected static CharArrayWriter SWRITER;
	
	/**
	 * This object is used with {@link SWRITER} to write string representations of the tree.
	 */
	protected static NewickWriter NWRITER;
	
	protected static final int UNCOUNTED = -1;
	public static final String NO_NAME = "";
	
	static {
		// initialize the objects used for writing trees
		SWRITER = new CharArrayWriter();
		NWRITER = new NewickWriter(SWRITER);
	}
	
	// fields
	protected String _name = Tree.NO_NAME;
	
	protected Hashtable<Integer,STINode<D>> _nodes = new Hashtable<Integer,STINode<D>>();
	protected HashSet<STINode<D>> _node_set = new HashSet<STINode<D>>();
	protected Hashtable<String,STINode<D>> _name2node = new Hashtable<String,STINode<D>>();
	protected STINode<D> _root;
	protected boolean _is_rooted;
	
	protected int _next_node_id = 0;
	
	protected int _leaf_count = UNCOUNTED;
	
	// constructor
	/**
	 * Create a rooted tree with only one node.
	 */
	public STITree() {
		this(NO_NAME, true);
	}
	
	/**
	 * Construct a tree from the newick string.
       */
	public STITree(String newick_string) throws IOException, ParseException {
		NewickReader nr = new NewickReader(new StringReader(newick_string));
		
		_is_rooted = nr.isNextRooted();
		
		// create a root
		_root = new STINode<D>(this,_next_node_id++, NO_NAME, null,null);

		
		nr.readTree(this);
	}
	
	/**
	 * Create a tree with a specified rootedness.
	 * 
	 * @param rooted indicates whether the tree is rooted or not.
	 */
	public STITree(boolean rooted) {
		this(NO_NAME, rooted);
	}
	
	/**
	 * Creates a new tree.
	 * 
	 * @param node_name is the name of the initial node created.
	 * @param rooted is <code>true</code> if the tree is rooted.  If so, the first node created is the root of the tree.
	 */
	public STITree(String node_name, boolean rooted) {
		_is_rooted = rooted;
		
		_root = new STINode<D>(this,_next_node_id++, node_name, null,null);
	}
	
	/**
	 * Make a copy of the specified tree.
	 * 
	 * @param t is the tree to copy.
	 */
	public STITree(Tree t) {
		this(t.getRoot(), t.isRooted());
	}
	
	/**
	 * Creates a new rooted tree with <code>root</code> as the root of the tree.
	 */
	public STITree(TNode root) {
		this(root, true);
	}
	
	/**
	 * Creates a new tree with <code>root</code> as the root of the tree.  All information is copied and preserved.
	 * ID values are preserved.
	 * 
	 * @param root specifies where the root should start.
	 */
	public STITree(TNode root, boolean rooted) {
		_is_rooted = rooted;
		
		_root = new STINode<D>(this, root.getID(), root.getName(), null, ((STINode<D>)root).getData());
		
		copyNode(root, _root);
	}
	
	// methods
	public int getHeight() {
		return _root.getHeight();
	}
	
	public String getName() {
		return _name;
	}
	
	public void setName(String name) {
		_name = name;
	}
	
	protected void copyNode(TNode src, STINode<D> cpy) {
		// copy the parent distance
		cpy.setParentDistance(src.getParentDistance());
		
		// copy the children
		for(TNode src_child : src.getChildren()) {
			STINode<D> cpy_child = new STINode<D>(this,src_child.getID(), src_child.getName(), null,((STINode<D>)src_child).getData());
			
			if(cpy_child.getID() >= _next_node_id) {
				_next_node_id = cpy_child.getID() + 1;
			}
			
			cpy.adoptChild(cpy_child);
			
			copyNode(src_child, cpy_child);
		}
	}
	
	public boolean isRooted() {
		return _is_rooted;
	}
	
	public STINode<D> getRoot() {
		return _root;
	}

	public int getNodeCount() {
		return _nodes.size();
	}

	public int getLeafCount() {
		
		if(_root == null) {
			return 0;
		}
		
		if(_leaf_count == UNCOUNTED) {
			_leaf_count = _root.getLeafCount();
		}
		
		return _leaf_count;
	}
	
	public String[] getLeaves(){
		String[] leaves = new String[getLeafCount()];
		int index = 0;
		for(TNode leaf: _root.getLeaves()){
			leaves[index++]=leaf.getName();
		}
		return leaves;
	}

	public STINode<D> getNode(int id) {
		return _nodes.get(new Integer(id));
	}
	
	public Iterable<STINode<D>> getNodes() {
		return new Iterable<STINode<D>>() {
			public Iterator<STINode<D>> iterator() {
				return _nodes.values().iterator();
			}
		};
	}
	
	protected void removeNodeRecord(STINode<D> node) {
		_nodes.remove(node.getID());
		_node_set.remove(node);
		_name2node.remove(node.getName());

		/*
		 * All recursion is handled by the removeNode method in STINode
		// remove all the node's children
		for(STINode<D> child : node.getChildren()) {
			removeNode(child);
		}
		*/
		
		if (_nodes.size() != _node_set.size()) {
			System.err.println("removeNode: Inconsistent _nodes and _node_set");
		}
		
		return;
	}
	
	public STINode<D> getNode(String name) {
		return _name2node.get(name);
	}
	
	public void removeNode(String name){
		getNode(name).removeNode();
	}
	
	/**
	 * This method prunes the tree so that it only contains leaves whose names
	 * are in the set <code>leaf_names</code>.  This method also removes unnecessary
	 * internal nodes from the tree. 
	 * 
	 * If this method constrains the tree such that no nodes can exist (i.e. the tree consists of
	 * one node that will be removed), the root will be set to <code>null</code> and the tree
	 * will be empty, and therefore, useless.
	 * 
	 * @param leaf_names specifies the group of leaves that should remain in the 
	 * tree once it is pruned.
	 */
	public void constrainByLeaves(Iterable<String> leaf_names) {
		Set<String> elim_leaves = new HashSet<String>(_name2node.keySet());
		
		for(String name : leaf_names) {
			elim_leaves.remove(name);
		}
		
		for(String node_name : elim_leaves) {
			STINode<D> node = _name2node.get(node_name);
			
			// skip named nodes that aren't leaves
			if(node == null || !node.isLeaf()) {
				continue;
			}
			
			STINode<D> parent = node._parent;
			
			// if this node is the root
			if(parent == null) {
				// destroy the node
				node.removeSelf();
				
				// set the root to null
				_root = null;
				
				// this tree is now useless
				return;
			}
			
			parent.removeChild(node, false);
		
			// prune the parent if it now has only one child
			if(parent._children.size() == 1) {
				
				if(parent != _root) {
					parent._parent.removeChild(parent, true);
				} else {
					STINode<D> nroot = parent._children.iterator().next();
					nroot._parent = null;
					
					_root._children.clear();
					_root.removeNode();
					
					_root = nroot;
				}
			}
		}
		
		return;
	}

	public boolean isEmpty() {
		return (getLeafCount() == 0);
	}

	public STINode<D> createRoot() {
		
		if(_root != null) {
			throw new RuntimeException("createRoot called on non-empty tree");
		}
		
		_root = new STINode<D>(this,_next_node_id++, NO_NAME, null,null);
		
		return _root;
	}
	
	public void setRooted(boolean rooted) {
		_is_rooted = rooted;
	}
	
	public List<Tree> getAllRootingTrees(){
		List<Tree> rerootingTrees = new ArrayList<Tree>();
		//int num_trees = 2*getLeafCount()-3;
		//int index=1;
		for(TNode node: postTraverse()){
			if(!node.isRoot() && !(node.getParent().isRoot())){
				Tree modified = new STITree<D>(this);
				/*for(TNode n:modified.postTraverse()){
					((STINode)n).setParentDistance(TMutableNode.NO_DISTANCE);
				}*/
				modified.rerootTreeAtEdge(node.getID());
				/*if(!Trees.isBinary(modified))
					System.out.println("Error");*/
				rerootingTrees.add(modified);
			}
		}
		rerootingTrees.add(new STITree<D>(this));
		return rerootingTrees;
	}
	
	public void rerootTreeAtEdge(int nodeID){
		rerootTreeAtEdge(this.getNode(nodeID));
	}
	
	public void rerootTreeAtEdge(String nodeName){
		rerootTreeAtEdge(this.getNode(nodeName));
	}
	
	public void rerootTreeAtEdge(TNode node){
		doRerooting(node.getParent());
		//_root.removeAllChildren()
		List<TNode> siblinglist = ((STINode)node).getSiblings();
		STINode newnode = _root.createChild();
		for(Object o: siblinglist){
			newnode.adoptChild((TMutableNode)o);
		}
	}
	
	public void rerootTreeAtNode(TNode node){
		if(!_node_set.contains(node)){
			throw new RuntimeException("node " + node + " is not in the tree "+ this.toNewick());
		}
		if(node.isRoot()){
			return;
		}
		if(node.isLeaf()){
			rerootTreeAtEdge(node);
		}
		else{
			doRerooting(node);
		}
	}


    private void doRerooting(TNode node){
        STINode parent = (STINode)(node.getParent());
        if(parent == null){
            return;
        }
        if(!parent.isRoot()){
            doRerooting(node.getParent());
        }
        parent._children.remove(node);
        if(parent.getChildCount()==1){
            ((STINode)node).adoptChild(parent);
            ((STINode)node).removeChild(parent, true);
        }
        else{
            ((STINode)node).adoptChild(parent);
        }

        _root=((STINode)node);
        _root._parent = null;
    }
	
	public STINode<D> selectRandomNode(boolean include_leaves, boolean include_root) {
		
		int idx = (int) Math.floor(Math.random() * _node_set.size());
		
		Iterator<STINode<D>> it = _node_set.iterator();
		
		for(int i = 0; i < idx; i++) {
			it.next();
		}
		
		STINode<D> n = it.next();
		
		if(!include_leaves && n.isLeaf()) {
			return selectRandomNode(include_leaves, include_root);
		} else if(!include_root && n.isRoot()) {
			return selectRandomNode(include_leaves, include_root);
		} else {
			return n;
		}
	}

	public String toNewick() {
		synchronized(SWRITER) {
			SWRITER.reset();
			NWRITER.writeTree(this, false);
		
			return SWRITER.toString();
		}
	}
	
	public String toNewickWD() {
		synchronized(SWRITER) {
			SWRITER.reset();
			NWRITER.writeTreeWD(this, false);
		
			return SWRITER.toString();
		}
	}
	
	public String toString() { return toNewick(); }
	
	public String toStringWD() { return toNewickWD(); }

	public String toString(int format) {
		switch(format) {
		case Tree.NEWICK_FORMAT:
			return toNewick();
		default:
			throw new RuntimeException("Unknown format " + format);
		}
	}
	
	
	/** 
	 * Compute all clusters in the tree.
	 * 
	 * @param leaves A list of leaves in the tree. If it is null, then the function will create a list of 
	 * leaves each time this function is called.
	 * 
	 * @return A list clusters induced by the tree.
	 */
	public List<STITreeCluster> getClusters(String leaves[], boolean gen) {
		PostTraversal<D> traversal = new PostTraversal<D>(_root);
		List<STITreeCluster> clusters = new LinkedList<STITreeCluster>();
		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		
		if (leaves == null) {
			int i = 0;
			leaves = new String[getLeafCount()];
			for (TNode node : getNodes()) {
				if (node.isLeaf()) {
					leaves[i++] = node.getName();
				}
			}
		}
		
		for (TNode node : traversal) {
			BitSet bs = new BitSet();
			
			if (node.isLeaf()) {
				for (int i = 0; i < leaves.length; i++) {
					if (node.getName().equals(leaves[i])) {
						bs.set(i);
						break;
					}
				}
				
				map.put(node, bs);
			}
			else {
				int childCount = node.getChildCount();
				BitSet[] childbslist = new BitSet[childCount];
				int index = 0;
				for (TNode child : node.getChildren()) {
					BitSet childCluster = map.get(child);
					bs.or(childCluster);
					childbslist[index++] = childCluster;
				}
				if(childCount>2 && gen){
					BitSet counter = new BitSet(childCount);
					boolean done = false;
					
					while (!done) {	// Repeat until all 2^n - 1 binary strings are generated.
						int i = 0;
						while (i < childCount && counter.get(i)) {
							counter.clear(i);
							i++;
						}
						if (i >= childCount) {
							done = true;	// Already generated all binary strings.
						}
						else {
							counter.set(i, true);
							if(counter.cardinality() > 1 && counter.cardinality() < childCount) {
								STITreeCluster tc = new STITreeCluster(leaves);
								BitSet tcbs = new BitSet(leaves.length);
								for (int j = counter.nextSetBit(0); j >= 0; j = counter.nextSetBit(j+1)) {
								     tcbs.or(childbslist[j]);
								 }			
								tc.setCluster(tcbs);
								clusters.add(tc);
							}
						}
					}
				}
				map.put(node, bs);
			}
			
			if (bs.cardinality() > 1 && bs.cardinality() < leaves.length) {
				STITreeCluster tc = new STITreeCluster(leaves);
				tc.setCluster(bs);
				clusters.add(tc);
			}
		}
		
		return clusters;
	}
	
	/*
	public List<STITreeCluster> getClusters(String leaves[]) {
		PostTraversal<D> traversal = new PostTraversal<D>(_root);
		List<STITreeCluster> clusters = new LinkedList<STITreeCluster>();
		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		
		if (leaves == null) {
			int i = 0;
			leaves = new String[getLeafCount()];
			for (TNode node : getNodes()) {
				if (node.isLeaf()) {
					leaves[i++] = node.getName();
				}
			}
		}
		
		for (TNode node : traversal) {
			BitSet bs = new BitSet();
			
			if (node.isLeaf()) {
				for (int i = 0; i < leaves.length; i++) {
					if (node.getName().equals(leaves[i])) {
						bs.set(i);
						break;
					}
				}
				
				map.put(node, bs);
			}
			else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = map.get(child);
					bs.or(childCluster);
				}
				
				map.put(node, bs);
			}
			
			if (bs.cardinality() > 1 && bs.cardinality() < leaves.length) {
				STITreeCluster tc = new STITreeCluster(leaves);
				tc.setCluster(bs);
				clusters.add(tc);
			}
		}
		
		return clusters;
	}
	*/
	
	public List<STITreeClusterWD> getClustersWD(String leaves[], boolean gen) {
		PostTraversal<D> traversal = new PostTraversal<D>(_root);
		List<STITreeClusterWD> clusters = new LinkedList<STITreeClusterWD>();
		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		
		if (leaves == null) {
			int i = 0;
			leaves = new String[getLeafCount()];
			for (TNode node : getNodes()) {
				if (node.isLeaf()) {
					leaves[i++] = node.getName();
				}
			}
		}
		
		for (TNode node : traversal) {
			BitSet bs = new BitSet();
			
			if (node.isLeaf()) {
				for (int i = 0; i < leaves.length; i++) {
					if (node.getName().equals(leaves[i])) {
						bs.set(i);
						break;
					}
				}
				
				map.put(node, bs);
			}
			else {
				int childCount = node.getChildCount();
				BitSet[] childbslist = new BitSet[childCount];
				int index = 0;
				for (TNode child : node.getChildren()) {
					BitSet childCluster = map.get(child);
					bs.or(childCluster);
					childbslist[index++] = childCluster;
				}
				if(childCount>2 && gen){
					BitSet counter = new BitSet(childCount);
					boolean done = false;
					
					while (!done) {	// Repeat until all 2^n - 1 binary strings are generated.
						int i = 0;
						while (i < childCount && counter.get(i)) {
							counter.clear(i);
							i++;
						}
						if (i >= childCount) {
							done = true;	// Already generated all binary strings.
						}
						else {
							counter.set(i, true);
							if(counter.cardinality() > 1 && counter.cardinality() < childCount) {
								STITreeClusterWD tc = new STITreeClusterWD(leaves);
								BitSet tcbs = new BitSet(leaves.length);
								for (int j = counter.nextSetBit(0); j >= 0; j = counter.nextSetBit(j+1)) {
								     tcbs.or(childbslist[j]);
								 }			
								tc.setCluster(tcbs);
								tc.setData(((STINode)node).getData());
								clusters.add(tc);
							}
						}
					}
				}
				map.put(node, bs);
			}
			
			if (bs.cardinality() > 1 && bs.cardinality() < leaves.length) {
				STITreeClusterWD tc = new STITreeClusterWD(leaves);
				tc.setCluster(bs);
				tc.setData(((STINode)node).getData());
				clusters.add(tc);
			}
		}
		
		return clusters;
	}
		
	/**
	 * Compute all bipartitions in the tree.
	 * 
	 * @param leaves: A list of leaves in the tree. If this is null, then the function will create a list of 
	 * leaves each time this function is called. 
	 * @return A list of bipartitions in the tree.
	 * @author cvthan
	 * 
	 * IMPORTANT: If you want to compare bipartitions of different trees (that have the same set of leaves),
	 * then the leaves passed to this function should be keep identical. The comparison functions in 
	 * STITreeBipartition assume the leaf assignment in bipartitions are the same.
	 */
	public List<STITreeBipartition> getBipartitions(String leaves[]) {
		PostTraversal<D> traversal = new PostTraversal<D>(_root);
		List<STITreeBipartition> bipartitions = new LinkedList<STITreeBipartition>();
		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		
		if (leaves == null) {
			int i = 0;
			
			leaves = new String[getLeafCount()];
			for (TNode node : getNodes()) {
				if (node.isLeaf()) {
					leaves[i] = node.getName();
					i++;
				}
			}
		}
		
		for (TNode node : traversal) {
			BitSet bs = new BitSet();
			
			if (node.isLeaf()) {	// Create a bipartition with a single leaf.
				// Find the index of this leaf.
				int i = 0;
				
				for (i = 0; i < leaves.length; i++) {
					if (node.getName().equals(leaves[i])) {
						break;
					}
				}
				
				// The leaf must always be found.
				assert(i < leaves.length);
				
				bs.set(i);				
				map.put(node, bs);
			}
			else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = map.get(child);
					bs.or(childCluster);
				}
				
				map.put(node, bs);
			}
			
			// Record only nontrivial, non-duplicate bipartitions.			
			if (bs.cardinality() > 1 && bs.cardinality() < leaves.length - 1) {			
				boolean duplicate = false;
				
				STITreeBipartition tb = new STITreeBipartition(leaves);
				tb.setBipartition(bs, null);
				
				for (int i = 0; i < bipartitions.size(); i++) {
					if (tb.isEqual(bipartitions.get(i))) {
						duplicate = true;
						break;
					}
				}
				
				if (!duplicate) {
					bipartitions.add(tb);
				}
			}
		}
		return bipartitions;		
	}
	
	public List<STITreeCluster> getBipartitionClusters(String leaves[], boolean gen){
		List<STITreeCluster> biClusters = new LinkedList<STITreeCluster>();
		Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
		int[] containingleaf = null;
		if (leaves == null) {
			int i = 0;
			
			leaves = new String[getLeafCount()];
			for (TNode node : getNodes()) {
				if (node.isLeaf()) {
					leaves[i] = node.getName();
					i++;
				}
			}
			containingleaf = new int[leaves.length];
			for(i=0;i<containingleaf.length;i++){
				containingleaf[i] = 1;
			}
		}
		else{
			containingleaf = new int[leaves.length];		
			for(int i=0;i<containingleaf.length;i++){
				containingleaf[i] = 0;
			}
			for (TNode node : this.postTraverse()) {
				if(node.isLeaf()){
					String name = node.getName();
					for(int i=0; i<leaves.length; i++){
						if(leaves[i].equals(name)){
							containingleaf[i] = 1;
							break;
						}
					}
				}
			}
		}
		//PostTraversal<D> traversal = new PostTraversal<D>(_root);
		//System.out.println(leaves.length);
		for (TNode node : this.postTraverse()) {
			BitSet bs = new BitSet(leaves.length);
			if (node.isLeaf()) {
				// Find the index of this leaf.
				int i = 0;
				
				for (i = 0; i < leaves.length; i++) {
					if (node.getName().equals(leaves[i])) {
						break;
					}
				}
				
				// The leaf must always be found.
				assert(i < leaves.length);
				
				bs.set(i);				
				map.put(node, bs);
			}
			else {
				int childCount = node.getChildCount();
				BitSet[] childbslist = new BitSet[childCount];
				int index = 0;
				for (TNode child : node.getChildren()) {
					BitSet childCluster = map.get(child);
					bs.or(childCluster);
					childbslist[index++] = childCluster;
				}
				//System.out.println(childCount);
				if(childCount>2 && gen){
					BitSet counter = new BitSet(childCount);
					boolean done = false;
					while (!done) {	// Repeat until all 2^n - 1 binary strings are generated.
						int i = 0;
						while (i < childCount && counter.get(i)) {
							counter.clear(i);
							i++;
						}
						if (i >= childCount) {
							done = true;	// Already generated all binary strings.
						}
						else {
							counter.set(i, true);
							//System.out.println(counter);
							if(counter.cardinality() > 1 && counter.cardinality() < childCount) {
								BitSet tcbs = new BitSet(leaves.length);
								for (int j=counter.nextSetBit(0); j>=0; j=counter.nextSetBit(j+1)) { 
									tcbs.or(childbslist[j]);
								 }			
								for(int j=0;j<2;j++){				
									if(tcbs.cardinality()<leaves.length && tcbs.cardinality()>0){
										STITreeCluster tb = new STITreeCluster(leaves);
										tb.setCluster((BitSet)tcbs.clone());
										if(!biClusters.contains(tb)){
											biClusters.add(tb);
										}
										
									}
									for(int k=0;k<leaves.length;k++){
										if(tcbs.get(k)){
											tcbs.set(k, false);
										}
										else{
											tcbs.set(k, true);
										}
									}
								}
							}
							
						}
					}
				}
				map.put(node, bs);
			}
			
			for(int i=0;i<2;i++){				
				if(bs.cardinality()<leaves.length && bs.cardinality()>0){
					STITreeCluster tb = new STITreeCluster(leaves);
					tb.setCluster((BitSet)bs.clone());
					if(!biClusters.contains(tb)){
						biClusters.add(tb);
					}
				}
				for(int j=0; j<leaves.length; j++){
					if(containingleaf[j]==1){
						bs.flip(j);
					}
				}
			}
		}
		
		return biClusters;				
	}
	
	
	public double gsi(String leaves[]){
		int groupNumber = leaves.length-1;
		if(groupNumber < 1){
			return -1;
		}
		List<String> taxa = Arrays.asList(this.getLeaves());
		BitSet gbs = new BitSet(taxa.size());
		for(String leaf: leaves){
			int index = taxa.indexOf(leaf);
			if(index == -1){
				return -1;
			}
			else{
				gbs.set(index);
			}
		}
		
		Map<TNode, BitSet> bsMap = new HashMap<TNode, BitSet>();
		boolean found = false;
		int total_degree = 0;
		int group_degree = 0;
		for (TNode node : this.postTraverse()) {
			if (node.isLeaf()) {
				int index = taxa.indexOf(node.getName());
				BitSet bs = new BitSet();
				bs.set(index);
				bsMap.put(node, bs);
			}
			else {
				BitSet bs = new BitSet();
				int degree = 0;
				for (TNode child : node.getChildren()) {
					BitSet v = bsMap.get(child);
					bs.or(v);
					degree++;
				}
				degree--;
				total_degree += degree;
				bsMap.put(node, bs);
				
				if(!found){
					BitSet temp = (BitSet) gbs.clone();
					temp.and(bs);
					found = temp.equals(gbs);
					if(!temp.isEmpty()){
						group_degree += degree; 
					}
				}				
			}
		}
		
		double observedgs = groupNumber*(1.0)/group_degree;
		double mings = groupNumber*(1.0)/total_degree;
		double gsi = (observedgs - mings)/(1-mings);
		return gsi;
	}
	
	/** 
	 * Traverse the tree in the post order.
	 * 
	 * @return: A list of nodes visited in the post order.
	 */
	public Iterable<TNode> postTraverse() {
		return new PostTraversal<D>(_root);
	}	
}
