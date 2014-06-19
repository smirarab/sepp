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

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/9/11
 * Time: 3:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class Collapse {

    protected static final String COLLAPSED_NAME_PREFIX = "__";

	// inner classes
	/**
	 * This class contains the description of what was collapsed by a
	 * call to {@link collapse}.  It can be used with the {@link expand}
	 * method to undo the collapse operation.
	 *
	 * @author Derek Ruths
	 */
	public static class CollapseDescriptor {
		public CollapseDescriptor() {}

		public Hashtable<String,STITree<Object>> _name2clade = new Hashtable<String,STITree<Object>>();
	}

	// collapse methods
	/**
	 * Collapse all clades that occur in both <code>t1</code> and <code>t2</code>.
	 *
	 * @return an object that can undo the collapse operation on trees that
	 * contain the collapsed clades.
	 */
	public static CollapseDescriptor collapse(MutableTree t1, MutableTree t2) {

		CollapseDescriptor cd = new CollapseDescriptor();

		runRule1(t1.getRoot(), t2, cd);
		return cd;
	}


	public static CollapseDescriptor collapse(List<Tree> trees){
		CollapseDescriptor cd = new CollapseDescriptor();
		Tree t1 = trees.get(0);
		trees.remove(0);
		runRule2((TMutableNode)t1.getRoot(), trees, cd);
		trees.add(0, t1);
		return cd;
	}




	/**
	 * Expand any collapsed clades in <code>tree</code> using the results of the
	 * collapse operation that created <code>cd</code>.
	 *
	 * @param cd is the object returned by the call to {@link collapse}.
	 * @param tree is the tree to expand the nodes of
	 */
	public static void expand(CollapseDescriptor cd, MutableTree tree) {

		boolean changed = true;

		while(changed) {
			changed = false;

			// go through all nodes in the collapse descriptor
			Iterator<Map.Entry<String,STITree<Object>>> i = cd._name2clade.entrySet().iterator();
			while(i.hasNext()) {
				Map.Entry<String,STITree<Object>> e = i.next();

				TMutableNode n = tree.getNode(e.getKey());
				// reconstitute the node
				if(n != null) {
					changed = true;
					reconstituteNode(n,e.getValue().getRoot());
				}
			}
		}
		Trees.removeBinaryNodes(tree);
	}

	private static void reconstituteNode(TMutableNode n, TNode in) {
		if(((STINode<Integer>)n).getData()== null || ((STINode<Integer>)n).getData()==0){
			n.setName(in.getName());
			n.setParentDistance(in.getParentDistance());

			// copy children
			for(TNode ichild : in.getChildren()) {
				TMutableNode child = n.createChild();
				reconstituteNode(child, ichild);
			}
		}
		else if(((STINode<Integer>)n).getData()==-1){
			TMutableNode parent = n.getParent();
			parent.removeChild(n, false);

			for(TNode ichild : in.getChildren()) {
				TMutableNode child = parent.createChild();
				reconstituteNode(child, ichild);
			}

		}
	}

	///////
	// Methods adapted from Cuong's implementation for RIATA-HGT
	/**
	 * The method collapses a subtree, rooted at peerNode1, and tr2. To refine the tree tr1 and tr2,
	 * calls this function with peerNode1 as the root of tr1.
	 *
	 * @param peerNode1
	 * @param tr2
	 *
	 * @return: true if there's actually a collapse; false if there's not.
	 */
	protected static boolean runRule1(TMutableNode peerNode1, MutableTree tr2, CollapseDescriptor cd)
	{
		// Base case: nothing to do with a single leaf.
		if (peerNode1.isLeaf())
			return false;

		// Non-trivial case. Recurse on peerNode1's childrent before collapse itself.
		boolean changed = false;

		// Collapse peerNode1's children.
		Iterator it = peerNode1.getChildren().iterator();
		TMutableNode currentNode = null;

		while (it.hasNext()) {
			TMutableNode oldNode = currentNode;
			currentNode =  (TMutableNode) it.next();

			if (oldNode != null) {
				boolean tmp = runRule1(oldNode, tr2, cd);
				changed = changed || tmp;
			}
		}
		if (currentNode != null) {
			boolean tmp = runRule1(currentNode, tr2, cd);
			changed = changed || tmp;
		}

		// Then, test if we can make the whole tr1 and tr2 as single leaves.
		STITree subtree1 = new STITree(peerNode1, true);
		if (isPendant(subtree1)) {
			TMutableNode peerNode2 = findPeerNode(tr2, subtree1);
			if (peerNode2 != null) {
				String newName = generateCollapsedName(peerNode2);

				// BEGIN ADDITION by Derek
				// store the clade that we're removing
				cd._name2clade.put(newName, new STITree<Object>(peerNode1));
				// END ADDITION

				fixNode(peerNode1, newName);	// Collapse in tr1.
				fixNode(peerNode2, newName);	// Collapse in tr2.

				changed = true;
			}
		}

		return changed;
	}

	/*
	protected static boolean runRule1(TMutableNode peerNode1, List<Tree> trees, CollapseDescriptor cd)
	{
		// Base case: nothing to do with a single leaf.
		if (peerNode1.isLeaf())
			return false;

		// Non-trivial case. Recurse on peerNode1's childrent before collapse itself.
		boolean changed = false;

		// Collapse peerNode1's children.
		Iterator it = peerNode1.getChildren().iterator();
		TMutableNode currentNode = null;

		while (it.hasNext()) {
			TMutableNode oldNode = currentNode;
			currentNode =  (TMutableNode) it.next();

			if (oldNode != null) {
				boolean tmp = runRule1(oldNode, trees, cd);
				changed = changed || tmp;
			}
		}
		if (currentNode != null) {
			boolean tmp = runRule1(currentNode, trees, cd);
			changed = changed || tmp;
		}

		// Then, test if we can make the whole tr1 and tr2 as single leaves.
		STITree subtree1 = new STITree(peerNode1, true);
		if (isPendant(subtree1)) {
			List<TMutableNode> peerNodelist = findPeerNode(trees, subtree1);
			if (peerNodelist != null) {
				String newName = generateCollapsedName(peerNodelist.get(0));

				// BEGIN ADDITION by Derek
				// store the clade that we're removing
				cd._name2clade.put(newName, new STITree<Object>(peerNode1));
				// END ADDITION

				fixNode(peerNode1, newName);	// Collapse in tr1.
				for(TMutableNode peerNode: peerNodelist){
					fixNode(peerNode, newName);	// Collapse in tr2.
				}

				changed = true;
			}
		}

		return changed;
	}
*/
	/*
	protected static boolean runRule2(TMutableNode peerNode1, List<Tree> trees, CollapseDescriptor cd)
	{
		// Base case: nothing to do with a single leaf.
		if (peerNode1.isLeaf())
			return false;
		// Non-trivial case. Recurse on peerNode1's childrent before collapse itself.
		boolean changed = false;

		// Collapse peerNode1's children.
		Iterator it = peerNode1.getChildren().iterator();
		TMutableNode currentNode = null;

		while (it.hasNext()) {
			TMutableNode oldNode = currentNode;
			currentNode =  (TMutableNode) it.next();

			if (oldNode != null) {
				boolean tmp = runRule2(oldNode, trees, cd);
				changed = changed || tmp;
			}
		}
		if (currentNode != null) {
			boolean tmp = runRule2(currentNode, trees, cd);
			changed = changed || tmp;
		}

		List<String> leafChildNames = new ArrayList<String>();
		for (TNode child: peerNode1.getChildren()) {
			if (child.isLeaf()){
				leafChildNames.add(child.getName());
			}
		}

		if (leafChildNames.size() >= 2) {
			List<BitSet> collapseLeavesList = new ArrayList<BitSet>();
			List<List<TMutableNode>> peerNodesList = findPeerNode(trees, leafChildNames, collapseLeavesList);

			if (peerNodesList!=null) {

				for(int i=1;i<collapseLeavesList.size();i++){
					BitSet bs1 = collapseLeavesList.get(i);
					for(int j=0;j<i;j++){
						BitSet bs2 = collapseLeavesList.get(j);
						if(bs1.cardinality() < bs2.cardinality()){
							collapseLeavesList.remove(bs1);
							collapseLeavesList.add(j,bs1);
							break;
						}
					}
				}
				String[] newNameList = new String[collapseLeavesList.size()];
				boolean[] visitedlist = new boolean[collapseLeavesList.size()];
				for(int i=0;i<collapseLeavesList.size();i++){
					visitedlist[i] = false;
				}

				for(int i=0; i<collapseLeavesList.size();i++){
					BitSet bs = (BitSet)collapseLeavesList.get(i).clone();
					List<String> collapseNames = new ArrayList<String>();
					for(int j=0; j<i; j++){
						if(visitedlist[j]){
							continue;
						}
						BitSet prebs = collapseLeavesList.get(j);
						BitSet temp = (BitSet)prebs.clone();
						temp.and(bs);
						if(temp.equals(prebs)){
							collapseNames.add(newNameList[j]);
							visitedlist[j] = true;
							bs.andNot(prebs);
						}
					}

					for (int j = bs.nextSetBit(0); j >= 0; j = bs.nextSetBit(j+1)) {
					     collapseNames.add(leafChildNames.get(j));
					 }
					String newName = generateCollapsedName(collapseNames);
					newNameList[i] = newName;

					// store the clade that we're removing
					STITree<Object> newTree = new STITree<Object>(true);
					for(String leaf: collapseNames){
						newTree.getRoot().createChild(leaf);
					}
					cd._name2clade.put(newName, newTree);

					fixNode(peerNode1, collapseNames, newName);

					for(List<TMutableNode> nodesInOneTree: peerNodesList){
						boolean found = false;
						for(TMutableNode peerNode: nodesInOneTree){
							for(TNode child: peerNode.getChildren()){
								if(child.isLeaf() && collapseNames.contains(child.getName())){
									fixNode(peerNode, collapseNames, newName);
									found = true;
									break;
								}
							}
						}
						if(!found){
							throw new RuntimeException("Error!");
						}
					}
				}
				changed = true;
			}
		}

		return changed;
	}
	*/

	protected static boolean runRule2(TMutableNode peerNode1, List<Tree> trees, CollapseDescriptor cd)
	{
		// Base case: nothing to do with a single leaf.
		if (peerNode1.isLeaf())
			return false;
		// Non-trivial case. Recurse on peerNode1's childrent before collapse itself.
		boolean changed = false;

		// Collapse peerNode1's children.
		Iterator it = peerNode1.getChildren().iterator();
		TMutableNode currentNode = null;

		while (it.hasNext()) {
			TMutableNode oldNode = currentNode;
			currentNode =  (TMutableNode) it.next();

			if (oldNode != null) {
				boolean tmp = runRule2(oldNode, trees, cd);
				changed = changed || tmp;
			}
		}
		if (currentNode != null) {
			boolean tmp = runRule2(currentNode, trees, cd);
			changed = changed || tmp;
		}

		List<String> leafChildNames = new ArrayList<String>();
		for (TNode child: peerNode1.getChildren()) {
			if (child.isLeaf()){
				leafChildNames.add(child.getName());
			}
		}

		if (leafChildNames.size() >= 2) {
			List<BitSet> collapseLeavesList = new ArrayList<BitSet>();
			List<List<TMutableNode>> peerNodesList = findPeerNode(trees, leafChildNames, collapseLeavesList);

			if (peerNodesList!=null) {

				for(int i=1;i<collapseLeavesList.size();i++){
					BitSet bs1 = collapseLeavesList.get(i);
					for(int j=0;j<i;j++){
						BitSet bs2 = collapseLeavesList.get(j);
						if(bs1.cardinality() < bs2.cardinality()){
							collapseLeavesList.remove(bs1);
							collapseLeavesList.add(j,bs1);
							break;
						}
					}
				}
				String[] newNameList = new String[collapseLeavesList.size()];
				boolean[] visitedlist = new boolean[collapseLeavesList.size()];
				for(int i=0;i<collapseLeavesList.size();i++){
					visitedlist[i] = false;
				}

				for(int i=0; i<collapseLeavesList.size();i++){
					BitSet bs = (BitSet)collapseLeavesList.get(i).clone();
					List<String> collapseNames = new ArrayList<String>();
					for(int j=0; j<i; j++){
						if(visitedlist[j]){
							continue;
						}
						BitSet prebs = collapseLeavesList.get(j);
						BitSet temp = (BitSet)prebs.clone();
						temp.and(bs);
						if(temp.equals(prebs)){
							collapseNames.add(newNameList[j]);
							visitedlist[j] = true;
							bs.andNot(prebs);
						}
					}

					for (int j = bs.nextSetBit(0); j >= 0; j = bs.nextSetBit(j+1)) {
					     collapseNames.add(leafChildNames.get(j));
					 }
					String newName = generateCollapsedName(collapseNames);
					newNameList[i] = newName;

					// store the clade that we're removing
					STITree<Object> newTree = new STITree<Object>(true);
					for(String leaf: collapseNames){
						newTree.getRoot().createChild(leaf);
					}
					cd._name2clade.put(newName, newTree);

					fixNode(peerNode1, collapseNames, newName);

					for(List<TMutableNode> nodesInOneTree: peerNodesList){
						boolean found = false;
						for(TMutableNode peerNode: nodesInOneTree){
							for(TNode child: peerNode.getChildren()){
								if(child.isLeaf() && collapseNames.contains(child.getName())){
									fixNode(peerNode, collapseNames, newName);
									found = true;
									break;
								}
							}
						}
						if(!found){
							throw new RuntimeException("Error!");
						}
					}
				}
				changed = true;
			}
		}

		return changed;
	}

	/**
	 * Tests if this tree is pendant, i.e. its immediate children are leaves.
	 *
	 * Require:
	 *     Tree tr is rooted.
	 *
	 * @param tr: The tree to be checked if it's pendant.
	 *
	 * @return
	 *     true if this tree is pedant; false if it's not.
	 */
	protected static boolean isPendant(Tree tr)
	{
		Iterator it = tr.getRoot().getChildren().iterator();
		boolean pendant = true;

		// Test if all of the tree's children are leaves.
		while (it.hasNext() && pendant) {
			TNode node = (TNode) it.next();

			if (!node.isLeaf())
				pendant = false;
		}

		return pendant;
	}

	/**
	 * Finds a peer node in the tree tr whose subtree is isormorphic to pendantStr.
	 * Require:
	 *     As this function is used by the collapse function, we only provide the implemetation
	 *     in the case where pendantStr is a pendant subtree.
	 * @param
	 *     pendantStr: A pendant subtree, i.e. its immediate children are leaves.
	 * @return
	 *     a node in this tree if pendantStr is actually its subtree; null otherwise.
	 *
	 * NOTE: This method is used in collpasing two trees.
	 */
	protected static TMutableNode findPeerNode(Tree tr, Tree pendantSubtree)
	{
		// Get the name of first child of pendantStr.
		Iterator it = pendantSubtree.getRoot().getChildren().iterator();
		String name = ((TNode) it.next()).getName();

		// Get the parent of this node in this tree.
		TMutableNode parent = (TMutableNode) tr.getNode(name).getParent();

		// Compare the pendant tree with the subtree rooted at node parent.
		STITree peerSubtree = new STITree(parent, true);
		if (isIdentical(pendantSubtree, peerSubtree))
			return parent;
		else
			return null;
	}

	/*
	protected static List<TMutableNode> findPeerNode(List<Tree> trees, Tree pendantSubtree)
	{
		// Get the name of first child of pendantStr.
		Iterator it = pendantSubtree.getRoot().getChildren().iterator();
		String name = ((TNode) it.next()).getName();
		List<TMutableNode> parentlist = new ArrayList<TMutableNode>();

		for(Tree tr: trees){
			// Get the parent of this node in this tree.
			TMutableNode parent = (TMutableNode) tr.getNode(name).getParent();

			// Compare the pendant tree with the subtree rooted at node parent.
			STITree peerSubtree = new STITree(parent, true);
			if (isIdentical(pendantSubtree, peerSubtree))
				parentlist.add(parent);
			else
				return null;
		}

		return parentlist;
	}

	*/


	protected static List<List<TMutableNode>> findPeerNode(List<Tree> trees, List<String> childNames, List<BitSet> collapseLeavesList)
	{
		List<List<TMutableNode>> collapseNodesList = new ArrayList<List<TMutableNode>>();
		BitSet unionbs = new BitSet(childNames.size());
		unionbs.set(0, childNames.size());
		collapseLeavesList.add(unionbs);

		for(Tree tr: trees){
			List<TMutableNode> collapseNodes = new ArrayList<TMutableNode>();
			int index = 0;
			List<BitSet> addlist = new ArrayList<BitSet>();
			List<BitSet> checkedlist = new ArrayList<BitSet>();
			while(index < collapseLeavesList.size()){
				if(index == 0){
					for(int i=1;i<collapseLeavesList.size();i++){
						BitSet bs1 = collapseLeavesList.get(i);
						for(int j=0;j<i;j++){
							BitSet bs2 = collapseLeavesList.get(j);
							if(bs1.cardinality() < bs2.cardinality()){
								collapseLeavesList.remove(bs1);
								collapseLeavesList.add(j,bs1);
								break;
							}
						}
					}
				}
				BitSet collapseLeaves = collapseLeavesList.get(index);
				if(checkedlist.contains(collapseLeaves)){
					index++;
					continue;
				}
				checkedlist.add(collapseLeaves);
				Map<TNode, BitSet> node2leaves = findSharedClade(tr, childNames, collapseLeaves);
				Collection<BitSet> bsset = node2leaves.values();
				for(TNode node: node2leaves.keySet()){
					if(!collapseNodes.contains(node)){
						collapseNodes.add((TMutableNode)node);
					}
				}
				boolean decompose = true;
				for(BitSet bs: bsset){
					if(bs.cardinality()==collapseLeaves.cardinality()){
						decompose = false;
					}
					else{
						addlist.add(bs);
					}
				}
				if(decompose){
					List<BitSet> rmlist = new ArrayList<BitSet>();
					rmlist.add(collapseLeaves);
					int size = collapseLeavesList.size();
					for(int i=index+1; i<size; i++){
						BitSet obs = collapseLeavesList.get(i);
						if(bs1Containbs2(obs, collapseLeaves)){
							BitSet temp = (BitSet)obs.clone();
							for(BitSet rmbs: rmlist){
								temp.andNot(rmbs);
							}
							if(!collapseLeavesList.contains(temp) && temp.cardinality()>=2){
								collapseLeavesList.add(temp);
								index = 0;
							}
							rmlist.add(obs);
						}
					}
					collapseLeavesList.removeAll(rmlist);
				}
				else{
					index++;
				}
			}
			for(BitSet bs: addlist){
				if(!collapseLeavesList.contains(bs)){
					collapseLeavesList.add(bs);
				}
			}
			collapseNodesList.add(collapseNodes);
			if(collapseLeavesList.size()==0){
				return null;
			}
		}

		return collapseNodesList;
	}



	private static boolean bs1Containbs2(BitSet bs1, BitSet bs2){
		BitSet temp = (BitSet)bs1.clone();
		temp.and(bs2);
		if(temp.equals(bs2)){
			return true;
		}else{
			return false;
		}
	}


	protected static Map<TNode, BitSet> findSharedClade(Tree tr, List<String> leafNamesList, BitSet peerbs){
		Map<TNode, BitSet> commonLeaves = new HashMap<TNode, BitSet>();
		Map<TNode, Boolean> isContained = new HashMap<TNode, Boolean>();
		int totalSize = leafNamesList.size();

		for (TNode node : tr.postTraverse()) {
			if (node.isLeaf()) {
				int index = leafNamesList.indexOf(node.getName());
				if(index==-1){
					isContained.put(node, false);
				}
				else{
					if(peerbs.get(index)){
						BitSet bs = new BitSet(totalSize);
						bs.set(index);
						commonLeaves.put(node, bs);
						isContained.put(node, true);
					}
					else{
						isContained.put(node, false);
					}
				}
			}
			else {
				boolean contain = true;
				BitSet bs = new BitSet(totalSize);
				int count = 0;
				for (TNode child : node.getChildren()) {
					boolean childContained = isContained.get(child);
					contain = contain && childContained;
					if(childContained){
						bs.or(commonLeaves.get(child));
						if(child.isLeaf()){
							commonLeaves.remove(child);
						}
						count++;
					}
					isContained.remove(child);
				}
				isContained.put(node, contain);
				if(count>=2){
					commonLeaves.put(node, bs);
				}
			}
		}

		return commonLeaves;
	}





	/**
	 * This method checks if the two trees tr1 and tr2 are identical. For the purpose of
	 * collapsing two trees, this function will check ONLY pendant trees. It'll always
	 * return false if either tr1 and tr2 is not pedant, even if they are really isomorphic.
	 *
	 * @param
	 *     tr1: Pendant tree 1.
	 * @param
	 *     tr2: Pendant tree 2.
	 * @return
	 *     true if tr1 and tr2 are pendant and isomorphic; false otherwise.
	 */
	protected static boolean isIdentical(Tree tr1, Tree tr2)
	{
		// Make sure that we only compare two trees with the same number of leaves and
		// the same number of internal nodes.
		if (tr1.getLeafCount() != tr2.getLeafCount() || tr1.getNodeCount() != tr2.getNodeCount())
			return false;

		// Forms the array of names of the root's children.
		TNode root = tr1.getRoot();
		String childNames[] = new String[root.getChildCount()];
		Iterator it = root.getChildren().iterator();

		for (int i = 0; it.hasNext(); i++) {
			TNode child = (TNode) it.next();
			if (!child.isLeaf())	// tr1 is not pendant.
				return false;
			else
				childNames[i] = child.getName();
		}

		// Is every child of pendantStr in this tree?
		boolean identical = true;
		it = tr2.getRoot().getChildren().iterator();

		while (it.hasNext() && identical) {
			TNode child = (TNode) it.next();
			String name = child.getName();
			int j;

			if (!child.isLeaf()) {	// tr2 is pendant.
				return false;
			}
			else {
				// Find this child in tr1's children.
				for (j = 0; j < childNames.length; j++) {
					if (name.equals(childNames[j]))
						break;
				}
				if (j >= childNames.length)
					identical = false;
			}
		}

		return identical;
	}


	/**
	 * Creates the new name for the TNode node. The new name is the concatenation of (sorted)
	 * names of node's children.
	 * Require:
	 *     The subtree rooted at node must be a pendant tree.
	 * @return
	 *     The new name of the node.
	 */
	protected static String generateCollapsedName(TNode node)
	{
		// Get the names of its children.
		String childNames[] = new String[node.getChildCount()];
		Iterator it = node.getChildren().iterator();

		for (int i = 0; it.hasNext(); i++) {
			TNode child = (TNode) it.next();
			childNames[i] = child.getName();
		}

		// Concatenate names from childNames to create a new name for the TNode node.
		String newName = COLLAPSED_NAME_PREFIX;
		Arrays.sort(childNames);
		for (int i = 0; i < childNames.length; i++) {
			if (i > 0) {
				newName += "_";
			}
			newName += childNames[i];
		}
		return newName;
	}

	protected static String generateCollapsedName(List<String> collapseNames)
	{
		// Get the names of its children.
		Object childNames[] = collapseNames.toArray();

		// Concatenate names from childNames to create a new name for the TNode node.
		String newName = COLLAPSED_NAME_PREFIX;
		Arrays.sort(childNames);
		for (int i = 0; i < childNames.length; i++) {
			if (i > 0) {
				newName += "_";
			}
			newName += childNames[i];
		}

		return newName;
	}

	/**
	 * This method actually deletes all children of the TNode node, and renames its name
	 * as the concatenation of the names of the children. Note that the names are sorted
	 * alphabetically before they are concatenated.
	 *
	 * Require:
	 *     Again, node must be the the root of a pedant subtree.
	 * @param
	 *     node: The nod to be fixed.
	 * @param
	 *     newName: Name for the new node.
	 *
	 * NOTE: This method is used in the collapse procedure.
	 */
	protected static void fixNode(TMutableNode node, String newName)
	{
		// Remove children of this node.
		Iterator it = node.getChildren().iterator();
		while (it.hasNext()) {
			TMutableNode child = (TMutableNode) it.next();
			node.removeChild(child, false);
			it = node.getChildren().iterator();
		}

		// Set the node's name as newName.
		node.setName(newName);
	}

	protected static void fixNode(TMutableNode node, List<String> collapseLeaves, String newName)
	{
		// Remove children of this node.
		Iterator it = node.getChildren().iterator();
		int leavesLeft = collapseLeaves.size();
		while (leavesLeft!=0) {
			TMutableNode child = (TMutableNode) it.next();
			if(collapseLeaves.contains(child.getName())){
				node.removeChild(child, false);
				it = node.getChildren().iterator();
				leavesLeft--;
			}
		}

		if(node.getChildCount()==0){
			node.setName(newName);
			((STINode<Integer>)node).setData(0);
		}
		else{
			TMutableNode child = node.createChild(newName);
			((STINode<Integer>)child).setData(-1);
		}
		// Set the node's name as newName.

	}

}
