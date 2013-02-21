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

import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;

import edu.rice.cs.bioinfo.programs.phylonet.structs.BitVector;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;



/**
 * This class contains methods related to the bipartitions contained within a tree.<BR>
 * Unless otherwise noted, bipartitions are represented using {@link BitVector}s in which each
 * bit represents a leaf in the tree.  If the bit is set, then the corresponding leaf is in the
 * bipartition.
 * 
 * @author Derek Ruths
 */
public class Bipartitions {

	/**
	 * This method computes all the bipartitions contained in a tree and stores them and the node
	 * below the edge that induces them in the map <code>h</code>.
	 * 
	 * @param tree
	 * @param leaf_assignment is the positional assignment of the leaves of the tree to positions
	 * in the bitvector.  If this hashtable is empty, this method will populate it with an
	 * assignment.
	 * @param h
	 */
	public static final void computeBipartitions(Tree tree, Map<String,Integer> leaf_assignment, Map<BitVector,TNode> h) {

        int leafAssignmentCount = leaf_assignment.size();

		if(leafAssignmentCount == 0) {
			assignLeafPositions(tree, leaf_assignment);
		}
        else
        {
            int treeLeafCount = tree.getLeafCount();
            if(leafAssignmentCount != tree.getLeafCount())
            {
			    throw new RuntimeException("leaf_assignment contains an incorrect number of entries");
            }
        }
		
		// skip the root
		for(TNode child : tree.getRoot().getChildren()) {
			addBipartitions(child, h, leaf_assignment);
		}
	}
	
	public static final void computeBipartitions(Tree tree, String[] taxa, Map<BitVector,TNode> h){
		Map<String, Integer> leaf_assignment = new HashMap<String, Integer>();
		for(int i=0;i<taxa.length;i++){
			leaf_assignment.put(taxa[i], i);
		}
		computeBipartitions(tree,leaf_assignment,h);
	}

	/**
	 * This method computes all the bipartitions contained in a tree and returns them in 
	 * an array.
	 * 
	 * @param tree
	 * @param leaf_assignment is the positional assignment of the leaves of the tree to positions
	 * in the bitvector.  If this hashtable is empty, this method will populate it with an
	 * assignment.
	 */
	public static final BitVector[] computeBipartitions(Tree tree, Map<String,Integer> leaf_assignment) {
		
		if(leaf_assignment.size() == 0) {
			assignLeafPositions(tree, leaf_assignment);
		} else if(leaf_assignment.size() != tree.getLeafCount()) {
			throw new RuntimeException("leaf_assignment contains an incorrect number of entries");
		}
		
		Hashtable<BitVector,TNode> hresult = new Hashtable<BitVector,TNode>();
		computeBipartitions(tree, leaf_assignment, hresult);
		
		BitVector[] bvarray = new BitVector[hresult.size()];
		
		int i = 0;
		for(BitVector bv : hresult.keySet()) {
			bvarray[i++] = bv;
		}
		
		return bvarray;
	}

	/**
	 * Generates an assignment of the leaves of the tree to positions in a {@link BitVector}.
	 * 
	 * @param tree
	 * @return
	 */
	public static final Map<String,Integer> assignLeafPositions(Tree tree) {
		Hashtable<String,Integer> assignments = new Hashtable<String,Integer>();
		
		assignLeafPositions(tree, assignments);
		
		return assignments;
	}
	
	/**
	 * Generate an assignment of leaves to positions in a {@link BitVector}.
	 * 
	 * @param tree is the tree whose leaves will be assigned positions.
	 * @param leaf_assignments is the object in which the position assignments will be stored.
	 * Because assignments are added to this map, the map should be empty when called by this
	 * method.
	 */
	public static final void assignLeafPositions(Tree tree, Map<String,Integer> leaf_assignments) {

		for(TNode n : tree.getNodes()) {
			if(n.isLeaf()) {
				leaf_assignments.put(n.getName(), leaf_assignments.size());
			}
		}
	}
	
	private static final BitVector addBipartitions(TNode n, Map<BitVector,TNode> h, Map<String,Integer> leaf_assignment) {
		
		// if this is a binary node, skip it - it doesn't partition anything
		if(n.getChildCount() == 1) {
			return addBipartitions(n.getChildren().iterator().next(), h, leaf_assignment);
		}
		
		BitVector bv = new BitVector(leaf_assignment.size());
		
		if(n.isLeaf()) {
			bv.setValue(leaf_assignment.get(n.getName()),true);
		} else {
			for(TNode child : n.getChildren()) {
				bv.or(addBipartitions(child, h, leaf_assignment));
			}
		}
		
		// only store the bipartition if it isn't in the array at present
		bv.not();
		/*
		if(!h.containsKey(bv)) {
			bv.not();
			h.put(bv, n);
		} else {
			if(bv.countOnes()==1){
				h.remove(bv);
				bv.not();
				h.put(bv, n);
			}
			else{
				bv.not();
			}
		}
		*/
		
		if(!h.containsKey(bv)) {
			bv.not();
			h.put(bv, n);
		} else {
			bv.not();
		}
		
		return bv;
	}
	
	/**
	 * @return a string array in which the string in position <code>i</code> is the leaf name that was
	 * assigned to position <code>i</code> in the partition bit vector. 
	 */
	public static final String[] getPositionToLeafNameMap(Map<String,Integer> leaf_assignment) {
		String[] r = new String[leaf_assignment.size()];
		
		for(Map.Entry<String,Integer> e : leaf_assignment.entrySet()) {
			r[e.getValue()] = e.getKey();
		}
		
		return r;
	}
	
	public static final boolean isCompatible(BitVector bv1,BitVector bv2){
		BitVector A1 = new BitVector(bv1);

		//A1&A2
		A1.and(bv2);
		if(A1.countOnes() == 0){  
			return true;
		}
		//A1&B2
		A1.setValue(bv1);
		bv2.not();
		A1.and(bv2);
		if(A1.countOnes() == 0){  
			bv2.not();
			return true;
		}
		//B1&B2
		bv1.not();
		A1.setValue(bv1);
		A1.and(bv2);
		if(A1.countOnes() == 0){  
			bv1.not();
			bv2.not();
			return true;
		}
		//B1&A2
		bv2.not();
		A1.setValue(bv1);
		A1.and(bv2);
		if(A1.countOnes() == 0){  
			bv1.not();
			return true;
		}
		bv1.not();
		return false;
	}
}
