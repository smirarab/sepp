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

import java.util.*;

/**
 * A class for representing a tree bipartition. The bipartition is represented by two bit sets--left set 
 * and right set--and denoted as LEFT | RIGHT. 
 * 
 * @author cvthan
 *
 */
public class STITreeBipartition {
	private String _taxa[];		// Pointer to the leaf assignment.
	private BitSet _left, _right;	// The two sets of leaves of a bipartition. They are complements to each other.
	
	/** 
	 * Create a trivial bipartition, that is, \emptyset | X.
	 */
	public STITreeBipartition(String taxa[]) {
		if (taxa == null || taxa.length == 0) {
			System.err.println("Invalid biparition");
			
			_taxa = null;
			_left = _right = null;
			return;
		}
		
		_taxa = taxa;
		_left = new BitSet(_taxa.length);
		_right = new BitSet(_taxa.length);
		for (int i = 0; i < _taxa.length; i++) {
			_right.set(i);
		}
	}
	
	public String[] getTaxa() {
		return _taxa;
	}
	
	/** 
	 * Add a leaf to the left half of the bipartition LEFT | RIGHT.
	 * 
	 * @param l
	 */
	public void addLeafToLeft(String l) {
		assert(_taxa != null && _taxa.length > 0);
		
		int i = 0;
		for (i = 0; i < _taxa.length; i++) {
			if (l.equals(_taxa[i])) {
				break;
			}
		}
		
		if (i < _taxa.length) {
			_left.set(i);
			_right.clear(i);
		}
	}
		
	/** 
	 * Add a leaf to the right half of the bipartition LEFT | RIGHT.
	 * 
	 * @param l
	 */
	public void addLeafToRight(String l) {
		assert(_taxa != null && _taxa.length > 0);
		
		int i = 0;
		for (i = 0; i < _taxa.length; i++) {
			if (l.equals(_taxa[i])) {
				break;
			}
		}
		
		if (i < _taxa.length) {
			_left.clear(i);
			_right.set(i);
		}
	}
	
	/** 
	 * Set left and right halves of this bipartition. Only one of the parameters is required; the other 
	 * is optional. However, if both are present (i.e., not null), then they must be the complement of each other.
	 * 
	 * @param left
	 * @param right
	 */
	public void setBipartition(BitSet left, BitSet right) {
		assert(_taxa != null && _taxa.length > 0);
		
		if (left != null && right != null) {
			BitSet temp = (BitSet) left.clone();
			temp.and(right);
			
			if (temp.cardinality() > 0) {
				System.err.println("Invalid arguments. The two bit vectors must be the complement of each other.");
				return;
			}

			if (left.length() > _taxa.length || right.length() > _taxa.length) {
				System.err.println("Invalid arguments. One of the bipartition contains leaves not in the set of taxa");
				return;
			}
			
			if (left.cardinality() + right.cardinality() != _taxa.length) {
				System.err.println("Invalid arguments. The two halves of the bipartition must make up the taxon set.");
				return;
			}
			
			_left = left;
			_right = right;
		}
		else {
			// Either one of them is null.
			if (left == null) {
				_right = right;
				for (int i = 0; i < _taxa.length; i++) {
					_left.set(i, !_right.get(i));
				}
			}
			else {
				// right is null.
				assert(right == null);
				_left = left;
				for (int i = 0; i < _taxa.length; i++) {
					_right.set(i, !_left.get(i));
				}
			}
		}
	}
	
	/** 
	 * 
	 * @return The bit vector representation of the left half.
	 */
	public BitSet getLeft() {
		return _left;
	}
	
	/** 
	 * @return The bit vector representation of the right half.
	 */
	public BitSet getRight() {
		return _right;
	}
	
	/** 
	 * Get the leaves in the left half of the bipartition LEFT | RIGHT.
	 * 
	 * @return
	 */
	public String[] getLeftLeaves() {
		assert(_taxa != null && _taxa.length > 0);
		
		String ll[] = new String[_left.cardinality()];
		int c = 0;
		for (int i = 0; i < _left.length(); i++) {
			if (_left.get(i)) {
				ll[c++] = _taxa[i];
			}
		}
		
		return ll;
	}
	
	/** 
	 * Get the leaves in the right half of the bipartition LEFT | RIGHT.
	 * 
	 * @return
	 */
	public String[] getRightLeaves() {
		assert(_taxa != null && _taxa.length > 0);
		
		String rl[] = new String[_right.cardinality()];
		int c = 0;
		for (int i = 0; i < _right.length(); i++) {
			if (_right.get(i)) {
				rl[c++] = _taxa[i];
			}
		}
		
		return rl;
	}
	
	/** 
	 * Checks if two bipartitions are equal or not. The function assumes that the two bipartitions are of 
	 * the same set of taxa, and they use the same ordering of the taxa.
	 * 
	 * @param tb
	 * @return
	 */
	public boolean isEqual(STITreeBipartition tb) {
		assert(_taxa != null && _left != null && _right != null);
		
		if (tb == null || tb._left == null || tb._right == null) {
			System.err.println("The bipartition is null. The function returns false.");
			return false;
		}
		
		if (_left == tb._left || _left == tb._right) {
			return true;
		}
		else {
			return false;
		}
	}
	
	/** 
	 * Checks if two bipartitions are compatible or not. This function assumes that two bipartitions in 
	 * comparison use the same set of taxa as well as the same taxon orderings.
	 * 
	 * @param tb
	 * @return
	 */
	public boolean isCompatible(STITreeBipartition tb) {
		assert(_taxa != null && _left != null && _right != null);
		
		if (tb == null || tb._left == null || tb._right == null) {
			System.err.println("The bipartition is null. The function returns false.");
			return false;
		}
		
		BitSet temp = (BitSet) _left.clone();
		temp.and(tb._left);
		if (temp.cardinality() == 0) {
			return true;
		}
		
		temp = (BitSet) _left.clone();
		temp.and(tb._right);
		if (temp.cardinality() == 0) {
			return true;
		}
		
		temp = (BitSet) _right.clone();
		temp.and(tb._left);
		if (temp.cardinality() == 0) {
			return true;
		}
		
		temp = (BitSet) _right.clone();
		temp.and(tb._right);
		if (temp.cardinality() == 0) {
			return true;
		}
		
		return false;
	}
	
	/** 
	 * Returns the string representation of this biapartition.
	 * 
	 */
	public String toString() {
		StringBuffer out = new StringBuffer();
		
		out.append("{");
		for (String s : getLeftLeaves()) {
			out.append(s + ", ");
		}
		out.delete(out.length() - 2, out.length());
		out.append("}");
		
		out.append(" | ");
		
		out.append("{");
		for (String s : getRightLeaves()) {
			out.append(s + ", ");
		}
		out.delete(out.length() - 2, out.length());
		out.append("}");
		
		return out.toString();
	}
}
