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

import java.util.BitSet;

/**
 * A class for representing a tree cluster.
 *
 * @author cvthan
 *
 */
public class STITreeCluster {
	protected String _taxa[];		// Set of taxa of a tree.
	protected BitSet _cluster;	// The cluster.


	//protected STITreeCluster(){}


	/**
	 * Create an empty cluster.
	 *
	 * @param taxa
	 */
	public STITreeCluster(String taxa[]) {
		if (taxa == null || taxa.length == 0) {
			System.err.println("Invalid cluster");

			_taxa = null;
			_cluster = null;
			return;
		}

		_taxa = taxa;
		_cluster = new BitSet(taxa.length);
	}

	/**
	 * Make the same copy of this cluster
	 *
	 * @param tc
	 */
	public STITreeCluster(STITreeCluster tc) {
		assert(tc._taxa != null && tc._taxa.length > 0);

		_taxa = tc._taxa;
		_cluster = new BitSet(_taxa.length);
		_cluster.or(tc._cluster);
	}

	public String[] getTaxa() {
		return _taxa;
	}

	/**
	 * Set a bit vector that represents this cluster.
	 *
	 * @param c
	 */
	public void setCluster(BitSet c) {
		if (c != null) {
			_cluster = c;
		}
		else {
			System.err.println("Null bit set.");
		}
	}

	public BitSet getCluster() {
		return _cluster;
	}

	public int getClusterSize() {
		return _cluster.cardinality();
	}

	/**
	 * Get the set of leaves in this cluster.
	 *
	 * @return
	 */
	public String[] getClusterLeaves() {
		assert(_taxa != null && _taxa.length > 0);

		String cl[] = new String[_cluster.cardinality()];
		int c = 0;
		for (int i = 0; i < _cluster.length(); i++) {
			if (_cluster.get(i)) {
				cl[c++] = _taxa[i];
			}
		}

		return cl;
	}

	/**
	 * Add a leaf to this cluster.
	 *
	 * @param l
	 */
	public void addLeaf(String l) {
		int i = 0;
		for (i = 0; i < _taxa.length; i++) {
			if (l.equals(_taxa[i])) {
				break;
			}
		}

		if (i < _taxa.length) {
			_cluster.set(i);
		}
	}


	/**
	 * Check if another cluster is equal to this cluster. This function assumes two clusters in
	 * comparison use the same set of taxa as well as taxon orderings.
	 *
	 * @param Another cluster.
	 * @return
	 */
	public boolean equals(Object o) {

		assert(_taxa != null && _taxa.length > 0);
		if(!(o instanceof STITreeCluster)){
			return false;
		}

		STITreeCluster tc = (STITreeCluster) o;
		if (tc == null || tc._cluster == null) {
			System.err.println("Cluster is null. The function returns false.");
			return false;
		}

		if (_cluster.equals(tc._cluster)) {
			return true;
		}
		else {
			return false;

		}
	}

	public int hashCode(){
		return _cluster.hashCode()+_taxa.hashCode();
	}

	/**
	 * Check if another cluster is compatible with this cluster. This function assumes that the two
	 * clusters in comparison use the same set of taxa as well as taxon orderings.
	 *
	 * @param tc
	 * @return
	 */
	public boolean isCompatible(STITreeCluster tc) {
		assert(_taxa != null && _taxa.length > 0);

		if (tc == null || tc._cluster == null) {
			System.err.println("Cluster is null. The function returns false.");
			return false;
		}

		BitSet temp = (BitSet) _cluster.clone();
		temp.and(tc._cluster);
		if (temp.equals(_cluster) || temp.equals(tc._cluster)) {
			// One cluster is a subset of the other.
			return true;
		}
		else if (temp.cardinality() == 0) {
			// Two clusters are disjoint.
			return true;
		}
		else {
			return false;
		}
	}

	/**
	 * Check if this cluster is disjoint with another cluster. This function assumes that the two clusters
	 * in comparison use the same set of tax and taxon orderings.
	 *
	 * @param tc
	 * @return
	 */
	public boolean isDisjoint(STITreeCluster tc) {
		assert(_taxa != null && _taxa.length >0);

		if (tc == null || tc._cluster == null) {
			System.err.println("Cluster is null. The function returns false.");
			return false;
		}

		return isDisjoint(tc._cluster);
	}

	public boolean isDisjoint(BitSet tc) {
		assert(_taxa != null && _taxa.length >0);

		if (tc == null) {
			System.err.println("Cluster is null. The function returns false.");
			return false;
		}

		BitSet temp = (BitSet) _cluster.clone();
		temp.and(tc);

		return (temp.cardinality() == 0);
	}


	/**
	 * Check if two clusters are complementary to each other (i.e., they are a partition of the taxon set).
	 *
	 * @param tc
	 * @return
	 */
	public boolean isComplementary(STITreeCluster tc) {
		assert(_taxa != null && _taxa.length > 0);

		if (tc == null || tc._cluster == null) {
			System.err.println("Cluster is null. The function returns false.");
			return false;
		}

		BitSet temp1 = (BitSet) _cluster.clone();
		temp1.and(tc._cluster);

		BitSet temp2 = (BitSet) _cluster.clone();
		temp2.or(tc._cluster);

		if (temp1.cardinality() == 0 && temp2.cardinality() == _taxa.length) {
			return true;
		}
		else {
			return false;
		}
	}

	/**
	 * Return true if the cluster contains leaf l.
	 *
	 * @param l
	 * @return
	 */
	public boolean containsLeaf(String l) {
		int i;
		for (i = 0; i < _taxa.length; i++) {
			if (_taxa[i].equals(l)) {
				break;
			}
		}

		return _cluster.get(i);
	}

	/**
	 * Return true if this cluster contains another cluster.
	 *
	 * @param tc
	 * @return
	 */
	public boolean containsCluster(STITreeCluster tc) {
		BitSet temp = (BitSet) tc._cluster.clone();
		temp.and(_cluster);

		return (temp.equals(tc._cluster));
	}

	/**
	 * Return true if this cluster contains cluster bs. Assume that bs use the same taxa-index reference.
	 *
	 * @param bs
	 * @return
	 */
	public boolean containsCluster(BitSet bs) {
		BitSet temp = (BitSet) bs.clone();
		temp.and(_cluster);

		return (temp.equals(bs));
	}

	/**
	 * Merge with another cluster.
	 *
	 * @param tc
	 * @return
	 */
	public STITreeCluster merge(STITreeCluster tc) {
		STITreeCluster temp = new STITreeCluster(this);
		temp._cluster.or(tc._cluster);

		return temp;
	}

	public STITreeCluster complementaryCluster(){
		STITreeCluster cc = new STITreeCluster(_taxa);
		BitSet bs = (BitSet)_cluster.clone();
		for(int i=0; i<_taxa.length; i++){
			if(bs.get(i)){
				bs.set(i,false);
			}
			else{
				bs.set(i,true);
			}
		}
		cc.setCluster(bs);
		return cc;
	}

	/**
	 * Returns the string representation of this cluster.
	 *
	 */
	public String toString() {
		StringBuffer out = new StringBuffer();

		out.append("{");
		for (String s : getClusterLeaves()) {
			out.append(s + ", ");
		}
		out.delete(out.length() - 2, out.length());
		out.append("}");

		return out.toString();
	}

//	public STITreeCluster() {
//		_leaves = null;
//		_cluster = null;
//	}
//
//	/**
//	 * Creates a new empty cluster.
//	 *
//	 * @param leaves
//	 */
//	public STITreeCluster(String leaves[]) {
//		_leaves = leaves;
//		_cluster = null;
//	}
//
//	public STITreeCluster(String leaves[], BitSet c) {
//		_leaves = leaves;
//		_cluster = c;
//	}
//
//	public void setLeaves(String leaves[]) {
//		_leaves = leaves;
//	}
//
//	public String[] getLeaves() {
//		return _leaves;
//	}
//
//	public void setCluster(BitSet c) {
//		_cluster = c;
//	}
//
//	public BitSet getCluster() {
//		return _cluster;
//	}
}