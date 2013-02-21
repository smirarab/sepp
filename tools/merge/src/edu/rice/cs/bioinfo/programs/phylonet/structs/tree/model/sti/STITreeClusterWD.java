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
import java.util.List;



public class STITreeClusterWD<D extends Object> extends STITreeCluster {
	//private STITreeClusterWT _cluster;
	private D _data;

	public STITreeClusterWD(STITreeCluster tc) {
		super(tc);
		_data = null;
	}
	
	public STITreeClusterWD(STITreeClusterWD<D> tc) {
		super(tc);
		_data = tc.getData();
	}
			
	public STITreeClusterWD(STITreeCluster tc,D data) {
		super(tc);
		_data = data;
	}

	public STITreeClusterWD(String[] taxa){
		super(taxa);
		_data = null;
	}
	
	public STITreeClusterWD<D> duplicate(){
		STITreeClusterWD<D> newCluster = new STITreeClusterWD<D>(_taxa);
		newCluster.setCluster(_cluster);
		newCluster.setData(_data);
		return newCluster;
	}
	
	//only work for binary tree
	public boolean canMakeBranch(STITreeClusterWD<D> c2, List<STITreeClusterWD<D>> clusters){
		boolean makeBranch = false;
		if(c2.containsCluster(this)){
			if(c2.getClusterSize()-_cluster.cardinality() == 1){
				makeBranch = true;
			}
			else{
				BitSet temp = (BitSet)(_cluster.clone());
				temp.xor(c2.getCluster());
				
				for(STITreeClusterWD<D> c: clusters){
					if(c.getCluster().equals(temp)){
						makeBranch = true;
						break;
					}
				}
			}
		}
		return makeBranch;
	}

	
	public boolean equals(Object o){
		assert(_taxa != null && _taxa.length > 0);
		if(!(o instanceof STITreeCluster)){
			return false;
		}
		STITreeClusterWD<D> tc = (STITreeClusterWD<D>) o;
		if (tc == null || tc._cluster == null) {
			//System.err.println("Cluster is null. The function returns false.");
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
	
	public void setData(D data){
		_data = data;
	}
	
	public D getData(){
		return _data;
	}
	
	public boolean containsCluster(STITreeClusterWD<D> tc) {
		return containsCluster(tc._cluster);
	}
	
	public String toString(){
		StringBuffer out = new StringBuffer();
		out.append(super.toString());
		out.append(_data.toString());
		return out.toString();
	}

}
