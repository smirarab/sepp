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

package edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeClusterWD;

import java.util.List;

/**
 * This interface defines the methods that all trees must implement.
 *
 * @author Derek Ruths
 */
public interface Tree  {

	/**
	 * When used with the {@link Printable.toString} method, the tree prints itself in
	 * newick format.
	 */
	public static final int NEWICK_FORMAT = 0;

	/**
	 * The value returned by {@link getName} when the tree does not have a name.
	 */
	public static final String NO_NAME = "";

	/**
	 * @return the name of the tree.
	 */
	public String getName();

	/**
	 * @return <code>true</code> if the tree contains no nodes.
	 */
	public boolean isEmpty();

	/**
	 * @return <code>true</code> if the tree is rooted.
	 */
	public boolean isRooted();

	/**
	 * If the tree is rooted, this returns the root.  If the tree is not rooted, then the
	 * method simply returns an internal node in the tree.  This will be the same tree always
	 * returned for this instance of the tree.  Parent/child relationships are configured with
	 * regard to this node.
	 */
	public TNode getRoot();

	/**
	 * Retrieve a node by its unique ID.
	 */
	public TNode getNode(int id);

	/**
	 * Retrieve a node by its name.
	 */
	public TNode getNode(String name);

	/**
	 * @return the number of nodes in the tree.
	 */
	public int getNodeCount();

	/**
	 * @return the number of leaves in the tree.
	 */
	public int getLeafCount();

	/**
	 * @return the set of leaf labels
	 */
	public String[] getLeaves();

	/**
	 * @return the an iterable list of nodes in the tree.
	 */
	public Iterable<? extends TNode> getNodes();

	/**
	 * @return the newick representation of this tree in a string.
	 */
	public String toNewick();

	/**
	 * @return the newick representation of this tree with the data stored for nodes in a string.
	 */
	public String toNewickWD();

	/**
	 * @return the newick representation of this tree with the data stored for nodes in a string.
	 */
	public String toStringWD();

	/**
	 * @return the list of tree nodes visited in the post order.
	 */
	public Iterable<TNode> postTraverse();

	/**
	 * Returns a string representation of this tree.  Valid arguments are {@link NEWICK_FORMAT}.
	 */
	public String toString(int format);

	/**
	 * @return the clusters of a given tree
	 */
	public List<STITreeCluster> getClusters(String leaves[], boolean gen);

	/**
	 * @return the bipartitions of a given tree
	 */
	public List<STITreeCluster> getBipartitionClusters(String leaves[], boolean gen);

	/**
	 * @return the clusters with data of a given tree
	 */
	public List<STITreeClusterWD> getClustersWD(String leaves[], boolean gen);

	/**
	 * @return all possible rooting trees
	 */
	public List<Tree> getAllRootingTrees();

	/**
	 * @return the gsi value of a given tree
	 */
	public double gsi(String leaves[]);

	/**
	 * Root a tree at an edge.
	 */
	public void rerootTreeAtEdge(int nodeID);

	/**
	 * Root a tree at a node.
	 */
	public void rerootTreeAtEdge(String nodeName);

	/**
	 * Root a tree at a node.
	 */
	public void rerootTreeAtNode(TNode node);
}
