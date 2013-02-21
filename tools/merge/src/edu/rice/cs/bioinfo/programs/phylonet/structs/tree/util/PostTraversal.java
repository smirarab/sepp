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

import java.util.Iterator;
import java.util.Stack;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

/**
 * An iterable class that allows us to visit the tree in the post order. To make the class more 
 * general, we pass it a start node, so that it will visit all nodes under the start node in the 
 * post order. If we want to visit all nodes in a tree, simply pass it the tree's root.
 */
public class PostTraversal<D> implements Iterable<TNode> {
	private TNode _start_node;
		
	public PostTraversal(TNode start) {
		_start_node = start;
	}
	
	public Iterator<TNode> iterator() {
		return new TreePostOrderIterator();
	}
	
	private class TreePostOrderIterator implements Iterator<TNode> {
		private Stack<TNode> _unvisited;		// Stores unvisited nodes.
		private Stack<Iterator<? extends TNode>> _iters;	// Stores pointers to children to be visited.

		public TreePostOrderIterator() {
			_unvisited = new Stack<TNode>();
			_iters = new Stack<Iterator<? extends TNode>>();
			
			if (_start_node != null) {
				_unvisited.push(_start_node);
				if (!_start_node.isLeaf()) {	// i.e., it has children.
					_iters.push(_start_node.getChildren().iterator());
				}
				else {
					_iters.push(null);
				}
			}
		}
		
		public boolean hasNext() {
			return !_unvisited.isEmpty();
		}
		
		public TNode next() {
			assert(!_unvisited.isEmpty());
			
			while (true) {
				Iterator<? extends TNode> it = _iters.peek();
				
				if (it != null && it.hasNext()) {	// Still has more children to visit.
					TNode child = it.next();
					
					_unvisited.push(child);
					if (!child.isLeaf()) {
						_iters.push(child.getChildren().iterator());
					}
					else {
						_iters.push(null);
					}					
				}
				else {
					break;	// Found next node in post order traveral.
				}
			}
			
			_iters.pop();
			return _unvisited.pop();
		}
		
		public void remove() {
			System.err.println("This method is currently not supported.");
			return;
		}
	}
}