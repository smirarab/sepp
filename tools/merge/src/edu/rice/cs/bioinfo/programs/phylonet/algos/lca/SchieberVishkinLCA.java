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

package edu.rice.cs.bioinfo.programs.phylonet.algos.lca;


import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.Stack;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.util.LIDAMath;


/**
 * This class implements the LCA algorithm first described by Harel and Tarjan in their 1984 paper <B>Fast algorithms for finding nearest common ancestors</B> and then
 * simplified by Schieber and Vishkin in 1988 in their paper <B>On finding lowest common ancestors: simplications and parallelization</B>.  This approach provides
 * an <code>O(1)</code> lookup time.
 * 
 * @author Derek Ruths
 */
public class SchieberVishkinLCA {
	
	// constants
	private int LAST_ONE = 0x80000000;
	private int ZERO = 0x00000000;
	private int UNSIGNED_INT_SIZE = 32;

	// inner classes
	private class TreeState {
		public TNode node;
		public Iterator<? extends TNode> child_iterator;

		public TreeState(TNode node, Iterator<? extends TNode> child_iterator) {
			this.node = node;
			this.child_iterator = child_iterator;
		}
	};
	
	// fields
	private Tree _tree;

	int _num_nodes;
	int _log_num_nodes;

	TNode[] _back_ref;
	int[] _dfs_num;
	int[] _dfs_num_exit;
	int[] _h_val;
	int[] _I_val;
	int[] _run_leader;
	int[] _Av_val;
	
	private Hashtable<TNode,Integer> _tnode_idx_lookup;

	// constructor
	/**
	 * Construct an object that will return the LCA for nodes in the specified tree.
	 * 
	 * @param tree is a rooted tree.
	 * 
	 * @throws RuntimeException if the tree is unrooted
	 */
	public SchieberVishkinLCA(Tree tree) {
		
		if(tree.isRooted() == false) {
			throw new RuntimeException("Tree must be rooted");
		}
		
		_tree = tree;
		
		preprocess();
	}
	
	// helper methods
	private int calc_h(int val) {
	  int pos = 1;

	  while((val & 1) == 0) {
	    val = val >> 1;
	    pos++;
	  } 

	  return pos;
	}

	private int calc_Iw(int j, int Ax, int Ix) {

	  int k = j - 1;  
	  int tmp = Ax;
	  tmp <<= (UNSIGNED_INT_SIZE - k);

	  while((tmp & LAST_ONE) == 0) {
	    tmp <<= 1;
	    k--;
	  }

	  int Iw;
	  Iw = (((int)Ix) >> k) << k;
	  Iw = Iw | (1 << (k - 1));
	  return (int) Iw;
	}

	private int calc_hIz(int b, int Az1, int Az2) {
	  int hb = calc_h(b);
	  int hIz = hb; 

	  Az1 >>= (hb - 1);
	  Az2 >>= (hb - 1);

	  while(!((Az1 & 1) != ZERO && (Az2 & 1) != ZERO)) {
	    Az1 >>= 1;
	    Az2 >>= 1;

	    hIz++;
	  }

	  return hIz;
	}

	private boolean is_binary_root(int n1, int log_n) {
	  
	  int p = (n1 >> (log_n - 1)) << (log_n - 1);

	  return (p == n1 && (n1 != 0));
	}

	int calc_lca_path_num(int n1, int n2, int log_n) {

	  if(is_binary_root(n1, log_n)) {
	    return n1;
	  } else if(is_binary_root(n2, log_n)) {
	    return n2;
	  }

	  int p = n1 ^ n2;

	  // find the left most one position
	  int lm_one = 0;

	  if(p == 0) {
	    lm_one = -1;
	  } else {
	    while((p & LAST_ONE) == 0) {
	      p = p << 1;
	      lm_one++;
	    }
	  }

	  // tack on another '1' in front of the left most 1 and zero it out
	  // before the one
	  p = n1 | n2;

	  int r = (UNSIGNED_INT_SIZE) - lm_one - 1;

	  p = (p >> r) << r;
	  p = p & 0xFFFFFFFE;

	  return p;
	}

	// give each node an index number
	private void numberTNodes() {
		
		_tnode_idx_lookup = new Hashtable<TNode,Integer>();
		
		int i = 0;
		for(TNode n : _tree.getNodes()) {
			_tnode_idx_lookup.put(n, new Integer(i++));
		}
		
		return;
	}
	
	private int getIndex(TNode n) {
		return _tnode_idx_lookup.get(n);
	}
	
	// main algorithm methods
	private void preprocess() {
		// index the nodes
		numberTNodes();
		
		// get the number of nodes in the tree
		_num_nodes = _tree.getNodeCount();
		_log_num_nodes = (int) Math.ceil(LIDAMath.log(2.0, (double) _num_nodes));

		// allocate a bunch of arrays
		_back_ref = new TMutableNode[_num_nodes];

		_dfs_num = new int[_num_nodes];
		_dfs_num_exit = new int[_num_nodes];
		_h_val = new int[_num_nodes];
		_I_val = new int[_num_nodes];
		_run_leader = new int[_num_nodes + 1]; // indexed by DFS node numbers
		_Av_val = new int[_num_nodes];

		// some useful variables
		//LinkedList stack = new LinkedList();
		Stack<TreeState> stack = new Stack<TreeState>();
		TMutableNode child;
		TreeState ts;
		
		// Step 1: assign depth-first search numbers and back references
		int dfs_idx = 1;

		/*
		_dfs_num[_tree->Get_Index()] = dfs_idx;
		_back_ref[_tree->Get_Index()] = _tree;
		_h_val[_tree->Get_Index()] = calc_h(dfs_idx);
		dfs_idx++;
		stack.Push(new TreeState(_tree, 0));
		 */
		_dfs_num[getIndex(_tree.getRoot())] = dfs_idx;
		_back_ref[getIndex(_tree.getRoot())] = _tree.getRoot();
		_h_val[getIndex(_tree.getRoot())] = calc_h(dfs_idx);
		dfs_idx++;
		stack.push(new TreeState(_tree.getRoot(), _tree.getRoot().getChildren().iterator()));
		
		while(!stack.isEmpty()) {
			ts = (TreeState) stack.peek();

			if(ts.child_iterator.hasNext()) {
				child = (TMutableNode) ts.child_iterator.next();
				
				int idx = getIndex(child);
				_dfs_num[idx] = dfs_idx;
				_back_ref[idx] = child;
				_h_val[idx] = calc_h(dfs_idx);

				dfs_idx++;

				stack.push(new TreeState(child, child.getChildren().iterator()));
			} else {
				ts = (TreeState) stack.pop();

				// record the last dfs value that was assigned under/in this
				// node
				_dfs_num_exit[getIndex(ts.node)] = dfs_idx - 1;
			}
		}

		// debugging
		/*
		if(_debug & LCA_DEBUG_DFS) {
			cout << "DFS\n";
			for(int i = 0; i < _num_nodes; i++) {
				cout << "\t" << _back_ref[i]->Get_Name() << " " << _dfs_num[i] << "\n";
			}
		}
		
		if(_debug & LCA_DEBUG_H) {
			cout << "H\n";
			for(int i = 0; i < _num_nodes; i++) {
				cout << "\t" << _back_ref[i]->Get_Name() << " " << _h_val[i] << "\n";
			}
		}
		 */
		
		// Step 2: calculate I(v) and the leader of runs
		stack.clear();

		stack.push(new TreeState(_tree.getRoot(), _tree.getRoot().getChildren().iterator()));

		while(!stack.isEmpty()) {
			ts = (TreeState) stack.peek();
	    
			if(ts.child_iterator.hasNext()) {
				child = (TMutableNode) ts.child_iterator.next();

				stack.push(new TreeState(child, child.getChildren().iterator()));
			} else {
				updateI(ts.node);
				ts = (TreeState) stack.pop();
			}
		}

		// debugging
		/*
		if(_debug & LCA_DEBUG_I) {
			cout << "I\n";
			for(int i = 0; i < _num_nodes; i++) {
				cout << "\t" << _back_ref[i]->Get_Name() << " " << _I_val[i] << "\n";
			}
		}
		
		if(_debug & LCA_DEBUG_L) {
			cout << "L\n";
			for(int i = 0; i < _num_nodes; i++) {
				cout << "\t" << _back_ref[i]->Get_Name() << " ";
				cout << _run_leader[_dfs_num[_I_val[i]]] << "\n";
			}
		}
		 */
		
		// Step 3: Calculate Av values
		int tmp_av;

		stack.clear();

		tmp_av = 1;
		tmp_av = tmp_av << (calc_h(_dfs_num[_I_val[getIndex(_tree.getRoot())]]) - 1);
		_Av_val[getIndex(_tree.getRoot())] = tmp_av;

		stack.push(new TreeState(_tree.getRoot(), _tree.getRoot().getChildren().iterator()));

		while(!stack.isEmpty()) {
			ts = (TreeState) stack.peek();

			if(ts.child_iterator.hasNext()) {
				child = (TMutableNode) ts.child_iterator.next();

				tmp_av = 1;
				tmp_av = tmp_av << (calc_h(_dfs_num[_I_val[getIndex(child)]]) - 1);
				tmp_av = tmp_av | _Av_val[getIndex(child.getParent())];
				_Av_val[getIndex(child)] = tmp_av;

				stack.push(new TreeState(child, child.getChildren().iterator()));
			} else {
				ts = (TreeState) stack.pop();
			}
		}

		/*
		// debugging
		if(_debug & LCA_DEBUG_A) {
			cout << "A\n";
			for(int i = 0; i < _num_nodes; i++) {
				cout << "\t" << _back_ref[i]->Get_Name() << " " << _Av_val[i] << "\n";
			}
		}
		 */
		
		// DONE!
		return;
	}
	 
	private void updateI(TNode t) {
	  
		if(t.isLeaf()) {
			_I_val[getIndex(t)] = getIndex(t);
			_run_leader[_dfs_num[_I_val[getIndex(t)]]] = getIndex(t);
		} else {
			int self_h = _h_val[getIndex(t)];
			int max_h = self_h;
			int idx = -1;

			// set I
			for(Iterator i = t.getChildren().iterator(); i.hasNext(); ) {
				TMutableNode child = (TMutableNode) i.next();

				if(_h_val[_I_val[getIndex(child)]] > max_h) {
					max_h = _h_val[_I_val[getIndex(child)]];
					idx = _I_val[getIndex(child)];
				}
			}

			if(max_h > self_h) {
				_I_val[getIndex(t)] = idx;
			} else {
				_I_val[getIndex(t)] = getIndex(t);
			}

			// set any necessary leaders
			int self_I = _I_val[getIndex(t)];

			for(TNode child : t.getChildren()) {
				if(_I_val[getIndex(child)] != self_I) {
					_run_leader[_dfs_num[_I_val[getIndex(child)]]] = getIndex(child);
				} else {
					_run_leader[_dfs_num[_I_val[getIndex(t)]]] = getIndex(t);
				}
			}
		}
	}

	private boolean isAncestor(TNode n1, TNode n2) {

		int dfs2 = _dfs_num[getIndex(n2)];
		int dfs1 = _dfs_num[getIndex(n1)];
		int dfs1_exit = _dfs_num_exit[getIndex(n1)];

		if(dfs1 > dfs2) {
			return false;
		}

		if(dfs2 <= dfs1_exit) {
			return true;
		}

		return false;
	}

	/**
	 * Find the LCA for a pair of nodes in the tree used to construct this object.
	 * 
	 * @param n1 is the first node
	 * @param n2 is the second node
	 * 
	 * @return the LCA of <code>n1</code> and <code>n2</code>
	 */
	public TNode getLCA(TNode n1, TNode n2) {

		// handle the trivial case
		if(n1 == n2) {
			return n1;
		}

		// if one node is an ancestor of another node
		if(isAncestor(n1, n2)) {
			return n1;
		} else if(isAncestor(n2, n1)) {
			return n2;
		}
	  
		// find the lca in conceptual binary tree for nodes I(x) and I(y)
		int b;
		int hIz;

		b = calc_lca_path_num(_dfs_num[_I_val[getIndex(n1)]], 
							 _dfs_num[_I_val[getIndex(n2)]], _log_num_nodes);

		/*
		if(_debug & LCA_DEBUG_B) {
			cout << "B " << b << "\n";
		}
		 */
		
		// Find h(I(z))
		hIz = calc_hIz(b, _Av_val[getIndex(n1)], 
					  _Av_val[getIndex(n2)]);

		/*
		if(_debug & LCA_DEBUG_HIZ) {
			cout << "hIz " << hIz << "\n";
		}
	     */
		
		// Find x_bar = closest node to x on the same run as z
		int x_bar = findNodeBar(n1, hIz);

		// Find y_bar
		int y_bar = findNodeBar(n2, hIz);

		// handle error case
		if(x_bar == -1 || y_bar == -1) {
			return null;
		}

		return _back_ref[((_dfs_num[x_bar] < _dfs_num[y_bar]) ? x_bar : y_bar)];
	}

	/**
	 * Find the LCA for a (non-empty) set of nodes.
	 * 
	 * @param nodes is the set of nodes.  This set must not be empty.
	 * 
	 * @return the LCA calculated for this set of nodes.
	 * 
	 * @throws RuntimeException if the set of nodes is empty.
	 */
	public TNode getLCA(Set<? extends TNode> nodes) {
		
		if(nodes.isEmpty()) {
			throw new RuntimeException("The set is empty");
		}
		
		Iterator<? extends TNode> i = nodes.iterator();
		
		TNode lca = i.next();
		
		while(i.hasNext()) {
			lca = getLCA(lca, i.next());
		}
		
		return lca;
	}
	
	private int findNodeBar(TNode t, int hIz) {
		
		int n_bar;
		int l = calc_h(_Av_val[getIndex(t)]);
		
		if(l == hIz) {
			n_bar = getIndex(t);
		} else if(l < hIz) {
			int I_w = calc_Iw(hIz, _Av_val[getIndex(t)], 
					_dfs_num[_I_val[getIndex(t)]]);
			
			if(_back_ref[_run_leader[I_w]] == _tree.getRoot()) {
				n_bar = getIndex(_tree.getRoot());
			} else {
				n_bar = getIndex(_back_ref[_run_leader[I_w]].getParent());
			}
		} else { // This should never happen
			System.err.println("FATAL ERROR: l > hIz in LCA:Get_LCA(...,...)");
			return -1;
		}
		
		return n_bar;
	}
	
	public Tree getTree() {
		return _tree;
	}
	
	// command line version
	public static void printUsage() {
		System.err.println();
		System.err.println("This is an implementation of the LCA algorithm described by");
		System.err.println("Schieber and Vishkin.  By default the tree is read from the");
		System.err.println("stdin and then the node list is read from the stdin.  The");
		System.err.println("tree should be rooted and specified in newick format. The");
		System.err.println("node list consists of the names of valid nodes in the tree.");
		System.err.println("Each line of the node list should contain a list of space");
		System.err.println("separated names of nodes in the tree.  The LCA will be");
		System.err.println("computed for each line.  The following properties can");
		System.err.println("be specified on the command line:");
		System.err.println();
		System.err.println("\t-h\tPrints this message");
		System.err.println("\t-o\tSpecifies a file that the output should be written to");
		System.err.println("\t-t\tSpecifies a file that contains the tree");
		System.err.println("\t-n\tSpecifies a file that contains the node list");
		System.err.println();
		System.err.println("The first line of output is the original tree with all");
		System.err.println("Internal nodes labeled.  Following this is the list of");
		System.err.println("computed LCAs.");
		System.err.println();
	}
	    /*
	public static void main(String[] args) {
		
		InputStream node_list_stream = System.in;
		PrintStream output = System.out;
		MutableTree tree = null;
		
		// read arguments
		for(int idx = 0; idx < args.length; idx++) {
			if(args[idx].charAt(0) != '-') {
				break;
			}

			if(args[idx].equals("-h")) {
				printUsage();
				return;
			} else if(args[idx].equals("-o")) {
				if(idx + 1 == args.length) {
					System.err.println("ERROR: Missing output file");
					printUsage();
					return;
				}
				
				File f = new File(args[++idx]);
				
				try {
					output = new PrintStream(new FileOutputStream(f));
				} catch(FileNotFoundException fnfe) {
					System.err.println("ERROR: Invalid output file");
					return;
				}
			} else if(args[idx].equals("-t")) {
				if(idx + 1 == args.length) {
					System.err.println("ERROR: Input tree property requires argument");
					printUsage();
					return;
				}
				
				NewickReader nr;
				
				try {
					nr = new NewickReader(new FileReader(args[++idx]));
				} catch(FileNotFoundException fnfe) {
					System.err.println("ERROR: Invalid input file " + args[idx]);
					return;
				}

				try {
					tree = nr.readTree();
				} catch(ParseException ioe) {
					System.err.println("ERROR: Invalid tree in file " + args[idx]);
					return;
				} catch(IOException ioe) {
					System.err.println("ERROR: Unable to read file " + args[idx]);
					return;
				}
			} else if(args[idx].equals("-n")) {
				if(idx + 1 == args.length) {
					System.err.println("ERROR: Node list property requires argument");
					printUsage();
					return;
				}
				
				try {
					node_list_stream = new FileInputStream(args[++idx]);
				} catch(FileNotFoundException fnfe) {
					System.err.println("ERROR: Unable to find file " + args[idx]);
					printUsage();
					return;
				}
			} else {
				System.err.println("ERROR: Unrecognized command-line argument " + args[idx]);
				printUsage();
				return;
			}
		}
		
		if(tree == null) {
			System.out.print("Enter rooted tree> ");
		
			NewickReader nr = new NewickReader(new InputStreamReader(System.in));
			try {
				tree = nr.readTree();
			} catch(IOException ioe) {
				System.err.println("ERROR: Unable to read file from stdin");
				return;
			} catch(ParseException pe) {
				System.err.println("ERROR: Invalid tree");
				return;
			}
		}
		
		if(tree.isRooted() == false) {
			System.err.println("ERROR: Tree must be rooted");
			return;
		}
		
		// auto label the tree
		Trees.autoLabelNodes(tree);
		
		// print the tree out
		NewickWriter nw = new NewickWriter(new OutputStreamWriter(output));
		nw.writeTree(tree, true);
		
		// compute LCAs!
		SchieberVishkinLCA lca_calculator = new SchieberVishkinLCA(tree);
		
		// read in a line
		LineNumberReader lr = new LineNumberReader(new InputStreamReader(node_list_stream));
		
		while(true) {
			if(node_list_stream == System.in) {
				System.out.print("Enter list of nodes> ");
			}
			
			String line;
			
			try {
				line = lr.readLine();
			} catch(IOException ioe) {
				System.err.println("ERROR: Unable to read node list");
				return;
			}
			
			if(line == null) {
				break;
			}
			
			line = line.trim();
			
			if(line.length() == 0) {
				continue;
			}
			
			String[] node_names = line.split("\\s+");
			
			if(node_names.length < 1) {
				System.err.println("ERROR: At least one node must be specified");
				continue;
			}
			
			TNode lca = tree.getNode(node_names[0]);
			
			if(lca == null) {
				System.err.println("ERROR: Node name " + node_names[0] + " does not name a node in the input tree");
				continue;
			}
			
			boolean error = false;
			
			for(int i = 1; i < node_names.length; i++) {

				if(lca == null) {
					System.err.println("ERROR: Unknown error occured while computing LCA");
					error = true;
					break;
				}
				
				TNode node = tree.getNode(node_names[i]);
				
				if(node == null) {
					System.err.println("ERROR: Node name " + node_names[i] + " does not name a node in the input tree");
					error = true;
					break;
				}
				
				lca = lca_calculator.getLCA(lca, node);
			}
			
			if(error) {
				continue;
			}
			
			// print the lca
			output.println(lca.getName());
		}
	} */
}
