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

package edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io;

import java.io.IOException;
import java.io.Reader;
import java.io.StreamTokenizer;
import java.util.LinkedList;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;



/**
 * <P>This class reads trees in newick format from a stream provided as a {@link Reader}.  The tree read is loaded
 * into a Tree object provided.</P>
 * 
 * <P>The exact format supported is as follows:</P>
 * 
 * 	[TREENAME =] [ROOTING] NEWICKSTRING[;]
 * 
 * <P>All components in square braces are optional.  TREENAME can be any legitimate string of characters.  ROOTING is
 * "[&U]" for unrooted trees and "[&R]" for rooted trees.  The semicolon is optionally used to indicate the end of the
 * tree.</P>
 * 
 * <P>Because of the way that this reader tokenizes, you cannot use a reader to read a tree and then use the same reader
 * to read other parts of a larger input stream.  Rather, it is necessary for you to read the tree containing section
 * out of the larger input stream and feed this string/stream to the tree by itself.  If you don't do this, the
 * tokenizer inside the NewickReader will potentially tokenize characters past the end of the trees.</P>
 * 
 * @author Derek Ruths
 */
public class NewickReader {

	// rooted constants
	private static final int ROOTED = 0;
	private static final int UNROOTED = 1;
	private static final int UNKNOWN = 2;
	
	// parsing constants
	private static final char LBRACE = '[';
	private static final char RBRACE = ']';
	private static final char AMPERSTAND = '&';
	private static final char COLON = ':';
	private static final char SEMICOLON = ';';
	private static final char COMMA = ',';
	private static final char LPAREN = '(';
	private static final char RPAREN = ')';
	private static final char EQUALS = '=';

	// fields
	private LookAheadStreamTokenizer _stok;
	
	private int _next_tree_rooted = UNKNOWN;
	
	// constructors
	public NewickReader(Reader r) {
		// configure the tokenizer
		StreamTokenizer stok = new StreamTokenizer(r);
		stok.resetSyntax();
		
		stok.wordChars('0','9');
		stok.wordChars('A','z');
		stok.wordChars('-','.');
		
		stok.whitespaceChars('\t',' ');
		
		stok.quoteChar('"');
		//_stok.parseNumbers();
		
		stok.ordinaryChar(COLON);
		stok.ordinaryChar(SEMICOLON);
		stok.ordinaryChar(COMMA);
		stok.ordinaryChar(LPAREN);
		stok.ordinaryChar(RPAREN);
		stok.ordinaryChar(LBRACE);
		stok.ordinaryChar(RBRACE);
		stok.ordinaryChar(AMPERSTAND);
		stok.ordinaryChar(EQUALS);
		
		_stok = new LookAheadStreamTokenizer(stok);
	}
	
	// methods
	/**
	 * This method reads in a tree and returns the tree read in.  It also handles the
	 * construction of a tree, using the {@link STITree} implementation of the {@link MutableTree}
	 * interface.  {@link readTree(Tree)} can be used to read the tree structure into
	 * a different tree implementation.
	 */
	public MutableTree readTree() throws IOException, ParseException {

		STITree<Object> tree = new STITree<Object>(isNextRooted());

		readTree(tree);

		return tree;
	}

	public String readTreeName() throws IOException, ParseException {

		// read the name
		if(_stok.nextToken() != StreamTokenizer.TT_WORD) {
			_stok.pushBack();
			return Tree.NO_NAME;
		}

		String name = _stok.sval;

		// try to read an equals sign
		if(_stok.nextToken() != EQUALS) {
			_stok.pushBack();
			_stok.pushBack();
			return Tree.NO_NAME;
		}

		return name;
	}

	/**
	 * @return <code>true</code> if the tree that will be parsed by the next call
	 * to {@link readTree} is rooted.
	 *
	 * @throws ParseException if the rootedness of this tree cannot be determined due to
	 * a syntax error.  Once this exception is thrown, this object may not parse remaining
	 * trees correctly.
	 */
	public boolean isNextRooted() throws IOException, ParseException {

		if(_next_tree_rooted == UNKNOWN) {
			_next_tree_rooted = (readRooted() == true)?ROOTED:UNROOTED;
		}

		return (_next_tree_rooted == ROOTED);
	}

	private boolean readRooted() throws IOException, ParseException {

		int token = _stok.nextToken();
		boolean is_rooted = true;

		// if unspecified, the tree is rooted
		if(token != LBRACE) {
			_stok.pushBack();
			is_rooted = true;
		} else {
			token = _stok.nextToken();

			if(token != AMPERSTAND) {

				// consume characters until we reach another brace
				while(token != RBRACE && token != StreamTokenizer.TT_EOF) {
					token = _stok.nextToken();
				}

				throw new ParseException("& character expected");
			}

			token = _stok.nextToken();

			if(token != StreamTokenizer.TT_WORD) {

				// consume characters until we reach another brace
				while(token != RBRACE && token != StreamTokenizer.TT_EOF) {
					token = _stok.nextToken();
				}

				throw new ParseException("R or U character expected");
			}

			if(_stok.sval.toUpperCase().equals("U")) {
				is_rooted = false;
			} else if(!_stok.sval.toUpperCase().equals("R")) {
				// consume characters until we reach another brace
				while(token != RBRACE && token != StreamTokenizer.TT_EOF) {
					token = _stok.nextToken();
				}

				throw new ParseException("R or U character expected");
			}

			// read the right brace
			token = _stok.nextToken();

			if(token != RBRACE) {

				// consume characters until we reach another brace
				while(token != RBRACE && token != StreamTokenizer.TT_EOF) {
					token = _stok.nextToken();
				}

				throw new ParseException("] character expected");
			}
		}

		return is_rooted;
	}

	/**
	 * This method reads a tree from the reader. The newick format is extended to include both bootstrap values
	 * and branch lengths.
	 *
	 * @param tree
	 * @throws IOException
	 * @throws ParseException
	 */
	public void readTree(STITree<Double> tree) throws IOException, ParseException {

		// get the tree name
		String tree_name = readTreeName();
		tree.setName(tree_name);

		// read the rooted-ness of the tree, if necessary
		readRooted();

		// clear information about the next tree
		_next_tree_rooted = UNKNOWN;

		// we'll need a root to start with
		if(tree.isEmpty()) {
			tree.createRoot();
		}

		readNode(tree.getRoot());

		int token = _stok.nextToken();

		// optionally read the semicolon
		if(token != SEMICOLON) {
			_stok.pushBack();
		}

		return;
	}

	/**
	 * This method reads a node in the extended newick format (with both bootstraps and branch lengths).
	 *
	 * @param node
	 * @throws IOException
	 * @throws ParseException
	 */
	private void readNode(STINode<Double> node) throws IOException, ParseException {

		// read the next token
		int token = _stok.nextToken();

		// read left paren
		if(token == LPAREN) {

			// read children
			token = COMMA;
			while(token == COMMA) {
				STINode<Double> child = node.createChild();

				readNode(child);

				token = _stok.nextToken();
			}

			if(token != RPAREN) {
				throw new ParseException("')' expected");
			}

			token = _stok.nextToken();
		}

		// read name
		if(token == StreamTokenizer.TT_WORD) {

			if(node.getTree().getNode(_stok.sval) != null) {
				throw new ParseException("Duplicate node " + _stok.sval);
			}

			node.setName(_stok.sval);

			token = _stok.nextToken();
		}

		// read length and bootstrap values.
		if(token == COLON) {
			token = _stok.nextToken();

			if (token != COLON) {
				// Either we have only branch lengths or both branch lengths and bootstraps
				if(token != StreamTokenizer.TT_WORD) {
					throw new ParseException("Branch length expected");
				}

				double dist = 0;

				try {
					dist = Double.parseDouble(_stok.sval);
				} catch(NumberFormatException nfe) {
					throw new ParseException("Number expected, found " + _stok.sval);
				}

				node.setParentDistance(dist);

				// Check if we have the bootstrap value.
				token = _stok.nextToken();
				if (token == COLON) {	// yes
					token = _stok.nextToken();
					if (token != StreamTokenizer.TT_WORD) {
						throw new ParseException("Bootstrap value expected");
					}

					double bootstrap = 0;

					try {
						bootstrap = Double.parseDouble(_stok.sval);
					}
					catch (NumberFormatException nfe) {
						throw new ParseException("Number expected, found " + _stok.sval);
					}

					node.setData(bootstrap);
				}
				else {	// no
					node.setData(new Double(0));

					_stok.pushBack();
				}
			}
			else {
				// Only bootstrap is available.
				token = _stok.nextToken();
				if (token != StreamTokenizer.TT_WORD) {
					throw new ParseException("Bootstrap value expected");
				}

				double bootstrap = 0;

				try {
					bootstrap = Double.parseDouble(_stok.sval);
				}
				catch (NumberFormatException nfe) {
					throw new ParseException("Number expected, found " + _stok.sval);
				}

				node.setData(bootstrap);
			}
		}
		else {
			_stok.pushBack();
		}

		// done!
		return;
	}

	/**
	 * This method reads a tree from the reader.  If the tree cannot be read, an exception will be thrown.
	 * <code>tree</code> should be a newly created tree, with no other structure in it.  The structure of
	 * the tree that is read will be pushed directly into the tree data structure underneath the root.
	 */
	public void readTree(MutableTree tree) throws IOException, ParseException {

		// get the tree name
		String tree_name = readTreeName();
		tree.setName(tree_name);

		// read the rooted-ness of the tree, if necessary
		readRooted();

		// clear information about the next tree
		_next_tree_rooted = UNKNOWN;

		// we'll need a root to start with
		if(tree.isEmpty()) {
			tree.createRoot();
		}

		readNode(tree.getRoot());

		int token = _stok.nextToken();

		// optionally read the semicolon
		if(token != SEMICOLON) {
			_stok.pushBack();
		}

		return;
	}

	private void readNode(TMutableNode node) throws IOException, ParseException {

		// read the next token
		int token = _stok.nextToken();

		// read left paren
		if(token == LPAREN) {

			// read children
			token = COMMA;
			while(token == COMMA) {
				TMutableNode child = node.createChild();

				readNode(child);

				token = _stok.nextToken();
			}

			if(token != RPAREN) {
				throw new ParseException("')' expected");
			}

			token = _stok.nextToken();
		}

		// read name
		if(token == StreamTokenizer.TT_WORD) {

			if(node.getTree().getNode(_stok.sval) != null) {
				throw new ParseException("Duplicate node " + _stok.sval);
			}

			node.setName(_stok.sval);

			token = _stok.nextToken();
		}

		// read length and bootstrap values.
		if(token == COLON) {
			token = _stok.nextToken();

			if(token != StreamTokenizer.TT_WORD) {
				throw new ParseException("Number expected");
			}

			double dist = 0;

			try {
				dist = Double.parseDouble(_stok.sval);
			} catch(NumberFormatException nfe) {
				throw new ParseException("Number expected, found " + _stok.sval);
			}

			node.setParentDistance(dist);
		} else {
			_stok.pushBack();
		}

		// done!
		return;
	}

	/**
	 * Indicates whether the end of the reader has been reached.
	 */
	public boolean reachedEOF() throws IOException {
		int token = _stok.nextToken();

		boolean at_eof = token == StreamTokenizer.TT_EOF;

		_stok.pushBack();

		return at_eof;
	}
}

class LookAheadStreamTokenizer {

	protected StreamTokenizer _stok;

	protected String sval;
	protected double nval;

	private static final char QUOTE = '"';
	
	protected int _curr_idx = -1;
	protected LinkedList<Integer> _past_tokens = new LinkedList<Integer>();
	protected LinkedList<String> _past_svals = new LinkedList<String>();
	protected LinkedList<Double> _past_nvals = new LinkedList<Double>();

	// constructors
	public LookAheadStreamTokenizer(StreamTokenizer stok) {
		_stok = stok;
	}

	// methods
	public int nextToken() throws IOException {

		int token = -1;

		if(_curr_idx < 0) {
			token = _stok.nextToken();
			sval = _stok.sval;
			nval = _stok.nval;

			_past_tokens.addFirst(token);
			_past_svals.addFirst(sval);
			_past_nvals.addFirst(nval);
		} else {
			sval = _past_svals.get(_curr_idx);
			nval = _past_nvals.get(_curr_idx);

			token = _past_tokens.get(_curr_idx);

			_curr_idx--;
		}
		
		if (token == QUOTE) {
		    return StreamTokenizer.TT_WORD;
		}
		return token;
	}

	public void pushBack() {
		_curr_idx++;
	}
}

