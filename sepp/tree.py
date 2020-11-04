"""SATe - Phylogenetic Tree Container, effectively a wrapper of
   dendropy.Tree"""

# This file is part of SATe

# SATe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Jiaye Yu and Mark Holder, University of Kansas

# This file is copied to SEPP and is used as an external library

from dendropy import Tree, Taxon, treecalc
from dendropy import DataSet as Dataset
from dendropy.datamodel.treemodel import _convert_node_to_root_polytomy as \
    convert_node_to_root_polytomy
from sepp import get_logger, sort_by_value
from sepp.alignment import get_pdistance
from sepp.decompose_tree import decompose_by_diameter


try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import copy
import string
import random
import re

_LOG = get_logger(__name__)


class PhylogeneticTree(object):
    """Data structure to store phylogenetic tree, wrapping dendropy.Tree."""
    def __init__(self, dendropy_tree, map_internal_node_names=True):
        self._tree = dendropy_tree
        assert isinstance(self._tree, Tree)
        self.n_leaves = self.count_leaves()
        self._tree.seed_node.edge.tail_node = None
        self._tree.seed_node.edge.length = None
        self._namemap = None
        self._revscript = None

        if map_internal_node_names:
            self.map_seq_names()

    def map_seq_names(self):
        i = 1
        tag = ''.join(random.choice(string.ascii_letters)
                      for letter in range(8))
        self._namemap = dict()
        revnamemap = dict()
        for n in self._tree.internal_nodes():
            if n.label:
                l1 = '%sN%06d%s' % (tag, i, re.sub(
                    r"[%s]" % string.punctuation.replace("_", " "),
                    repl="_", string=n.label))
                self._namemap[n] = l1
                revnamemap[l1] = n.label
                i = i + 1
        self._revscript = '''import ast, re, string, sys
revnamemap = ast.literal_eval("%s")
def relabel_newick(newick_string):
    pattern = re.compile("(%sN[^(,:)<>]+)")
    invalidChars = set(string.punctuation).union(set(string.whitespace))
    def replace_func(m):
        repl = m.group(1)
        if m.group(1) in revnamemap:
            repl = revnamemap[m.group(1)]
            if any(char in invalidChars for char in repl):
                repl = "'%%s'" %%repl
        else:
            repl = m.group(1)

        return repl
    t = pattern.sub(replace_func,newick_string)
    return t
for l2 in sys.stdin.readlines():
        sys.stdout.write(relabel_newick(l2))''' % (str(revnamemap), tag)
        # _LOG.info("Use this code to map back internal nodes::\n %s" %
        #  str(self._revscript))

    def rename_script(self):
        return self._revscript

    def write_newick_node(self, node, out):
        child_nodes = node.child_nodes()
        if child_nodes:
            out.write('(')
            f_child = child_nodes[0]
            for child in child_nodes:
                if child is not f_child:
                    out.write(',')
                self.write_newick_node(child, out)
            out.write(')')

        if node.is_leaf() or not self._namemap:
            out.write("%s" % node._get_node_token())
        elif node in self._namemap:
            out.write("%s" % self._namemap[node])
        e = node.edge
        if e:
            sel = e.length
            if sel is not None:
                s = ""
                try:
                    s = float(sel)
                    s = str(s)
                except ValueError:
                    s = str(sel)
                if s:
                    out.write(":%s[%s]" % (s, e.label))

    def get_tree(self):
        return self._tree
    den_tree = property(get_tree)

    def count_leaves(self):
        return len(self._tree.leaf_nodes())

    def count_nodes(self):
        return len(self._tree.nodes())

    def calc_splits(self):
        for i in self._tree.postorder_edge_iter():
            nd = i.head_node
            if nd.is_leaf():
                i.num_leaves_below = 1
            else:
                i.num_leaves_below = sum([j.edge.num_leaves_below
                                          for j in nd.child_nodes()])

    def get_clade_edge(self, minSize):
        root = self._tree.seed_node
        # should only be 2 children, but in unlikely event that clade-based
        # decomp results in 1 child, then we should handle in some way
        root_children = root.child_nodes()
        assert len(root_children) == 2

        # Always break first child edge
        clade_edge = root.child_nodes()[0].edge
        assert clade_edge is not None
        return clade_edge

    def get_centroid_edge(self, minSize):
        """Get centroid edge"""
        root = self._tree.seed_node
        root_children = root.child_nodes()
        if root_children and (not hasattr(root_children[0].edge,
                                          "num_leaves_below")):
            self.calc_splits()
            n_leaves = self.count_leaves()
        else:
            if root.edge:
                n_leaves = root.edge.num_leaves_below
            else:
                n_leaves = sum([j.edge.num_leaves_below
                                for j in root_children])
        centroid_edge = None
        centroid_imbalance = n_leaves
        if (n_leaves <= minSize):
            return None
        half_taxa = n_leaves / 2
        for edge in self._tree.postorder_edge_iter():
            if edge.tail_node is None:
                continue
            n_descendants = edge.num_leaves_below
            if n_descendants > 1:
                pass
            if (
                (minSize is not None) and
                ((n_descendants < minSize) or
                 (self.n_leaves - n_descendants < minSize))
               ):
                continue
            imbalance = abs(half_taxa - n_descendants)
            if (imbalance < centroid_imbalance):
                centroid_edge = edge
                centroid_imbalance = imbalance
            assert centroid_edge is not None
        return centroid_edge

    def get_longest_edge(self, minSize):
        longest_edge = None
        longest_len = -1.0
        for edge in self._tree.postorder_edge_iter():
            if edge.tail_node is None:
                continue
            onesideSize = len(edge.head_node.leaf_nodes())
            if (
                (minSize is not None) and
                ((onesideSize < minSize) or
                 (self.n_leaves - onesideSize < minSize))
               ):
                continue
            if edge.length is not None and edge.length > longest_len:
                longest_edge = edge
                longest_len = edge.length
        assert longest_edge is not None
        return longest_edge

    def get_adjacent_edges(self, e):
        he = [i for i in e.head_node.get_incident_edges() if i is not e]
        te = [i for i in e.tail_node.get_incident_edges() if i is not e]
        he.extend(te)
        return he

    def get_breaking_edge(self, option, minSize):
        if option.lower() == 'centroid':
            return self.get_centroid_edge(minSize)
        elif option.lower() == 'longest':
            return self.get_longest_edge(minSize)
        elif option.lower() == 'clade':
            return self.get_clade_edge(minSize)
        else:
            raise ValueError('Unknown break strategy "%s"' % option)

    def bipartition_by_root(self):
        if (self.n_leaves == 1):
            return (None, None, None)
        root = self._tree.seed_node
        t1_root = root._child_nodes[0]
        t = self._tree
        t.prune_subtree(t1_root, update_splits=True, delete_outdegree_one=True)
        t1 = PhylogeneticTree(t)
        t2 = PhylogeneticTree(Tree(t1_root))
        # Reroot if there's more than node left
        if (t2.n_leaves > 1):
            t2._tree.reroot_at_node(t1_root)
        return t1, t2, root

    def bipartition_by_edge(self, e):
        """Prunes the subtree that attached to the head_node of edge e and
           returns them as a separate tree."""

        t = self._tree
        nr = e.head_node
        assert e.tail_node is not None
        assert e.head_node is not None
        assert nr.parent_node is e.tail_node
        is_valid_tree(t)

        potentially_deleted_nd = e.tail_node
        grandparent_nd = potentially_deleted_nd.parent_node
        e.tail_node.remove_child(nr, suppress_unifurcations=True)

        nr.edge.length = None
        nr.parent_node = None
        convert_node_to_root_polytomy(nr)
        t1 = PhylogeneticTree(Tree(seed_node=nr))
        # temp we could speed this up,
        # by telling the Phylogenetic tree how many leaves it has
        n1 = t1.n_leaves

        if hasattr(e, "num_leaves_below"):
            if grandparent_nd is None:
                old_root = potentially_deleted_nd
                if old_root.edge:
                    old_root.edge.num_leaves_below -= n1
            else:
                if potentially_deleted_nd in grandparent_nd.child_nodes():
                    potentially_deleted_nd.edge.num_leaves_below -= n1
                old_root = grandparent_nd
                if old_root.edge:
                    old_root.edge.num_leaves_below -= n1
                while old_root.parent_node:
                    old_root = old_root.parent_node
                    if old_root.edge:
                        old_root.edge.num_leaves_below -= n1
        else:
            old_root = grandparent_nd or potentially_deleted_nd
            while old_root.parent_node:
                old_root = old_root.parent_node

        t2 = PhylogeneticTree(Tree(seed_node=old_root))

        is_valid_tree(t1._tree)
        is_valid_tree(t2._tree)
        return t1, t2

    def leaf_node_names(self):
        leaves = self._tree.leaf_nodes()
        return [i.taxon.label for i in leaves]

    def compose_newick(self, labels=False):
        if not labels:
            return self._tree.as_string(
                schema="newick", suppress_rooting=True,
                suppress_internal_node_labels=not labels)
        else:
            stringIO = StringIO()
            self.write_newick_node(self._tree.seed_node, stringIO)
            ret = stringIO.getvalue()
            stringIO.close()
            return ret

    def write_newick_to_path(self, path):
        tree_handle = open(path, "w")
        tree_handle.write(self.compose_newick())
        tree_handle.write("")
        tree_handle.close()

    def read_tree_from_file(self, treefile, file_format):
        dataset = Dataset()
        dataset.read(open(treefile, 'rU'), schema=file_format)
        dendropy_tree = dataset.trees_blocks[0][0]
        self._tree = dendropy_tree
        self.n_leaves = self.count_leaves()

    def get_subtree(self, taxa):
        if len(taxa) == 0:
            return None
        tree = Tree(self._tree)
        if isinstance(taxa[0], str):
            tree.prune_taxa_with_labels(taxa)
        elif isinstance(taxa[0], Taxon):
            tree.prune_taxa(taxa)
        return PhylogeneticTree(tree)

    def bisect_tree(self, breaking_edge_style='centroid', minSize=None):
        """Partition 'tree' into two parts
        """
        snl = self.n_leaves
        if (breaking_edge_style != 'clade'):
            e = self.get_breaking_edge(breaking_edge_style, minSize)
            if (e is None):
                return None, None, None
            _LOG.debug("breaking_edge length = %s, %s" % (
                e.length, breaking_edge_style))
            tree1, tree2 = self.bipartition_by_edge(e)
        else:
            tree1, tree2, e = self.bipartition_by_root()

        _LOG.debug("Tree 1 has %s nodes, tree 2 has %s nodes" % (
            tree1.n_leaves, tree2.n_leaves))
        assert snl == tree1.n_leaves + tree2.n_leaves
        return tree1, tree2, e

#    def decompose_tree(self, maxSize, strategy, minSize = None, tree_map={},
#                       decomp_strategy='normal', pdistance=1, distances=None):
    def decompose_tree(self, maxSize, strategy, minSize=None, tree_map={},
                       decomp_strategy='normal', pdistance=1, distances=None,
                       maxDiam=None):
        """
        This function decomposes the tree until all subtrees are smaller than
        the max size, but does not decompose below min size.
        Two possible decompositions strategies can used: "centroid" and
        "longest".
        Returns a map containing the subtrees, in an ordered fashion.

        SIDE EFFECT: deroots the tree (TODO: necessary?)
        """
        # uym2 added #
        if (decomp_strategy in ["midpoint", "centroid"]):
            T = decompose_by_diameter(
                self._tree, strategy=decomp_strategy, max_size=maxSize,
                max_diam=maxDiam, min_size=minSize)
            for i, t in enumerate(T):
                tree_map[i] = PhylogeneticTree(t)
            return tree_map
        ##############

        # Don't deroot if doing clade-based decomposition
        if (strategy != 'clade'):
            self._tree.deroot()
        else:
            # If doing clade-based decomp and it's not rooted, root it!
            if self._tree.is_rooted is False:
                self._tree.reroot_at_midpoint()

        if (
            (decomp_strategy == 'hierarchical') and
            (self.count_leaves() > maxSize)
           ):
            tree_map[len(tree_map)] = copy.deepcopy(self)
        if (
            (self.count_leaves() > maxSize) or
            ((pdistance != 1) and
             (get_pdistance(distances, self.leaf_node_names()) > pdistance))
           ):
            (t1, t2, e) = self.bisect_tree(strategy, minSize)
            if e is not None:
                t1.decompose_tree(maxSize, strategy, minSize, tree_map,
                                  decomp_strategy, pdistance, distances)
                t2.decompose_tree(maxSize, strategy, minSize, tree_map,
                                  decomp_strategy, pdistance, distances)
            else:
                tree_map[len(tree_map)] = self
                _LOG.warning(
                    ("It was not possible to break-down the following tree "
                     "according to given subset sizes: %d , %d:\n %s") % (
                     minSize, maxSize, self._tree))
        else:
            tree_map[len(tree_map)] = self
        return tree_map

    def lable_edges(self):
        en = 0
        for e in self._tree.postorder_edge_iter():
            e.label = en
            en += 1

    '''
    Returns a given number of taxa that are closest to a given leaf
    '''
    def branchOut(self, centerTaxon, subsetSize, **kwargs):
        dist = {}
        pdm = treecalc.PatristicDistanceMatrix(self.den_tree)
        for i, s in enumerate(self.den_tree.taxon_set):  # @UnusedVariable
            if "filterTaxon" in kwargs:
                if not kwargs["filterTaxon"](s):
                    continue
            dist[s.label] = pdm(centerTaxon, s)
        incircle = sort_by_value(dist)[0:subsetSize]
        return [node[0] for node in incircle]


def node_formatter(n):
    return str(id(n))


def edge_formatter(e):
    return "%s %f " % (str(id(e)), e.length)


def is_valid_tree(t):
    if (t.is_rooted):
        return True
    assert t and t
    rc = t.seed_node.child_nodes()
    num_children = len(rc)
    if num_children == 0:
        return True
    if num_children == 1:
        assert not rc[0].child_nodes()
        return True
    if num_children == 2:
        # What is with this code?  Why do we check the same variable twice?
        # Bug?  NN
        assert((not rc[0].child_nodes()) and (not rc[0].child_nodes()))
    return True
