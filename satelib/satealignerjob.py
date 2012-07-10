#!/usr/bin/env python

"""Main script of SATe in command-line mode
"""

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


import math
import os
import copy
from threading import Lock
from satelib import get_logger
_LOG = get_logger(__name__)

def bisect_tree(tree, breaking_edge_style='centroid'):
    """Partition 'tree' into two parts
    """
    e = tree.get_breaking_edge(breaking_edge_style)
    _LOG.debug("breaking_edge length = %s, %s" % (e.length, breaking_edge_style) )
    snl = tree.n_leaves
    tree1, tree2 = tree.bipartition_by_edge(e)
    _LOG.debug("Tree 1 has %s nodes, tree 2 has %s nodes" % (tree1.n_leaves, tree2.n_leaves) )
    assert snl == tree1.n_leaves + tree2.n_leaves
    return tree1, tree2