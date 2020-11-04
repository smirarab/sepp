# uym2 added
# June 2017
# utils for tree decomposition


from dendropy import Tree
try:
    from queue import Queue  # python 3
except ImportError:
    from Queue import Queue  # python 2
# from tree import PhylogeneticTree
from sepp import get_logger

_LOG = get_logger(__name__)


def decompose_by_diameter(a_tree, strategy, max_size=None, min_size=None,
                          max_diam=None):
    def __ini_record__():
        for node in a_tree.postorder_node_iter():
            __update_node__(node)

    def __find_midpoint_edge__(tre):
        u = tre.seed_node.bestLCA.anchor
        uel = u.edge_length if u.edge_length else 0
        d = 0
        while d + uel < tre.seed_node.diameter / 2:
            d += uel
            u = u.parent_node
            uel = u.edge_length if u.edge_length else 0
        return u.edge

    def __find_centroid_edge__(tre):
        u = tre.seed_node
        product = 0
        acc_nleaf = 0

        while not u.is_leaf():
            max_child = None
            max_child_nleaf = 0
            for ch in u.child_node_iter():
                if ch.nleaf > max_child_nleaf:
                    max_child_nleaf = ch.nleaf
                    max_child = ch
            acc_nleaf += (u.nleaf-max_child.nleaf)
            new_product = max_child.nleaf * acc_nleaf
            if new_product <= product:
                break
            product = new_product
            u = max_child

        return u.edge

    def __bisect__(tre, edg):
        # e = __find_centroid_edge__(t)

        u = edg.tail_node
        v = edg.head_node

        u.remove_child(v)
        tr1 = Tree(seed_node=v)

        if u.num_child_nodes() == 1:
            p = u.parent_node
            v = u.child_nodes()[0]
            l_v = v.edge_length if v.edge_length else 0
            u.remove_child(v)
            # u is the seed_node; this means the tree runs out of all but one
            # side
            if p is None:
                tre.seed_node = v
                return tre, tr1
            l_u = u.edge_length if u.edge_length else 0
            p.remove_child(u)
            p.add_child(v)
            v.edge_length = l_u + l_v
            u = p

        while u is not None:
            __update_node__(u)
            u = u.parent_node

        return tre, tr1

    def __clean_up__(tre):
        for node in tre.postorder_node_iter():
            delattr(node, "nleaf")
            delattr(node, "anchor")
            # delattr(node,"maxheight")
            delattr(node, "maxdepth")
            delattr(node, "diameter")
            # delattr(node,"topo_diam")
            delattr(node, "bestLCA")

    def __update_node__(node):
        if node.is_leaf():
            node.anchor = node
            # node.maxheight = 0
            node.maxdepth = 0
            node.diameter = 0
            # node.topo_diam = 0
            node.bestLCA = node
            node.nleaf = 1
            return

        # n1 = -1
        # n2 = -1
        d1 = -1
        d2 = -1
        anchor1 = None
        node.diameter = 0
        # node.topo_diam = 0
        node.bestLCA = None
        node.nleaf = 0

        for ch in node.child_node_iter():
            node.nleaf += ch.nleaf
#               n = ch.maxheight + 1
            d = ch.maxdepth + ch.edge_length if ch.edge_length else 0
#               if n > n1:
#                   n2 = n1
#                   n1 = n
#                   anchor2 = anchor1
#                   anchor1 = ch.anchor
#               elif n > n2:
#                   n2 = n
#                   anchor2 = ch.anchor
            if d > d1:
                d2 = d1
                d1 = d
                anchor1 = ch.anchor
            elif d > d2:
                d2 = d
            if ch.diameter > node.diameter:
                node.diameter = ch.diameter
                node.bestLCA = ch.bestLCA
#               node.diameter = max(ch.diameter,node.diameter)

#        node.diameter = max(d1+d2, node.diameter)
        node.maxdepth = d1
#        node.maxheight = n1
        node.anchor = anchor1
        if d1+d2 > node.diameter:
            node.diameter = d1+d2
            node.bestLCA = node

    def __get_breaking_edge__(tre, edge_type):
        if tre.seed_node.nleaf <= max_size and \
                tre.seed_node.diameter <= max_diam:
            return None
        if edge_type == 'midpoint':
            ed = __find_midpoint_edge__(tre)
        elif edge_type == 'centroid':
            ed = __find_centroid_edge__(tre)
        else:
            _LOG.warning(("Invalid decomposition type! Please use either "
                          "'midpoint' or 'centroid'"))
            return None

        n = ed.head_node.nleaf
        if (n < min_size) or (tre.seed_node.nleaf - n) < min_size:
            return None
        return ed

    def __check_stop__(tre):
        return ((tre.seed_node.nleaf <= max_size and
                 tre.seed_node.diameter <= max_diam) or
                (tre.seed_node.nleaf // 2 < min_size))

    def __break_by_MP_centroid__(tre):
        ed = __get_breaking_edge__(tre, 'midpoint')
        if ed is None:
            # print("Midpoint failed. Trying centroid decomposition...")
            ed = __get_breaking_edge__(tre, 'centroid')
        # else:
        #    print("Successfully splitted by midpoint")
        return ed

    def __break(tre):
        if strategy == "centroid":
            return __get_breaking_edge__(tre, 'centroid')
        elif strategy == "midpoint":
            return __break_by_MP_centroid__(tre)
        else:
            raise Exception("strategy not valid: %s" % strategy)

    tqueue = Queue()

    _LOG.debug("Starting brlen decomposition ...")
    __ini_record__()
    min_size = min_size if min_size else 0
    max_size = max_size if max_size else a_tree.seed_node.nleaf
    max_diam = max_diam if max_diam else a_tree.seed_node.diameter

    _LOG.debug(
        "Now breaking by %s with min %d and max %d sizes and diameter %f ..." %
        (strategy, min_size, max_size, max_diam))
    # try using midpoint
    e = __break(a_tree)

    if e is None:
        __clean_up__(a_tree)
        return [a_tree]

    tree_map = []
    tqueue.put((a_tree, e))
    while not tqueue.empty():
        t, e = tqueue.get()
        t1, t2 = __bisect__(t, e)
        e1 = __break(t1)
        if e1 is None:
            __clean_up__(t1)
            tree_map.append(t1)
        else:
            tqueue.put((t1, e1))
        e2 = __break(t2)
        if e2 is None:
            __clean_up__(t2)
            tree_map.append(t2)
        else:
            tqueue.put((t2, e2))

    return tree_map
