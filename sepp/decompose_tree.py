from dendropy import Tree, Node

def decompose_by_diameter(a_tree,min_diam=None,nsubtree=1024):
    def __ini_record__():
        self.bestLCA = None
        self.nleaf = 0

        a_tree.bestLCA = None

        for node in a_tree.postorder_node_iter():
               __updateNode__(node)
               if (a_tree.bestLCA is None or node.record[3] > a_tree.bestLCA.record[3]):
                   a_tree.bestLCA = node

        min_diam = min_diam if min_diam is not None else a_tree.seed_node.diameter/nsubtree

    def __find_centroid_edge__():
        h = 0
        u = a_tree.bestLCA.record[0]
        while (h < a_tree.bestLCA.record[3]):
            h += 1
            u = u.parent_node()
        return u.edge

    def __updateNode__(node):
        if node.is_leaf():
            node.record = [node,node,0,0]
            node.maxpath = 0
            return

        n1 = -1
        n2 = -1
        d1 = -1
        d2 = -1
        anchor1 = None
        anchor2 = None
        node.diameter = -1

        for ch in node.child_node_iter():
               ch_record = ch.record
               n = ch_record[2] + 1
               d = ch.maxdepth + ch.edge_length()
               if n > n1:
                   n2 = n1
                   n1 = n
                   anchor2 = anchor1
                   anchor1 = ch_record[0]
               elif n > n2:
                   n2 = n
                   anchor2 = ch_record[0]
               if d > d1:
                   d2 = d1
                   d1 = d
               elif d > d2:
                   d2 = d
               node.diameter = max(ch.diameter,node.diameter)

        node.diameter = max(d1+d2, node.diameter)
        node.maxdepth = d1

        node.record = [anchor1,anchor2,n1,n1+n2] 


