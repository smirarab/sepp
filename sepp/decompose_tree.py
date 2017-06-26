# uym2 added
# June 2017
# utils for tree decomposition


from dendropy import Tree, Node
try:
    from queue import Queue # python 3
except:
    from Queue import Queue # python 2
#from tree import PhylogeneticTree


def decompose_by_diameter(a_tree,max_diam=None,nsubtree=50):
    def __ini_record__():
        for node in a_tree.postorder_node_iter():
               __updateNode__(node)
    def __find_centroid_edge__(t):
        h = 1
        u = t.seed_node.bestLCA.anchor
        while (h < t.seed_node.topo_diam//2):
            h += 1
            u = u.parent_node
        return u.edge

    def __bisect__(t):
        e = __find_centroid_edge__(t)
        
        u = e.tail_node
        v = e.head_node

        u.remove_child(v)
        t1 = Tree(seed_node = v)

        if u.num_child_nodes() == 1:
            p = u.parent_node
            v = u.child_nodes()[0]
            l_v = v.edge_length
            u.remove_child(v)
            if p is None: # u is the seed_node; this means the tree runs out of all but one side
                t.seed_node = v
                return t,t1
            l_u = u.edge_length
            p.remove_child(u)
            p.add_child(v)
            v.edge_length = l_u+l_v
            u = p

        while u is not None:
            __updateNode__(u)
            u = u.parent_node

        return t,t1

    def __clean_up__(t):
        for node in t.postorder_node_iter():
            delattr(node,"anchor")
            delattr(node,"maxheight")
            delattr(node,"maxdepth")
            delattr(node,"diameter")
            delattr(node,"topo_diam")
            delattr(node,"bestLCA")

    def __updateNode__(node):
        if node.is_leaf():
            node.anchor = node
            node.maxheight = 0
            node.maxdepth = 0
            node.diameter = 0
            node.topo_diam = 0
            node.bestLCA = node
            return

        n1 = -1
        n2 = -1
        d1 = -1
        d2 = -1
        anchor1 = None
        anchor2 = None
        node.diameter = 0
        node.topo_diam = 0
        node.bestLCA = None

        for ch in node.child_node_iter():
               n = ch.maxheight + 1
               d = ch.maxdepth + ch.edge_length
               if n > n1:
                   n2 = n1
                   n1 = n
                   anchor2 = anchor1
                   anchor1 = ch.anchor
               elif n > n2:
                   n2 = n
                   anchor2 = ch.anchor
               if d > d1:
                   d2 = d1
                   d1 = d
               elif d > d2:
                   d2 = d
               if ch.topo_diam > node.topo_diam:
                   node.top_diam = ch.topo_diam
                   node.bestLCA = ch.bestLCA
               node.diameter = max(ch.diameter,node.diameter)

        node.diameter = max(d1+d2, node.diameter)
        node.maxdepth = d1
        node.maxheight = n1
        node.anchor = anchor1
        if n1+n2 > node.topo_diam:
            node.topo_diam = n1+n2
            node.bestLCA = node

    tqueue = Queue()
    __ini_record__()
    max_diam = max_diam if max_diam is not None else a_tree.seed_node.diameter/nsubtree
    if a_tree.seed_node.diameter <= max_diam:
        return [a_tree]
        
    treeMap = [] 
    tqueue.put(a_tree)
    while not tqueue.empty():
        t = tqueue.get()
        print(t.seed_node.diameter)
        t1,t2 = __bisect__(t)
        print(t1.seed_node.diameter)
        print(t2.seed_node.diameter)
        print("\n")
        if t1.seed_node.diameter <= max_diam:
             __clean_up__(t1)            
             #treeMap[i] = PhylogeneticTree(t1)
             treeMap.append(t1)
        else:
            tqueue.put(t1)
        if t2.seed_node.diameter <= max_diam:
             __clean_up__(t2)
             #treeMap[i] = PhylogeneticTree(t2)
             treeMap.append(t2)
        else:
            tqueue.put(t2)

    return treeMap
