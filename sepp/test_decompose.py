from sepp.tree import PhylogeneticTree
import dendropy

def bipart_by_length(treelist,maxLength):    
    good_trees = []
    while (len(treelist) != 0):
        tree = treelist.pop()
        edge=get_long_edge(tree,maxLength,3)
        if (edge is not None):
            (t1,t2) = tree.bipartition_by_edge(edge)
            treelist.append(t1)
            treelist.append(t2)
        else:
            good_trees.append(tree)
    return good_trees

def get_long_edge(tree, maxLength, minSize):
    for edge in tree._tree.postorder_edge_iter():
        if edge.tail_node is None:
            continue
        onesideSize = len(edge.head_node.leaf_nodes())
        if edge.length is not None and edge.length > maxLength:
            return edge
            

if __name__ == '__main__':   
    tree = PhylogeneticTree(dendropy.Tree(stream=open('/home/n-z/namphuon/no_backup/data/silva/silva_nr99_reference.tre'), schema="newick", preserve_underscores=True))    
    trees = bipart_by_length([tree],0.15577181378590388)
    for x in xrange(len(trees)):
        trees[x].write_newick_to_path("/home/n-z/namphuon/no_backup/data/silva/trees/%d.tree" % x)            