import dendropy
import sys

# Arguments: 1: a "reference" tree with internal taxon labels
#               example: 99_otus.tree
#            2: a "new" tree that is a refinement of the refrence tree, but
#               lacks internal nodes
#               example: reference-gg-raxml-bl-rooted.tre
#            3: the output: the new tree with internal nodes added
#               example: reference-gg-raxml-bl-rooted-relabelled.tre
# Assumes there is a single node called k__Archaea to root the output tree on

taxa = dendropy.TaxonNamespace()

t = dendropy.Tree.get_from_path(
    src=sys.argv[1], schema='newick', taxon_namespace=taxa)
t2 = dendropy.Tree.get_from_path(
    src=sys.argv[2], schema='newick', taxon_namespace=taxa)

t.encode_bipartitions()

missing = []
mapped = 0
for n in t.postorder_node_iter():
    if n.is_internal() and n.label is not None:
        # print (n.label)
        n.edge.bipartition.is_mutable = False
        if n.edge.bipartition not in t2.bipartition_edge_map:
            missing.append(n.edge.bipartition)
        else:
            t2n = t2.bipartition_edge_map[n.edge.bipartition].head_node
            t2n.label = n.label
            mapped = mapped + 1

a = t2.find_node_with_label(label="k__Archaea")

t2.reroot_at_node(a, update_bipartitions=False)

t2.write_to_path(dest=sys.argv[3], schema='newick', suppress_rooting=True,
                 suppress_internal_node_labels=False)
