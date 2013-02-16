import urllib2
import Bio.Phylo as bp
from Bio.Phylo import Newick
import os
import tarfile


col_delimiter = '\t|\t'
row_delimiter = '\t|\n'
url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
tree_filename = 'ncbi_taxonomy.newick'
tree_format = 'newick'

# download the taxonomy archive
filename = url.split('/')[-1]
if os.path.exists(filename):
    print 'Using existing copy of %s' % filename
else:
    print 'Downloading %s...' % filename
    r = urllib2.urlopen(urllib2.Request(url))
    assert r.geturl() == url
    with open(filename, 'wb') as output_file:
        output_file.write(r.read())
    r.close()

# extract the text dump
for filename in ('nodes.dmp', 'names.dmp'):
    if os.path.exists(filename):
        print 'Using existing copy of %s' % filename
    else:
        print 'Extracting %s...' % filename
        archive = tarfile.open(name=filename, mode='r:gz')
        archive.extract(filename)

# get names for all tax_ids from names.dmp
print 'Getting names...'
scientific_names = {}
common_names = {}
with open('names.dmp') as names_file:
    for line in names_file:
        line = line.rstrip(row_delimiter)
        values = line.split(col_delimiter)
        tax_id, name_txt, _, name_type = values[:4]
        if name_type == 'scientific name':
            scientific_names[tax_id] = name_txt
        elif name_type == 'common name':
            common_names[tax_id] = name_txt

# read all node info from nodes.dmp
print 'Reading taxonomy...'
nodes = {}
with open('nodes.dmp') as nodes_file:
    for line in nodes_file:
        line = line.rstrip(row_delimiter)
        values = line.split(col_delimiter)
        tax_id, parent_id = values[:2]
        this_node = Newick.Clade(name=scientific_names[tax_id])
        nodes[tax_id] = this_node
        this_node.parent = parent_id

# create tree from nodes dictionary
print 'Building tree...'
for node_id, this_node in nodes.items():
    if node_id == this_node.parent:
        root_node = this_node
        print 'Found root.'
    else:
        parent_node = nodes[this_node.parent]
        parent_node.clades.append(this_node)
    del this_node.parent

tree = Newick.Tree(root=root_node)

# write tree to file
print 'Writing %s tree to %s...' % (tree_format, tree_filename)
bp.write([tree], tree_filename, tree_format)

print 'Done!'