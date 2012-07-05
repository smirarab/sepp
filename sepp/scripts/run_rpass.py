#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A test script
   Created: June 7. 2011
   Author: Nam Nguyen
"""

import sys
import getopt, tempfile
from satelib.tree import *
import dendropy
from sepp.taxonomic import *
from Bio import AlignIO  
import copy
import ConfigParser

def global_alignment(tree_file, alignment_file, fragment_file, output, logger, tempdir="/scratch/cluster/namphuon/test2/"):

	tempfile.tempdir = tempdir
	f = tempfile.NamedTemporaryFile()	
        prefix = f.name        
	f.close()
	tempfiles = [];
	
	initTime = os.times()
	sTimer = os.times()
	alignment = read_fasta(alignment_file);
	eTimer = os.times()
	logger.write("Time: Reading alignment: %s \n" %  (eTimer[4]-sTimer[4]))
	
	sTimer = os.times()
	fragments = read_fasta(fragment_file);
	eTimer = os.times()
	logger.write("Time: Reading fragments: %s \n" %  (eTimer[4]-sTimer[4]))
	
	sTimer = os.times()	
	write_alignment(prefix + ".original.alignment.stockholm", alignment, 'stockholm')
	hmmr_profile(prefix + ".original.alignment.stockholm", prefix + ".full.profile", "stockholm")
	tempfiles.append(prefix + ".full.profile")
	eTimer = os.times()
	logger.write("Time: Profile entire alignment: %s \n" %  (eTimer[4]-sTimer[4]))
	
	sTimer = os.times()
	hmmr_align(fragment_file, prefix + ".full.profile", 
		   "%s.combined.alignment.sto" % output, original_file = prefix + ".original.alignment.stockholm", trim=True)		   	
	tempfiles.append(prefix + ".original.alignment.stockholm")	
	eTimer = os.times()
	logger.write("Time: Global alignment: %s \n" %  (eTimer[4]-sTimer[4]))
	logger.flush()
	[os.remove(d) for d in tempfiles]

def local_align_local_place_combined_align(tree_file, alignment_file, fragment_file, output, logger, 
			     size=100, super_size=1000, strategy="centroid", filters=True, elim=0.01,
			     global_align=True, tempdir="/scratch/cluster/namphuon/test2/"):
	initTime = os.times()
	sTimer = os.times()
	alignment = read_fasta(alignment_file);
	eTimer = os.times()
	logger.write("Time: Reading alignment: %s \n" %  (eTimer[4]-sTimer[4]))
	
		
	sTimer = os.times()
	fragments = read_fasta(fragment_file);
	eTimer = os.times()
	logger.write("Time: Reading fragments: %s \n" %  (eTimer[4]-sTimer[4]))
				
	sTimer = os.times()
	tree = PhylogeneticTree(dendropy.Tree(stream=open(tree_file), schema="newick"))
	en = 0
	for e in tree._tree.postorder_edge_iter():
		e.label = en
		en += 1

	handle = open(output + ".labeled.tree", "w")
	#print tree._tree.as_ascii_plot(show_internal_node_labels=True)
	write_newick(tree._tree, handle)
	handle.write(";\n")	
	handle.close()		
	
	eTimer = os.times()
	logger.write("Time: Decompose tree: %s \n" %  (eTimer[4]-sTimer[4]))

	#Decompose the tree into the supersets, copy each one since decomposing further destroys the trees in memory
	#TODO:  We don't need to keep the tree structure, we should just create a wrapper class rather than all
	#	this pointless copying
	super_tree_map = (decompose_tree(tree, super_size, {}, strategy=strategy))
	
	for key in super_tree_map.keys(): 
	    handle = open(output + "." + str(key) + ".labeled.tree", "w")
	    write_newick(super_tree_map[key]._tree, handle)
	    handle.write(";\n")
	    handle.close()
	
	super_tree_copy = {}
	for key in super_tree_map:
	      super_tree_copy[key] = copy.deepcopy(super_tree_map[key])
	
	#Now for each tree, decompose into the smaller set, keeping track of which tree the final subtrees
	#came from in original map, keeping all subtrees in tree map.  
	tree_map = {}
	original_map = {}
	counter = 0
	for key in super_tree_copy.keys():
	      temp_map = (decompose_tree(super_tree_copy[key], size, {}, strategy=strategy))
	      for temp_key in temp_map.keys():
		   tree_map[counter] = temp_map[temp_key]
		   original_map[counter] = key
		   counter = counter + 1
		   
	#Find the set membership for the smallest subsets
	print "Finding membership"
	nammy = NammyClass(tree_map, alignment, fragments, temp_directory=tempdir, logger=logger)
		
	sTimer = os.times()
	(local_set) = nammy.find_membership(filters=filters, elim=elim)	
	eTimer = os.times()
	logger.write("Time: Total membership time: %s \n" %  (eTimer[4]-sTimer[4]))
	
	handle = open("%s.meta" % output, "w")
	write_membership(handle, local_set, tree_map)
	handle.close()	
		
	print "Align locally using membership results"
	nammy.subset_align(local_set, super_tree_map, original_map, output_name = ("%s" % (output)), ext = "combined.sto")
	
	eTimer = os.times()
	logger.write("Time: Total alignment everything: %s \n" %  (eTimer[4]-initTime[4]))
	nammy.clean_files()

def local_align_local_place_combined_tree(tree_file, raxml_file, output, logger, size=100, 
				 super_size=1000, strategy="centroid", global_align=True, merge=True, clean=False,
				 pckg=""):
	
	initTime = os.times()
	tree = PhylogeneticTree(dendropy.Tree(stream=open(tree_file), schema="newick"))
	en = 0
	for e in tree._tree.postorder_edge_iter():
		e.label = en
		en += 1
	super_tree_map = (decompose_tree(tree, super_size, {}, strategy=strategy))
	
	eTimer = os.times()
	logger.write("Time: Decompose tree: %s \n" %  (eTimer[4]-initTime[4]))
	
	print "Run Pplacer locally on each subtree"
	sTimer = os.times()
	for tree_idx in super_tree_map:	      
	      tree_handle = open("%s.tree.%s" % (output, str(tree_idx)), "w")
	      tree_handle.write(super_tree_map[tree_idx].compose_newick() + ";")
	      tree_handle.close()
	      	      
	      run_pplacer("%s.tree.%s" % (output, str(tree_idx)), "%s.%s.combined.sto" % (output, str(tree_idx)),
			  raxml_file, "%s.%s.json" % (output, str(tree_idx)), pckg=pckg)


	eTimer = os.times()	
	
	logger.write("Time: Total Pplacer: %s \n" %  (eTimer[4]-sTimer[4]))
	
	if (merge):
	    sTimer = os.times()	    
	    idx = output.rfind("/")
	    directory = ".";
	    if (idx != -1):
		directory = output[:output.rfind("/")]
	    merge_trees(directory, output + ".labeled.tree")
	    eTimer = os.times()
	    logger.write("Time: Total merge: %s \n" %  (eTimer[4]-sTimer[4]))
	logger.write("Time: Total time: %s \n" %  (eTimer[4]-initTime[4]))

def local_align_local_place_align(tree_file, alignment_file, fragment_file, output, logger, 
			     size=100, strategy="centroid", filters=True, elim=0.01,
			     global_align=True, tempdir="/scratch/cluster/namphuon/test2/"):
			       
	initTime = os.times()
	sTimer = os.times()
	alignment = read_fasta(alignment_file);
	eTimer = os.times()
	logger.write("Time: Reading alignment: %s \n" %  (eTimer[4]-sTimer[4]))
	
	sTimer = os.times()
	fragments = read_fasta(fragment_file);
	eTimer = os.times()
	logger.write("Time: Reading fragments: %s \n" %  (eTimer[4]-sTimer[4]))
				
	sTimer = os.times()
	tree = PhylogeneticTree(dendropy.Tree(stream=open(tree_file), schema="newick"))
	en = 0
	for e in tree._tree.postorder_edge_iter():
		e.label = en
		en += 1

	handle = open(output + ".labeled.tree", "w")
	#print tree._tree.as_ascii_plot(show_internal_node_labels=True)
	write_newick(tree._tree, handle)
	handle.write(";\n")	
	handle.close()		
	
	eTimer = os.times()
	logger.write("Time: Decompose tree: %s \n" %  (eTimer[4]-sTimer[4]))

	tree_map = (decompose_tree(tree, size, {}, strategy=strategy))
	
	for key in tree_map.keys(): 
	    handle = open(output + "." + str(key) + ".labeled.tree", "w")
	    write_newick(tree_map[key]._tree, handle)
	    handle.write(";\n")
	    handle.close()	
	
	print "Run local alignment"
	nammy = NammyClass(tree_map, alignment, fragments, temp_directory=tempdir, logger=logger)
	
	sTimer = os.times()
	(local_set) = nammy.find_membership(filters=filters, elim=elim)

	handle = open("%s.meta" % output, "w")
	write_membership(handle, local_set, tree_map)
	handle.close()
	
	nammy.locally_align(local_set, global_align=global_align)
	nammy.move_alignments(output, extension="combined.sto")
	eTimer = os.times()
	logger.write("Time: Total local align: %s \n" %  (eTimer[4]-initTime[4]))
	nammy.clean_files()
  
def local_align_local_place_tree(tree_file, raxml_file, output, logger, size=100, strategy="centroid", global_align=True, 
				 merge=True, clean=False, pckg=""):
	initTime = os.times()
	tree = PhylogeneticTree(dendropy.Tree(stream=open(tree_file), schema="newick"))
	en = 0
	for e in tree._tree.postorder_edge_iter():
		e.label = en
		en += 1

	handle = open(output + ".labeled.tree", "w")
	#print tree._tree.as_ascii_plot(show_internal_node_labels=True)
	write_newick(tree._tree, handle)
	handle.write(";\n")	
	handle.close()		
	
	eTimer = os.times()
	logger.write("Time: Decompose tree: %s \n" %  (eTimer[4]-initTime[4]))

	tree_map = (decompose_tree(tree, size, {}, strategy=strategy))		
	
	print "Run Pplacer locally on each subtree"
	sTimer = os.times()
	for tree_idx in tree_map:
	      tree_handle = open("%s.tree.%s" % (output, str(tree_idx)), "w")
	      tree_handle.write(tree_map[tree_idx].compose_newick() + ";")
	      tree_handle.close()
	      	      
	      run_pplacer("%s.tree.%s" % (output, str(tree_idx)), "%s.%s.combined.sto" % (output, str(tree_idx)),
			  raxml_file, "%s.%s.json" % (output, str(tree_idx)), pckg=pckg)

	
	if (global_align and os.path.isfile("%s.unmatched.combined.sto" % output)):
	    print "Run Pplacer globally on all remaining fragmnets"
	    run_pplacer(tree_file, "%s.unmatched.combined.sto" % output, raxml_file, "%s.json.unmatched" % (output), pckg=pckg)
	eTimer = os.times()
	logger.write("Time: Total Pplacer: %s \n" %  (eTimer[4]-sTimer[4]))
	
	if (merge):
	    sTimer = os.times()	    
	    idx = output.rfind("/")
	    directory = ".";
	    if (idx != -1):
		directory = output[:output.rfind("/")]		
	    merge_trees(directory, output + ".labeled.tree")			
	    eTimer = os.times()
	    logger.write("Time: Total merge: %s \n" %  (eTimer[4]-sTimer[4]))
	logger.write("Time: Total tree: %s \n" %  (eTimer[4]-initTime[4]))
  
def local_align_global_place_align(tree_file, alignment_file, fragment_file, output, logger, 
			     size=100, strategy="centroid", filters=True, elim=0.01,
			     global_align=True, tempdir="/scratch/cluster/namphuon/test2/"):
	initTime = os.times()
	sTimer = os.times()
	alignment = read_fasta(alignment_file);
	eTimer = os.times()
	logger.write("Time: Reading alignment: %s \n" %  (eTimer[4]-sTimer[4]))
	
		
	sTimer = os.times()
	fragments = read_fasta(fragment_file);
	eTimer = os.times()
	logger.write("Time: Reading fragments: %s \n" %  (eTimer[4]-sTimer[4]))
				
	sTimer = os.times()
	tree = PhylogeneticTree(dendropy.Tree(stream=open(tree_file), schema="newick"))
	eTimer = os.times()
	logger.write("Time: Decompose tree: %s \n" %  (eTimer[4]-sTimer[4]))

	tree_map = (decompose_tree(tree, size, {}, strategy=strategy))
	
	print "Run local alignment, but will global place"
	nammy = NammyClass(tree_map, alignment, fragments, temp_directory=tempdir, logger=logger)
	
	sTimer = os.times()
	(local_set) = nammy.find_membership(filters=filters, elim=elim)
	
	handle = open("%s.meta" % output, "w")
	write_membership(handle, local_set, tree_map)	
	handle.close()

	nammy.locally_align(local_set, global_align=global_align)
	final_alignment = nammy.get_merged_alignments(local_set, global_align=global_align)
	nammy.clean_files()
	eTimer = os.times()
	logger.write("Time: Total local align: %s \n" %  (eTimer[4]-initTime[4]))

	write_alignment("%s.alignment.sto" % output, final_alignment, "stockholm")
	nammy.clean_files()
	logger.flush()
	
def papara_align(tree_file, alignment_file, fragment_file, output, logger, tempdir="/scratch/cluster/namphuon/test2/", threads=4):
	initTime = os.times()	
	sTimer = os.times()	
	run_papara(tree_file, alignment_file, fragment_file, "%s.combined.alignment.sto" % output, tempdir=tempdir, file_format="stockholm", threads=threads)
	eTimer = os.times()
	logger.write("Time: Papara alignment: %s \n" %  (eTimer[4]-initTime[4]))	
	logger.flush()  
	
def global_placement(tree_file, raxml_file, output, logger, pckg="", suffix=""):			   
	sTimer = os.times()
	print "Run Pplacer"
	sTimer = os.times()
	run_pplacer(tree_file, "%s.%s" % (output, suffix), raxml_file, "%s.json" % output, pckg=pckg)	
	eTimer = os.times()
	logger.write("Time: Total Pplacer: %s \n" %  (eTimer[4]-sTimer[4]))
	logger.flush()	

def main(argv):
	tree_file = None
	alignment_file = None
	raxml_file = ""
	log = None
	output = None
	fragment_file = None
	size = 25
	super_size = None
	global_align = True
	pckg = ""
	
	#Basically run pplacer on results
	make_tree = True
	
	#Basically only make alignment
	make_align = True
	tempdir = "/scratch/cluster/namphuon/test2/"
	elim = "0.01"
	filters = True
	merge = True
	label_tree = None
	method = "nammy"
	threads = 4
	
	strategy = "centroid"
	clean = False
	
    	try:
        	opts, args = getopt.getopt(argv, "t:a:r:L:l:o:f:T:s:g:d:x:e:S:m:M:A:B:X:c:R:", 
        	["tree=", "alignment=", "raxml=", "label=",
        	"log=", "output=", "fragment=", "size=", "no=", "global=", "directory=", 
        	"filter=", "elim=", "super=", "merge=", "method=", "make_align=" "break=", "threads=", "clean=", "pckg="])
    	except getopt.GetoptError:
		print "Error!"
        	sys.exit(2)	
	for opt, arg in opts:
		if opt in ("-t", "--tree"):
			tree_file = arg
		elif opt in ("-a", "--alignment"):
			alignment_file = arg
		elif opt in ("-r", "--raxml"):
			raxml_file = arg			
		elif opt in ("-l", "--log"):
			log = arg
		elif opt in ("-o", "--output"):
			output = arg			
		elif opt in ("-f", "--fragment"):			
			fragment_file = arg  
		elif opt in ("-s", "--size"):			
			size = int(arg)
		elif opt in ("-T", "--make_tree"):
			if (arg.upper() == "FALSE"):		  
			    make_tree = False
		elif opt in ("-g", "--global"):
			if (arg.upper() == "FALSE"):
			    global_align = False
		elif opt in ("-c", "--clean"):
			if (arg.upper() != "FALSE"):
			    clean = True
		elif opt in ("-d", "--directory"):
			tempdir = arg
		elif opt in ("-x", "--filter"):			
			if (arg.upper() == "FALSE"):		
			    filters = False
		elif opt in ("-e", "--elim"):
			elim = arg		
		elif opt in ("-S", "--super"):
			super_size = int(arg)
		elif opt in ("-m", "--merge"):
			if (arg.upper() == "FALSE"):		  
			    merge = False
		elif opt in ("-L", "--label"):
			label_tree = arg
		elif opt in ("-M", "--method"):
			method = arg
		elif opt in ("-B", "--break"):
			strategy = arg			
		elif opt in ("-A", "--make_align"):
			if (arg.upper() == "FALSE"):
			    make_align = False
		elif opt in ("-X", "--threads"):
			threads = int(arg)
		elif opt in ("-R", "--pckg"):
			pckg = arg
			
			    
	if (clean):
	    idx = output.rfind("/")
	    outdir = ".";
	    if (idx != -1):
		outdir = output[:output.rfind("/")]
	    os.system("rm %s/*merged*json*" % (outdir));  	      			
			
	if log is None:
	      logger = sys.stderr
	else:
	      logger = open(log, 'w')
	      
	if (method == "nammy"):
	      if (make_align):
		  local_align_global_place_align(tree_file, alignment_file, fragment_file, output, logger, 
			     size=size, strategy=strategy, filters=filters, elim=elim,
			     global_align=global_align, tempdir=tempdir)
	      if (make_tree):		
		  global_placement(tree_file, raxml_file, output, logger, pckg=pckg, suffix="alignment.sto")
	elif (method == "local"):
	      if (make_align):
		  local_align_local_place_align(tree_file, alignment_file, fragment_file, output, logger, 
			     size=size, strategy=strategy, filters=filters, elim=elim,
			     global_align=global_align, tempdir=tempdir)
	      if (make_tree):
		  local_align_local_place_tree(tree_file, raxml_file, output, logger, size=size, strategy=strategy, global_align=global_align, merge=merge, clean=clean, pckg=pckg)
	elif (method == "combined"):
	      if (make_align):
		  local_align_local_place_combined_align(tree_file, alignment_file, fragment_file, output, logger, 
			     size=size, super_size=super_size, strategy=strategy, filters=filters, elim=elim,
			     global_align=global_align, tempdir=tempdir)
	      if (make_tree):
		  local_align_local_place_combined_tree(tree_file, raxml_file, output, logger, size=size, super_size=super_size, 
					       strategy=strategy, global_align=global_align, merge=merge, clean=clean, pckg=pckg)
	elif (method == "papara"):
	      if (make_align):		  
		  papara_align(tree_file, alignment_file, fragment_file, output, logger, tempdir=tempdir, threads=threads)			      
	      if (make_tree):
		  global_placement(tree_file, raxml_file, output, logger, pckg=pckg, suffix="combined.alignment.sto")		  
	elif (method == "pplacer"):
	      if (make_align):
		  global_alignment(tree_file, alignment_file, fragment_file, output, logger, tempdir=tempdir)
	      if (make_tree):
		  global_placement(tree_file, raxml_file, output, logger, pckg=pckg, suffix="combined.alignment.sto")

	if (not log is None):
	    logger.close()
	    
if __name__ == "__main__":
	main(sys.argv[1:])
	
	
