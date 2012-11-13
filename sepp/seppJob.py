##########################################################################, a#
##    Copyright 2012 Siavash Mirarab, Nam Nguyen, and Tandy Warnow.
##    This file is part of SEPP.
##
##    SEPP is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    SEPP is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with SEPP.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################
import math
import sys
import re
import shutil


'''
Created on Jul 13, 2012

@author: Nam Nguyen
'''


from config import options
import dendropy
from sepp.taxonomic import tempfile, _read_fasta, write_newick, decompose_tree,\
    NammyClass, write_membership, run_pplacer, merge_trees, write_alignment,\
    hmmr_profile, hmmr_align
from Bio import AlignIO  
from satelib.tree import *
from sepp.filemgr import remove_temp
import os

def global_alignment(tree_file, alignment_file, fragment_file, output, logger, tempdir="/scratch/cluster/namphuon/test2/"):

    tempfile.tempdir = tempdir
    f = tempfile.NamedTemporaryFile()    
    prefix = f.name        
    f.close()
    tempfiles = [];
    
    initTime = os.times()
    sTimer = os.times()
    alignment = _read_fasta(alignment_file);
    eTimer = os.times()
    logger.write("Time: Reading alignment: %s \n" %  (eTimer[4]-sTimer[4]))
    #print alignment.keys()
    
    sTimer = os.times()
    fragments = _read_fasta(fragment_file);
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
    write_alignment("%s.aligned.fasta" % output, alignment, "fasta")               
    tempfiles.append(prefix + ".original.alignment.stockholm")    
    eTimer = os.times()
    logger.write("Time: Global alignment: %s \n" %  (eTimer[4]-sTimer[4]))
    logger.flush()
    [remove_temp(d) for d in tempfiles]

def local_align_local_place_combined_align(tree_file, alignment_file, fragment_file, output, logger, merge_tem_dir,
                 size=100, super_size=1000, strategy="centroid", filters=True, elim=0.01,
                 global_align=True, tempdir="/scratch/cluster/namphuon/test2/"):
    prefix = os.path.basename(output)

    initTime = os.times()
    sTimer = os.times()
    alignment = _read_fasta(alignment_file);
    eTimer = os.times()
    logger.write("Time: Reading alignment: %s \n" %  (eTimer[4]-sTimer[4]))
    
        
    sTimer = os.times()
    fragments = _read_fasta(fragment_file);
    eTimer = os.times()
    logger.write("Time: Reading fragments: %s \n" %  (eTimer[4]-sTimer[4]))
                
    sTimer = os.times()
    tree = PhylogeneticTree(dendropy.Tree(stream=open(tree_file), schema="newick", preserve_underscores=True))
    en = 0
    for e in tree._tree.postorder_edge_iter():
        e.label = en
        en += 1

    handle = open(os.path.join(merge_tem_dir,prefix + ".labeled.tree"), "w")
    #print tree._tree.as_ascii_plot(show_internal_node_labels=True)
    write_newick(tree._tree, handle)
    handle.write(";\n")    
    handle.close()        
    
    eTimer = os.times()
    logger.write("Time: Decompose tree: %s \n" %  (eTimer[4]-sTimer[4]))

    #Decompose the tree into the supersets, copy each one since decomposing further destroys the trees in memory
    #TODO:  We don't need to keep the tree structure, we should just create a wrapper class rather than all
    #    this pointless copying
    super_tree_map = (decompose_tree(tree, super_size, {}, strategy=strategy))
    
    for key in super_tree_map.keys(): 
        handle = open(os.path.join(merge_tem_dir,prefix + "." + str(key) + ".labeled.tree"), "w")
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

def local_align_local_place_combined_tree(tree_file, raxml_file, output, logger,merge_temp, 
                size=100, super_size=1000, strategy="centroid", global_align=True, merge=True, clean=False,
                 pckg=""):
    
    prefix = os.path.basename(output)
    initTime = os.times()
    tree = PhylogeneticTree(dendropy.Tree(stream=open(tree_file), schema="newick", preserve_underscores=True))
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
                    "%s.%s.aligned.fasta" % (output, str(tree_idx)),
            raxml_file, os.path.join(merge_temp,"%s.%s.json" % (prefix, str(tree_idx))), pckg=pckg)


    eTimer = os.times()    
    
    logger.write("Time: Total Pplacer: %s \n" %  (eTimer[4]-sTimer[4]))
        
    if (merge):
        sTimer = os.times()        
        merge_trees(merge_temp, os.path.join(merge_temp,prefix + ".labeled.tree"), "%s.json" % output)
        eTimer = os.times()
        logger.write("Time: Total merge: %s \n" %  (eTimer[4]-sTimer[4]))
    logger.write("Time: Total time: %s \n" %  (eTimer[4]-initTime[4]))

def local_align_local_place_align(tree_file, alignment_file, fragment_file, output, logger,
                  merge_tem_dir = None,tempdir=None, 
                 size=100, strategy="centroid", filters=True, elim=0.01,
                 global_align=True):
                   
    prefix = os.path.basename(output)
    
    initTime = os.times()
    sTimer = os.times()
    alignment = _read_fasta(alignment_file);
    eTimer = os.times()
    logger.write("Time: Reading alignment: %s \n" %  (eTimer[4]-sTimer[4]))
    #logger.write("Alignment keys:\n %s\n" %alignment.keys())   

    sTimer = os.times()
    fragments = _read_fasta(fragment_file);
    eTimer = os.times()
    logger.write("Time: Reading fragments: %s \n" %  (eTimer[4]-sTimer[4]))
    #logger.write("Fragment keys:\n %s\n" %fragments.keys())   
                
    sTimer = os.times()
    tree = PhylogeneticTree(dendropy.Tree(stream=open(tree_file), schema="newick", preserve_underscores=True))
    en = 0
    for e in tree._tree.postorder_edge_iter():
        e.label = en
        en += 1
    #logger.write("Tree :\n %s\n" %[x.taxon.label for x in tree._tree.leaf_nodes()])   
    
    handle = open(os.path.join(merge_tem_dir,prefix + ".labeled.tree"), "w")
    #print tree._tree.as_ascii_plot(show_internal_node_labels=True)
    write_newick(tree._tree, handle)
    handle.write(";\n")    
    handle.close()        
    
    eTimer = os.times()
    logger.write("Time: Decompose tree: %s \n" %  (eTimer[4]-sTimer[4]))

    tree_map = (decompose_tree(tree, size, {}, strategy=strategy))
    
    for key in tree_map.keys(): 
        handle = open(os.path.join(merge_tem_dir,prefix + "." + str(key) + ".labeled.tree"), "w")
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
  
def local_align_local_place_tree(tree_file, raxml_file, output, logger, merge_temp, 
                 size=100, strategy="centroid", global_align=True, 
                 merge=True, clean=False, pckg=""):
    prefix = os.path.basename(output)
    initTime = os.times()
    tree = PhylogeneticTree(dendropy.Tree(stream=open(tree_file), schema="newick", preserve_underscores=True))
    en = 0
    temp_files = []
    for e in tree._tree.postorder_edge_iter():
        e.label = en
        en += 1
    
    labeled_tree_path = os.path.join(merge_temp, prefix + ".labeled.tree")
    handle = open(labeled_tree_path,"w")
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
        temp_files.append("%s.tree.%s" % (output, str(tree_idx)))
        temp_files.append("%s.%s.combined.sto" % (output, str(tree_idx)))
        run_pplacer("%s.tree.%s" % (output, str(tree_idx)), "%s.%s.combined.sto" % (output, str(tree_idx)),
                    "%s.%s.aligned.fasta" % (output, str(tree_idx)),
            raxml_file, os.path.join(merge_temp,"%s.%s.json" % (prefix, str(tree_idx))), pckg=pckg)

    
    if (global_align and os.path.isfile("%s.unmatched.combined.sto" % output)):
        print "Run Pplacer globally on all remaining fragmnets"
        run_pplacer(tree_file, "%s.unmatched.aligned.fasta" % output,
                    "%s.unmatched.aligned.fasta" % output,
                     raxml_file, os.path.join(merge_temp,"%s.json.unmatched" % (prefix)), pckg=pckg)
    eTimer = os.times()
    logger.write("Time: Total Pplacer: %s \n" %  (eTimer[4]-sTimer[4]))
    
    if (merge):
        sTimer = os.times() 
        merge_trees(merge_temp, labeled_tree_path, "%s.json" % output)            
        eTimer = os.times()
        logger.write("Time: Total merge: %s \n" %  (eTimer[4]-sTimer[4]))
    
    for x in temp_files:
        remove_temp(x)
    logger.write("Time: Total tree: %s \n" %  (eTimer[4]-initTime[4]))

def local_align_global_place_align(tree_file, alignment_file, fragment_file, output, logger, 
                 size=100, strategy="centroid", filters=True, elim=0.01,
                 global_align=True, tempdir=""):
    initTime = os.times()
    sTimer = os.times()
    alignment = _read_fasta(alignment_file);
    eTimer = os.times()
    logger.write("Time: Reading alignment: %s \n" %  (eTimer[4]-sTimer[4]))
    
        
    sTimer = os.times()
    fragments = _read_fasta(fragment_file);
    eTimer = os.times()
    logger.write("Time: Reading fragments: %s \n" %  (eTimer[4]-sTimer[4]))
                
    sTimer = os.times()
    tree = PhylogeneticTree(dendropy.Tree(stream=open(tree_file), schema="newick", preserve_underscores=True))
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
    write_alignment("%s.aligned.fasta" % output, final_alignment, "fasta")
    nammy.clean_files()
    logger.flush()
    
    
def global_placement(tree_file, raxml_file, output, logger, pckg="", suffix="", ref_suffix=""):               
    sTimer = os.times()
    print "Run Pplacer"
    sTimer = os.times()
    run_pplacer(tree_file, "%s.%s" % (output, suffix), "%s.%s" % (output, ref_suffix), raxml_file, "%s.json" % output, pckg=pckg)    
    eTimer = os.times()
    logger.write("Time: Total Pplacer: %s \n" %  (eTimer[4]-sTimer[4]))
    logger.flush()    

#    if len(args) != 1:
#        parser.error("Incorrect number of arguments. Try the -h flag for help.")

#    input = args[0]

    return (_options_singelton)


def checkOptions(_options_singelton):
    supply = ""
    if (_options_singelton.tree_file is None):
        supply = supply + " tree file,"
    if (_options_singelton.alignment_file is None):
        supply = supply + " alignment file,"
    if (_options_singelton.raxml_file is None):
        supply = supply + " raxml file,";
    if (_options_singelton.fragment_file is None):
        supply = supply + " fragment file"
    if (supply != ""):
        print  >>sys.stderr, "Failed to supply: %s\nRun with -h option to see a list of _options_singelton." % supply
        exit(1)    

  
def run_with_arguments():
    
    #Increase recursion limit, can hit on very large datasets
    os.sys.setrecursionlimit(1000000)

    #Check all key info given
    checkOptions(options)
    
    #Get taxa sizes, inefficient, can be done differently
    alignment = alignment.read_fasta(options['alignment_file']);
    total = len(alignment)
    _options_singelton.method = None
    
    
    
    #If sizes are not set, then use 10% rule
    if (_options_singelton.size is None and _options_singelton.superSize is None):
        _options_singelton.size = int(total*.10)
        _options_singelton.superSize = int(total*.10)
        _options_singelton.method = "local"        
    if (_options_singelton.superSize is None):
        _options_singelton.superSize = total;
        
    size = int(_options_singelton.size)
    superSize = int(_options_singelton.superSize)
    #Determine algorithm based on size setting, will recode this since this is silly
    if (_options_singelton.size > _options_singelton.superSize):
        print  >>sys.stderr, "Alignment size must be smaller than or equal placment size\n"
        exit(1)
    elif ((_options_singelton.superSize >= total)):
        _options_singelton.method = "nammy"
    elif ((_options_singelton.size < _options_singelton.superSize)):
        _options_singelton.method = "combined"
    elif ((_options_singelton.size == _options_singelton.superSize)):
        _options_singelton.method = "local"
    else:
        #Should not reach here
        print  >>sys.stderr, "Error in size setting\n"
        exit(1)                              

    tree_file = _options_singelton.tree_file
    alignment_file = _options_singelton.alignment_file
    raxml_file = _options_singelton.raxml_file
    output = _options_singelton.output
    fragment_file = _options_singelton.fragment_file
    ''' Don't place sequences that HMMER cannot mathc to any subset'''
    global_align = False
    pckg = None #_options_singelton.package

    tempdir = _options_singelton.tempdir    
    #Hard coded filter e cutoff so that all fragments are scored
    elim = "99999999"
    #Filters are off by default
    filters = False
    #always merge
    merge = True
    #always output labeled tree
    #label_tree = True
    method = _options_singelton.method
    #threads = 4
    
    strategy = "centroid"
    clean = True
    
    #Getting the output directory
    outdir = os.curdir if _options_singelton.outdir is None else _options_singelton.outdir
    if os.path.split(output)[0] != '':
        print >>sys.stderr, "Output prefix cannot be a directory: %s" %output
        exit(1)
    outfile = output
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        files = os.listdir(outdir)
        #print files 
        for f in files:
            if (re.search('^' + outfile + '.*', f) is not None):
                print >>sys.stderr, "Output directory [%s] already contains files with prefix [%s]...\nTerminating to avoid loss of existing files." % (outdir,outfile)
                exit(1)
    
    output = os.path.join(outdir,output)
    
    logger = sys.stderr
    
    print "Running %s %s %s" % (_options_singelton.method, _options_singelton.size, _options_singelton.superSize)
    #Need to change it from alignment/tree in separate steps back into one, right now a little inefficient
    #clean the output directory.  This can get very messy, will fix later how files are kept track with later
    if (method == "nammy"):
        local_align_global_place_align(tree_file, alignment_file, fragment_file, output, logger, 
            size=size, strategy=strategy, filters=filters, elim=elim,
            global_align=global_align, tempdir=tempdir)
        global_placement(tree_file, raxml_file, output, logger, pckg=pckg, suffix="alignment.sto", ref_suffix="aligned.fasta")
        if (os.path.isfile("%s.meta" % output)):
            remove_temp("%s.meta" % output)
        if (True):# _options_singelton.keep_align == False):
            if (os.path.isfile("%s.alignment.sto" % output)):
                remove_temp("%s.alignment.sto" % output)
          
    elif (method == "local"):
        merge_temp = tempfile.mkdtemp()
        local_align_local_place_align(tree_file, alignment_file, fragment_file, output, logger,merge_temp, 
            size=size, strategy=strategy, filters=filters, elim=elim,
            global_align=global_align, tempdir=tempdir)
        local_align_local_place_tree(tree_file, raxml_file, output, logger, merge_temp,
            size=size, strategy=strategy, global_align=global_align, merge=merge, 
            clean=clean, pckg=pckg)
        if (os.path.isfile("%s.meta" % output)):
            remove_temp("%s.meta" % output)

    elif (method == "combined"):
        merge_temp = tempfile.mkdtemp()
        local_align_local_place_combined_align(tree_file, alignment_file, fragment_file, output, logger, merge_temp,
            size=size, super_size=superSize, strategy=strategy, filters=filters, elim=elim,
            global_align=global_align, tempdir=tempdir)
        local_align_local_place_combined_tree(tree_file, raxml_file, output, logger, merge_temp, size=size, super_size=superSize, 
            strategy=strategy, global_align=global_align, merge=merge, clean=clean, pckg=pckg)
        #Now clean files
        if (os.path.isfile("%s.meta" % output)):
            remove_temp("%s.meta" % output)
        if (os.path.isfile("%s.labeled.tree" % output)):
            remove_temp("%s.labeled.tree" % output)
        
        files = os.listdir(outdir)        
        for f in files:
            if (re.search(outfile + '\.\d+(\.labeled\.tree|\.merged\.json|\.json)', f) is not None):
                remove_temp(os.path.join(outdir,f))
            elif (re.search(outfile + '\.tree\.' + '\d+', f) is not None):
                remove_temp(os.path.join(outdir,f))
            elif (re.search(outfile + '\.\d+\.combined\.sto', f) is not None and True): # _options_singelton.keep_align == False):
                remove_temp(os.path.join(outdir,f))                  
          
    #elif (method == "papara"):
        #if (make_align):          
            #papara_align(tree_file, alignment_file, fragment_file, output, logger, tempdir=tempdir, threads=threads)                  
        #if (make_tree):
            #global_placement(tree_file, raxml_file, output, logger, pckg=pckg, suffix="combined.alignment.sto")          
    elif (method == "pplacer"):
        global_alignment(tree_file, alignment_file, fragment_file, output, logger, tempdir=tempdir)
        global_placement(tree_file, raxml_file, output, logger, pckg=pckg, suffix="combined.alignment.sto", ref_suffix="aligned.fasta")
    
