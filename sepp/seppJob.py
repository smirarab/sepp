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


'''
Created on Jul 13, 2012

@author: Nam Nguyen
'''


import ConfigParser
from optparse import OptionParser, OptionGroup
import dendropy
from sepp.taxonomic import tempfile, read_fasta, write_newick, decompose_tree,\
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
    write_alignment("%s.aligned.fasta" % output, alignment, "fasta")               
    tempfiles.append(prefix + ".original.alignment.stockholm")    
    eTimer = os.times()
    logger.write("Time: Global alignment: %s \n" %  (eTimer[4]-sTimer[4]))
    logger.flush()
    [remove_temp(d) for d in tempfiles]

def local_align_local_place_combined_align(tree_file, alignment_file, fragment_file, output, logger, merge_tem_dir,
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

    handle = open(os.path.join(merge_tem_dir,output + ".labeled.tree"), "w")
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
        handle = open(os.path.join(merge_tem_dir,output + "." + str(key) + ".labeled.tree"), "w")
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
                    "%s.%s.aligned.fasta" % (output, str(tree_idx)),
            raxml_file, os.path.join(merge_temp,"%s.%s.json" % (output, str(tree_idx))), pckg=pckg)


    eTimer = os.times()    
    
    logger.write("Time: Total Pplacer: %s \n" %  (eTimer[4]-sTimer[4]))
        
    if (merge):
        sTimer = os.times()        
        merge_trees(merge_temp, os.path.join(merge_temp,output + ".labeled.tree"), "%s.json" % output)
        eTimer = os.times()
        logger.write("Time: Total merge: %s \n" %  (eTimer[4]-sTimer[4]))
    logger.write("Time: Total time: %s \n" %  (eTimer[4]-initTime[4]))

def local_align_local_place_align(tree_file, alignment_file, fragment_file, output, logger,
                  merge_tem_dir = None,tempdir=None, 
                 size=100, strategy="centroid", filters=True, elim=0.01,
                 global_align=True):
                   
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

    handle = open(os.path.join(merge_tem_dir,output + ".labeled.tree"), "w")
    #print tree._tree.as_ascii_plot(show_internal_node_labels=True)
    write_newick(tree._tree, handle)
    handle.write(";\n")    
    handle.close()        
    
    eTimer = os.times()
    logger.write("Time: Decompose tree: %s \n" %  (eTimer[4]-sTimer[4]))

    tree_map = (decompose_tree(tree, size, {}, strategy=strategy))
    
    for key in tree_map.keys(): 
        handle = open(os.path.join(merge_tem_dir,output + "." + str(key) + ".labeled.tree"), "w")
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
    initTime = os.times()
    tree = PhylogeneticTree(dendropy.Tree(stream=open(tree_file), schema="newick"))
    en = 0
    temp_files = []
    for e in tree._tree.postorder_edge_iter():
        e.label = en
        en += 1
    
    labeled_tree_path = os.path.join(merge_temp, output + ".labeled.tree")
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
            raxml_file, os.path.join(merge_temp,"%s.%s.json" % (output, str(tree_idx))), pckg=pckg)

    
    if (global_align and os.path.isfile("%s.unmatched.combined.sto" % output)):
        print "Run Pplacer globally on all remaining fragmnets"
        run_pplacer(tree_file, "%s.unmatched.aligned.fasta" % output,
                    "%s.unmatched.aligned.fasta" % output,
                     raxml_file, os.path.join(merge_temp,"%s.json.unmatched" % (output)), pckg=pckg)
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

def parseOptions (commandLine = None):
    '''Parse command line for options'''

    desc = ' '.join(['This script runs the SEPP algorithm on an input tree, alignment, fragment file, and RAxML info file.'])

    parser = OptionParser(usage = "usage: %prog [options] -t tree_file -a alignment_file -f fragment_file -r raxml_file -o output", 
                          version = "%prog 1.0", 
                          description = desc)

    parser.set_defaults(size = None,
                        superSize = None,
                        output = "output",
                        keep_align = False,
                        tempdir = tempfile.tempdir)    

    group4InfoString = ' '.join(["These options determine the alignment decomposition size and", 
                                 "taxon insertion size.  If None is given, then the default",
                                 "is to align/place at 10% of total taxa.  The alignment decomosition size must be",
                                 "less than the taxon insertion size."])
    group4 = OptionGroup(parser, "Decomposition Options".upper(), group4InfoString)                                 
    
    group4.add_option("-A", "--alignmentSize", type = "int", 
                      dest = "size", metavar = "N", 
                      help = "max alignment subset size of N"
                             "[default: %default]")    
    group4.add_option("-P", "--placementSize", type = "int", 
                      dest = "superSize", metavar = "N", 
                      help = "max placement subset size of N"
                             "[default: %default]")    
    parser.add_option_group(group4)                             
    
    group5InfoString = ' '.join(["These options control output."])
    group5 = OptionGroup(parser, "Output Options".upper(), group5InfoString)                                 
    
    group5.add_option("-d", "--tempdir", 
                      dest = "tempdir", metavar = "DIR", 
                      help = "Tempfile files written to DIR"
                             "[default: %default]")    
    group5.add_option("-o", "--output", 
                      dest = "output", metavar = "OUTPUT", 
                      help = "output with prefix OUTPUT"
                             "[default: %default]")    
    group5.add_option("-k", "--keep_align", 
                      dest = "keep_align", 
                      help = "Flag to keep alignment files"                      
                             "[default: %default]")
    parser.add_option_group(group5)                             
                             
    group6InfoString = ' '.join(["These options control input."])
    group6 = OptionGroup(parser, "Input Options".upper(), group6InfoString)                                 
    
    group6.add_option("-c", "--config", 
                      dest = "config_file", metavar = "CONFIG", 
                      help = "Input config file"
                             "[default: %default]")    
    group6.add_option("-t", "--tree", 
                      dest = "tree_file", metavar = "TREE", 
                      help = "Input tree file"
                             "[default: %default]")    
    group6.add_option("-r", "--raxml", 
                      dest = "raxml_file", metavar = "RAXML", 
                      help = "RAxML_info file"
                             "[default: %default]")    
    group6.add_option("-a", "--alignment", 
                      dest = "alignment_file", metavar = "ALIGN", 
                      help = "Aligned fasta file"
                             "[default: %default]")    
    group6.add_option("-f", "--fragment", 
                      dest = "fragment_file", metavar = "FRAG", 
                      help = "fragment file"
                             "[default: %default]")                                                          
    group6.add_option("-p", "--package", 
                      dest = "package", metavar = "PKG", 
                      help = "package directory"
                             "[default: %default]")                                                          
                             

    parser.add_option_group(group6)
                                 
    if commandLine:
        (options, args) = parser.parse_args(commandLine)
    else:
        (options, args) = parser.parse_args()

#    if len(args) != 1:
#        parser.error("Incorrect number of arguments. Try the -h flag for help.")

#    input = args[0]

    return (options)


def checkOptions(options):
    supply = ""
    if (options.tree_file is None):
        supply = supply + " tree file "
    if (options.alignment_file is None):
        supply = supply + " alignment file "
    if (options.raxml_file is None):
        supply = supply + " raxml file ";
    if (options.fragment_file is None):
        supply = supply + " fragment file "
    if (supply != ""):
        print "Failed to supply: %s\n" % supply
        exit()    

def parseConfig(options):
    config = ConfigParser.ConfigParser()
    config.read(options.config_file)
    #sections = ['input','output','algorithm']
    for (k,v) in options.__dict__.items():    
        try:
            if (k.endswith("file")):
                value = config.get('input', k)
                if ((not value is None and value != "None") and (v is None)):
                    options.__dict__[k] = value
            elif (k.endswith("ize")):    
                if (k == "size"):
                    value = config.get('algorithm', 'alignment_size')
                else:
                    value = config.get('algorithm', 'placement_size')
                if ((not value is None and value != "None") and (v is None)):
                    options.__dict__[k] = value
            else:    
                value = config.get('output', k)
                if ((not value is None and value != "None") and (v is None)):
                    options.__dict__[k] = value
                #Special case for output file
                if ((not value is None and value != "None") and (v == "output") and (k == "output")):
                    options.__dict__[k] = value
        except:      
            continue
    return options

  
def run_with_arguments():
    
    #Increase recursion limit, can hit on very large datasets
    os.sys.setrecursionlimit(1000000)

    (options) = parseOptions()

    #Now parse the config file, only overwriting changes that are 
    #not specified in the options
    if not options.config_file is None:
        options = parseConfig(options)
    
    #Check all key info given
    checkOptions(options)
    
    #Get taxa sizes, inefficient, can be done differently
    alignment = read_fasta(options.alignment_file);
    total = len(alignment)
    max_subset = math.log(total, 2)
    options.method = None
    
    
    
    #If sizes are not set, then use 10% rule
    if (options.size is None and options.superSize is None):
        options.size = int(total*.10)
        options.superSize = int(total*.10)
        options.method = "local"        
    if (options.superSize is None):
        options.superSize = total;
        
    size = int(options.size)
    superSize = int(options.superSize)
    #Determine algorithm based on size setting, will recode this since this is silly
    if (options.size > options.superSize):
        print "Alignment size must be smaller than or equal placment size\n"
        exit()
    elif ((options.superSize >= total)):
        options.method = "nammy"
    elif ((options.size < options.superSize)):
        options.method = "combined"
    elif ((options.size == options.superSize)):
        options.method = "local"
    else:
        #Should not reach here
        print "Error in size setting\n"
        exit()                              

    tree_file = options.tree_file
    alignment_file = options.alignment_file
    raxml_file = options.raxml_file
    output = options.output
    fragment_file = options.fragment_file
    ''' Don't place sequences that HMMER cannot mathc to any subset'''
    global_align = False
    pckg = options.package

    tempdir = options.tempdir
    print tempdir
    #Hard coded filter e cutoff so that all fragments are scored
    elim = "99999999"
    #Filters are off by default
    filters = False
    #always merge
    merge = True
    #always output labeled tree
    label_tree = True
    method = options.method
    threads = 4
    
    strategy = "centroid"
    clean = True
    
    #Getting the output directory
    idx = output.rfind("/")
    outdir = "."
    outfile = output
    if (idx != -1):
        outdir = output[:output.rfind("/")]
    outfile = output[(output.rfind("/")+1):]
              
    logger = sys.stderr
    
    print "Running %s %s %s" % (options.method, options.size, options.superSize)
    #Need to change it from alignment/tree in separate steps back into one, right now a little inefficient
    #clean the output directory.  This can get very messy, will fix later how files are kept track with later
    if (method == "nammy"):
        local_align_global_place_align(tree_file, alignment_file, fragment_file, output, logger, 
            size=size, strategy=strategy, filters=filters, elim=elim,
            global_align=global_align, tempdir=tempdir)
        global_placement(tree_file, raxml_file, output, logger, pckg=pckg, suffix="alignment.sto", ref_suffix="aligned.fasta")
        if (os.path.isfile("%s.meta" % output)):
            remove_temp("%s.meta" % output)
        if (options.keep_align == False):
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
            if (re.search(outfile + '\.\d+' + '(\.labeled\.tree|\.merged\.json|\.json)', f) is not None):
                remove_temp(outdir+"/"+f)
            elif (re.search(outfile + '\.tree\.' + '\d+', f) is not None):
                remove_temp(outdir+"/"+f)
            elif (re.search(outfile + '\.\d+\.combined\.sto', f) is not None and options.keep_align == False):
                remove_temp(outdir+"/"+f)                  
          
    #elif (method == "papara"):
        #if (make_align):          
            #papara_align(tree_file, alignment_file, fragment_file, output, logger, tempdir=tempdir, threads=threads)                  
        #if (make_tree):
            #global_placement(tree_file, raxml_file, output, logger, pckg=pckg, suffix="combined.alignment.sto")          
    elif (method == "pplacer"):
        global_alignment(tree_file, alignment_file, fragment_file, output, logger, tempdir=tempdir)
        global_placement(tree_file, raxml_file, output, logger, pckg=pckg, suffix="combined.alignment.sto", ref_suffix="aligned.fasta")
    
