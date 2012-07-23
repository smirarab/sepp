#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script bisects a tree into smaller parts.
   Created: June 7. 2011
   Author: Nam Nguyen
"""

import sys, re, os, tempfile, array
from dendropy import Tree
from satelib.satealignerjob import bisect_tree
from satelib.tree import PhylogeneticTree
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sepp.alignment import Alignment
from sepp.utilities import read_sto_alignment
from sepp import get_setup_path
import shutil
from sepp.filemgr import remove_temp

LIB_PATH = os.path.join(get_setup_path(), "lib")
MERGE_JAR_FILE = os.path.join(LIB_PATH, "merge.jar")

def merge_trees(directory, tree, output):    
    print ("java -jar %s %s %s %s" % (MERGE_JAR_FILE, directory, tree, output))
    os.system("java -jar %s %s %s %s" % (MERGE_JAR_FILE, directory, tree, output))    

def write_newick(tree, out):
    """
    This returns the Node as a NEWICK with edge labels
    """
    if type(tree) is PhylogeneticTree:
        t = tree._tree
    elif type(tree) is Tree:
        t = tree
    else: 
        raise "Object type %s not recognized." % str(type(tree))
        
    node = t.seed_node
    write_newick_node(node, out)

def write_newick_node(node, out):
    child_nodes = node.child_nodes()
    if child_nodes:
        out.write('(')
        f_child = child_nodes[0]
        for child in child_nodes:
            if child is not f_child:
                out.write(',')
            write_newick_node(child, out)
        out.write(')')

    out.write(node.get_node_str())
    e = node.edge
    if e:
        sel = e.length
        if sel is not None:
            s = ""
            try:
                s = float(sel)
                s = str(s)
            except ValueError:
                s = str(sel)
            if s:
                out.write(":%s[%s]" % (s , e.label)) 


def read_membership(membership_file):
    file = open(membership_file, 'r');
    
    set_pattern = re.compile(r"([^:]*):(.*)")
    local_sets = {}
    ordered_sets = {}
    
    for line in file:
        if (line.strip() == "Splits"):
            my_set = ordered_sets
        elif (line.strip() == "Fragment assignments"):
            my_set = local_sets
        else:
            matches = set_pattern.search(line)
            my_set[int(matches.group(1).strip())] = matches.group(2).split()
    return (local_sets, ordered_sets)
      
  
def write_membership(output_stream, local_sets, tree_map):
    output_stream.write("Splits\n")
    for tree_idx in tree_map.keys():
        subset = [i.taxon.label for i in tree_map[tree_idx]._tree.leaf_nodes()]
        output_stream.write("%s: %s\n" % (str(tree_idx), " ".join(subset)))
    output_stream.write("Fragment assignments\n")
    for item in local_sets.keys():
        output_stream.write("%s: %s\n" % (str(item), " ".join(local_sets[item])))

"""
    This function returns the file name, minus the extension.
"""
def get_filename(file_name):
    end_idx = len(file_name)
    start_idx = 0
    if (file_name.rfind(".") != -1):
        end_idx = file_name.rfind(".")
    if (file_name.rfind("/") != -1):    
        start_idx = file_name.rfind("/") + 1
    return file_name[start_idx:end_idx]


"""
    This function decomposes the tree until all subtrees are smaller than the max size.  Two
    possible decompositions strategies can used: "centroid" and "longest".  Returns a
    map containing the subtrees, in an ordered fashion.
"""
def decompose_tree(tree, maxSize=25, tree_map={}, strategy="centroid"):          
    tree._tree.deroot()
    if (tree.count_leaves() > maxSize):    
        (t1, t2) = bisect_tree(tree, strategy)
        decompose_tree(t1, maxSize, tree_map, strategy)
        decompose_tree(t2, maxSize, tree_map, strategy)
    else:
        tree_map[len(tree_map)] = tree
    return tree_map
    
"""
    This function runs pplacer on a newick format tree, a fasta alignment file containing
    full and fragmentary sequences, a raxml info file.  the output is a json file.
    TODO Fix input arguments
"""
def run_pplacer(tree="", alignment="",ref_alignment="", raxml="", output="", tempdir=None, pckg="", version="pplacer", options=""):

    executalbe = version #os.path.join(LIB_PATH,version)
    ##Get full path to files
    if (tempdir is None):
        temp_directory = tempfile.tempdir         
        
    full_tree_path = os.path.abspath(tree)
    full_raxml_path = os.path.abspath(raxml)
    full_alignment_path = os.path.abspath(alignment)
    full_ref_alignment_path = os.path.abspath(ref_alignment)                          

    full_output_path = os.path.abspath(output)        

    #Change to temporary directory to make sure no name collision crash occurs
    current_directory = os.getcwd()
    
    if (tree != ""):
        tree = "-t %s" % full_tree_path
    elif (tree is None):
        tree = ""         
    if (raxml != ""):
        raxml = "-s %s" % full_raxml_path
    elif (raxml is None):
        raxml = ""         
    if (ref_alignment != ""):
        ref_alignment = "-r %s" % full_ref_alignment_path
    if (pckg is None):
        pckg = ""
    elif (pckg != ""):
        full_pckg_path = os.path.abspath(pckg)
        pckg = "-c %s" % full_pckg_path
        
        
    new_temp = tempfile.mkdtemp()
    os.chdir(new_temp)
        
    print("%s --out-dir %s %s %s %s %s %s %s" % (executalbe, new_temp , pckg, tree, raxml, ref_alignment, options, full_alignment_path));
    os.system("%s --out-dir %s %s %s %s %s %s %s" % (executalbe, new_temp, pckg, tree, raxml, ref_alignment, options, full_alignment_path));
    if (os.path.exists("%s.json" % get_filename(alignment))):
        print "%s.json to be moved to %s " % (get_filename(alignment), full_output_path)
        shutil.move("%s.json" % get_filename(alignment), full_output_path)
    if (os.path.exists("%s.jplace" % get_filename(alignment))):
        print "%s.jplace to be moved to %s " % (get_filename(alignment), full_output_path)
        shutil.move("%s.jplace" % get_filename(alignment), full_output_path)
    os.chdir(current_directory)
    remove_temp(new_temp)


def read_alignment(fn, format="fasta"):
    try:
        ioA = AlignIO.read(fn, format)
        alignment = Alignment()
        for row in ioA:
            alignment[row.id] = row.seq.data.upper()
    except Exception as inst:
        print type(inst)   
        print inst.args    
        print inst
        raise RuntimeError("Could not read %s alignment" % format)
        
    return alignment  

"""     
    This function reads file in fasta format and returns
    a dictionary containing the alignment
"""
def read_fasta(input_file):
    file = open(input_file, 'r')
    alignment = {}
    counter = 0;
    name_pattern = re.compile(r">([^\n\r\f]*)")    
    line = file.readline()
    while (not (line == "") and line.find(">") != -1):
        m = name_pattern.search(line)
        name = m.group(1).strip()
        line = file.readline()
        dna = "";
        while (not (line == "") and line.find(">") == -1):
            line = line.strip()
            dna = dna + re.sub("\s+", "", line)
            line = file.readline()
        alignment[name] = dna
    return alignment

"""
    This function reads file in non-interleaved phylip format and returns
    a dictionary containing the alignment
"""     
def read_phylip(input_file):         
    file = open(input_file, 'r')
    alignment = {}
    counter = 0;
    pattern = re.compile(r"([^\s]*)\s+([^\s]*)")    
    line = file.readline()
    m = pattern.search(line)
    (num_taxa, num_sites) = (int(m.group(1).strip()), int(m.group(2).strip()))        
    line = file.readline()
    while (not (line == "" or line is None)):        
        m = pattern.search(line)
        alignment[m.group(1).strip()] = m.group(2).strip()
        assert len(m.group(2)) == num_sites
        line = file.readline()
    assert len(alignment.keys()) == num_taxa
    return alignment

"""    
    This function writes alignments out to a file
"""
def write_alignment(output_file, alignment, format="fasta"):
    handle = open(output_file, "w")
    records = alignment_to_records(alignment)
    SeqIO.write(records, output_file, format)
    handle.flush()
    handle.close()
    
def write_phylip(output_file, alignment):    
    header = False
    fid = open(output_file, 'w')
    for (key, value) in alignment.items():
        if (not header):
            fid.write("%d %d\n" % (len(alignment.items()), len(value)))
            header = True
        fid.write("%s\t%s\n" % (key, value))
    fid.flush()
    fid.close()
        
        

"""
    This function reads an hmmr_search output file
    and returns a dictionary that contains the e-values
    of the searched fragments
"""
def hmmr_read_search(hmmr_file):
    file = open(hmmr_file, 'r');
    start_reading = False
    results = {}
    
    #Group 1 (e-value) 2 (bitscore) and 9 (taxon name) contain the relevant information, other ones can be ignored unless we plan
    #to do something later
    pattern = re.compile(r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)")
    
    for line in file:
        line = line.strip()
        if (not start_reading and line.startswith("E-value") == True):
            start_reading = True
        elif (start_reading and line == ""):
            start_reading = False
            break
        elif (start_reading):             
            matches = pattern.search(line)
            if (matches is not None and matches.group(0).find("--") == -1):
                results[matches.group(9).strip()] = (float(matches.group(1).strip()), float(matches.group(2).strip()))
    return results


"""
    This function finds the e-score for a set of sequences
    input_file in fasta format    
"""
def hmmr_search(input_file, hmmr_file, output_file, elim=None, filters=True):
    executalbe = "hmmsearch" # os.path.join(LIB_PATH,"hmmsearch")
    elim_str = ""    
    filter_string = ""
    if (elim is not None):
        elim_str = "-E %s" % str(elim)        
    if (not filters):
        filter_string = "--max"
    print "%s %s --noali %s -o %s %s %s" % (executalbe, filter_string, elim_str, output_file, hmmr_file, input_file)
    os.system("%s %s --noali %s -o %s %s %s" % (executalbe, filter_string, elim_str, output_file, hmmr_file, input_file));
    
"""
    This function aligns fragments against an alignment.  input_file is the fragment file in
    fasta format, the original_file is the reference alignment in stockholm format.  The output
    is the aligned fragments against the reference alignment, in stockholm format
    input_file is in fasta format, original_file is in stockholm format
"""
def hmmr_align(input_file, hmmr_file, output_file, original_file=None, trim=None):
    executalbe = "hmmalign" #os.path.join(LIB_PATH,"hmmalign")
    options = ""
    if (original_file is not None):
        options = "--mapali %s" % original_file
    if (trim):
        options = options + " --trim";
    print("%s %s -o %s %s %s" % (executalbe, options, output_file, hmmr_file, input_file));
    os.system("%s %s -o %s %s %s" % (executalbe, options, output_file, hmmr_file, input_file));
    
"""
    This function generates hmmr profile for an alignment
    input_file in fasta format
    
    NOTE --enone to turn "entropy weighting" off in hmmbuild step, suggested by Eddy
    Or $format = "--informat afa";
"""
def hmmr_profile(input_file, output_file, format="fasta", options=""):
    executable =  "hmmbuild" # os.path.join(LIB_PATH,"hmmbuild")
    special = ""
    if (format == "fasta"):
        special = "--informat afa"
    os.system("%s --symfrac 0.0 --dna %s %s %s %s" % (executable, special, options, output_file, input_file))


"""
   This function converts dictionary to a record class usable by BioPython
"""
def alignment_to_records(alignment):
    return [SeqRecord(Seq(alignment[key]), id=key, description="") for key in alignment.keys()]    
    

def temp_filter(item):
    if (item.find("_") == -1):
        return True    
        
def remove_all_gap_columns(alignment):
    w = len(alignment.values()[0])
    gapcols = set(range(0, w))
    for line in alignment.values():
        toremove = []
        for p in gapcols:
            if line[p] != "-":
                toremove.append(p)
        gapcols = gapcols ^ set(toremove)    
    
    included = [x for x in range(0, w) if x not in gapcols]
    
    for s, line in alignment.items():        
        l = [line[i] for i in included]
        l = array.array('c', l).tostring()
        alignment[s] = l
         
"""
    This class implements the local alignment of the fragments and generates temp files containing
    the local alignments, profiles, etc.... In order to clean the files afterward, clean_files needs
    to be called after the class runs.  The temp files will be overwritten if the functions are ran
    twice in a row.  Should break this up into more managable chunks.
    TODO:  Come up with a better name
    TODO:  Build database of profiles rather than iterate through each one individually
    TODO:  Clean up code
    TODO:  Parallelize the code, should be easy to do
    TODO:  The temporary files thing is very messy, should do something about it
    TODO:  Probably change the tempfile dictionary keys to be prefaced with the prefix so we don't ever leave
        hanging files
"""      
class NammyClass:
    def __init__(self, tree_map, alignment, fragments, prefix=None, temp_directory=None, logger=None):
        self.tree_map = tree_map
        self.alignment = alignment
        self.fragments = fragments
        self.logger = sys.stderr
        self.tempfiles = {}
        self.prefix = prefix
        self.temp_directory = temp_directory
        self.logger = logger
        
        if (self.prefix is None):
            tempfile.tempdir = self.temp_directory
            f = tempfile.NamedTemporaryFile()     
            self.prefix = f.name         
            f.close()
        
    def find_membership(self, filters=True, elim=".01"):
        #Write fragments to fasta files
        write_alignment(self.prefix + ".fragment.alignment", self.fragments)     
        self.tempfiles["fragment.alignment"] = (self.prefix + ".fragment.alignment")
             
        #Now keep track of max e-value, and which set gave the max e-value
        max_evalues = dict([(name, (100000000, -1)) for name in self.fragments.keys()])     
        
        #TODO: Very easy to parallelize this step
        #Iterate through all trees to generate profiles
        sTimer = os.times()     
        for tree_idx in self.tree_map.keys():
            #Get all the leaves in the tree
            subset = [i.taxon.label for i in self.tree_map[tree_idx]._tree.leaf_nodes()]
            
            #First generate induced alignment on subset           
            new_alignment = dict([(name, self.alignment[name]) for name in subset])
            write_alignment(self.prefix + ".alignment.sto." + str(tree_idx), new_alignment, format="stockholm")
            self.tempfiles["alignment.sto." + str(tree_idx)] = self.prefix + ".alignment.sto." + str(tree_idx)
            
            #Profile induced alignmented           
            hmmr_profile(self.prefix + ".alignment.sto." + str(tree_idx), self.prefix + ".profile." + str(tree_idx), format="stockholm")
            self.tempfiles["profile.%s" % str(tree_idx)] = self.prefix + ".profile." + str(tree_idx)
            
            #Compute e-values of fragments against induced profile
            hmmr_search(self.prefix + ".fragment.alignment", self.prefix + ".profile." + str(tree_idx), self.prefix + ".search." + str(tree_idx), elim=elim, filters=filters)
            self.tempfiles["search.%s" % str(tree_idx)] = (self.prefix + ".search." + str(tree_idx))
            
            #Read e-values
            evalues = hmmr_read_search(self.prefix + ".search." + str(tree_idx))
            for key in evalues.keys():
                self.logger.write("FRAGMENT %s %s %s %d\n" % (key, str(evalues[key][1]), str(evalues[key][0]), tree_idx))
                (best_value, set_idx) = max_evalues[key]
                if (best_value > evalues[key][0]):
                    max_evalues[key] = (evalues[key][0], tree_idx)
    
        eTimer = os.times()
        self.logger.write("Time: Find subset membership: %s \n" % (eTimer[4] - sTimer[4]))
    
        #print "\n".join([key + ": " + str(max_evalues[key]) for key in max_evalues.keys()])     
        
        #Now build the sequences for local alignment, first finding out which ones to align to which set
        local_sets = dict([(i, []) for i in xrange(len(self.tree_map))])          
        local_sets[-1] = []  #In case nothing aligned
        [local_sets[max_evalues[key][1]].append(key) for key in max_evalues]
             
        return local_sets          
    
    #Align the fragments according to the local sets  
    def locally_align(self, local_sets=None, global_align=True):
        initTimer = os.times()
        if (local_sets == None):
            local_sets = self.find_membership()          
          
        #Write orignal alignment stockholm files
        write_alignment(self.prefix + ".original.alignment.stockholm", self.alignment, 'stockholm')
        self.tempfiles["original.alignment.stockholm"] = self.prefix + ".original.alignment.stockholm"
             
        #Now go through and align the fragments locally to the subsets.  If global_align option
        #is selected, after fragments whose membership exists have been aligned, perform a global
        #alignment on the remaining           
        sTimer = os.times()     
        for key in local_sets.keys():                  
            if (len(local_sets[key]) == 0 or key == -1):
                continue           
            frag_alignment = dict([(frag, self.fragments[frag]) for frag in local_sets[key]])     
            write_alignment(self.prefix + ".fragments.alignment." + str(key), frag_alignment)           
            self.tempfiles["fragments.alignment." + str(key)] = (self.prefix + ".fragments.alignment." + str(key))
            
            hmmr_align(self.prefix + ".fragments.alignment." + str(key), self.prefix + ".profile." + str(key), self.prefix + ".fragments.aligned." + str(key), original_file=self.prefix + ".alignment.sto." + str(key), trim=False)
            self.tempfiles["fragments.aligned." + str(key)] = self.prefix + ".fragments.aligned." + str(key)                      
            write_alignment(self.prefix + "." + str(key) + ".ref.fasta" , frag_alignment)           
            self.tempfiles[str(key) + ".ref.fasta"] = (self.prefix + str(key) + ".ref.fasta")
    
        #If we global align, then align fragments that did not match up to any subset against the entire alignment      
        if (global_align and (len(local_sets[-1]) > 0)):
            sTimer = os.times()
            #Write entire alignment and profile entire alignment         
            hmmr_profile(self.prefix + ".original.alignment.stockholm", self.prefix + ".profile", format="stockholm")
            self.tempfiles["profile"] = self.prefix + ".profile"
          
            #Write unmatch fragment files and align to profile
            frag_alignment = dict([(frag, self.fragments[frag]) for frag in local_sets[-1]])
            write_alignment(self.prefix + ".fragments.unmatched", frag_alignment)
            self.tempfiles["fragments.unmatched"] = self.prefix + ".fragments.unmatched"
            hmmr_align(self.prefix + ".fragments.unmatched", self.prefix + ".profile", self.prefix + ".fragments.unmatched.aligned", original_file=self.prefix + ".original.alignment.stockholm", trim=False)     
            self.tempfiles["fragments.unmatched.aligned"] = self.prefix + ".fragments.unmatched.aligned"
             
            eTimer = os.times()
            self.logger.write("Time: Global alignment: %s \n" % (eTimer[4] - sTimer[4]))
    
        eTimer = os.times()
        self.logger.write("Time: Total alignment: %s \n" % (eTimer[4] - initTimer[4]))
        return
    
    def get_merged_alignments(self, local_sets, global_align=True):
        final_alignment = Alignment()     
        final_alignment.set_alignment(self.alignment)
              
        sTimer = os.times()     
        for key in local_sets.keys():                  
            if (len(local_sets[key]) == 0 or key == -1):
                continue
            if (not os.path.isfile(self.prefix + ".fragments.aligned." + str(key))):
                self.logger.write("ERROR: %s does not exist for merging! \n" % (self.prefix + ".fragments.aligned." + str(key)))
            else:
                temp_align = read_sto_alignment(self.prefix + ".fragments.aligned." + str(key))           
                final_alignment.merge_alignment_in(temp_align)
              
        if (global_align and os.path.isfile(self.prefix + ".fragments.unmatched.aligned")):
            temp_align = read_sto_alignment(self.prefix + ".fragments.unmatched.aligned")
            final_alignment.merge_alignment_in(temp_align)
             
        eTimer = os.times()
        remove_all_gap_columns(final_alignment)
        self.logger.write("Time: Global alignment: %s \n" % (eTimer[4] - sTimer[4]))
        
        return (final_alignment)
         
    def subset_align(self, local_sets, super_tree_map, local_to_super_map, output_name=None, ext="combined.sto"):
        #Generate all super tree induced alignments.  These will be merged with the local alignments
        alignment_map = {}
        for tree_idx in super_tree_map.keys():
            subset = [i.taxon.label for i in super_tree_map[tree_idx]._tree.leaf_nodes()]
            a = Alignment()
            a.set_alignment(dict([(s, self.alignment[s]) for s in subset]))
            alignment_map[tree_idx] = a

        
      
        #First generate fragment alignment against profile/induced alignment of tree set
        #Align to just the induced alignment, merge into the global set
        sTimer = os.times()     
        for key in local_sets.keys():                  
            if (len(local_sets[key]) == 0 or key == -1):
                continue
            if (not os.path.isfile(self.prefix + ".fragments.alignment." + str(key))):
                frag_alignment = dict([(frag, self.fragments[frag]) for frag in local_sets[key]])
                write_alignment(self.prefix + ".fragments.alignment." + str(key), frag_alignment)
            hmmr_align(self.prefix + ".fragments.alignment." + str(key), self.prefix + ".profile." + str(key), self.prefix + ".fragments.aligned." + str(key), original_file=self.prefix + ".alignment.sto." + str(key), trim=False)
            temp_align = read_sto_alignment(self.prefix + ".fragments.aligned." + str(key))
            temp_full_length_align = read_sto_alignment(self.prefix + ".alignment.sto." + str(key))
              
            #Quick test, see if full fragments in new alignment are all found in the correct tree, if not, report error
            full_length = temp_full_length_align.keys()
            for item in full_length:
                if (not (item in alignment_map[local_to_super_map[key]])):
                    print "%s is not found in tree %s" % (item, str(local_to_super_map[key]))
    
            self.tempfiles["fragments.alignment." + str(key)] = (self.prefix + ".fragments.alignment." + str(key))           
            self.tempfiles["fragments.aligned." + str(key)] = (self.prefix + ".fragments.aligned." + str(key))
              
            alignment_map[local_to_super_map[key]].merge_alignment_in(temp_align)             
    
        
        #Write back to file, if output name is none, use prefix
        if (output_name is None):
            output_name = self.prefix
             
        for aln_idx in alignment_map.keys():
            remove_all_gap_columns(alignment_map[aln_idx])
            write_alignment(("%s.%s.%s" % (output_name, str(aln_idx), ext)), alignment_map[aln_idx],
               format="stockholm")
            write_alignment(("%s.%s.%s" % (output_name, str(aln_idx), "ref.fasta")), alignment_map[aln_idx],
               format="fasta")            
            if (output_name is None):
                self.tempfiles[("%s.%s.%s" % (output_name, str(aln_idx), ext))] = ("%s.%s.%s" % (output_name, str(aln_idx)))
    
        eTimer = os.times()
        self.logger.write("Time: Done subset align and merge: %s \n" % (eTimer[4] - sTimer[4]))
        
    def clean_files(self):
        for t in self.tempfiles.values():           
            if (os.path.isfile(t)):
                remove_temp(t)
            
    def move_alignments(self, output, extension=None, alignment_name="fragments.aligned"):
        period = "."
        if (extension is None):
            period = ""
        for i in xrange(len(self.tree_map)):
            if (os.path.isfile(self.prefix + (".%s." % alignment_name) + str(i))):          
                shutil.move(self.prefix + (".%s." % alignment_name) + str(i), "%s.%s%s%s" % (output, str(i), period, extension))
                a = read_alignment( "%s.%s%s%s" % (output, str(i), period, extension), "stockholm" )
                write_alignment( "%s.%s%s%s" % (output, str(i), period, "ref.fasta"), a) 
        if (os.path.isfile(self.prefix + ".fragments.unmatched.aligned")):
            shutil.move(self.prefix + ".fragments.unmatched.aligned", "%s.unmatched%s%s" % (output, period, extension))
