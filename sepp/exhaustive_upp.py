'''
Created on Oct 10, 2012

@author: smirarab
'''
import sys,random,argparse,os,shutil
from argparse import ArgumentParser, Namespace
from sepp import get_logger
from sepp.alignment import MutableAlignment, ExtendedAlignment,_write_fasta
from sepp.exhaustive import JoinAlignJobs, ExhaustiveAlgorithm
from sepp.jobs import PastaAlignJob
from sepp.filemgr import get_temp_file
from sepp.config import options,valid_decomp_strategy
import sepp.config
from sepp.math_utils import lcm
from sepp.problem import SeppProblem
from sepp.scheduler import JobPool
from multiprocessing import Pool, Manager
from sepp.alignment import ExtendedAlignment

_LOG = get_logger(__name__)

class UPPJoinAlignJobs(JoinAlignJobs):
    '''
    After all alignments jobs for a placement subset have finished, 
    we need to build those extended alignments. This join takes care of that step. 
    '''
    def __init__(self):
        JoinAlignJobs.__init__(self)
    
    def perform(self):            
        pp = self.placement_problem        
        
        assert isinstance(pp, SeppProblem)
        pp.annotations["search_join_object"] = self                    

# Useful for multi-core merging if ever needed
#def mergetwo(x):
#    ((i,j),extended) = x
#    a=extended[i]
#    b=extended[j]
#    a.merge_in(b,convert_to_string=True)
#    extended[j] = None
#    extended[i] = a
#    del b
#    return "Success"
 

class UPPExhaustiveAlgorithm(ExhaustiveAlgorithm):
    '''
    This implements the exhaustive algorithm where all alignments subsets
    are searched for every fragment. This is for UPP, meaning that no placement
    is performed, and that there is always only one placement subset (currently).
    '''
    def __init__(self):
        ExhaustiveAlgorithm.__init__(self)     
               
    def generate_backbone(self):
        _LOG.info("Reading input sequences: %s" %(self.options.sequence_file))
        sequences = MutableAlignment()
        sequences.read_file_object(self.options.sequence_file)
        fragments = MutableAlignment()
        if (options().median_full_length is not None):
          if (options().median_full_length == -1):
            seq_lengths = sorted([len(seq) for seq in sequences.values()])              
            lengths = len(seq_lengths)
            if lengths % 2:
              options().median_full_length = (seq_lengths[lengths / 2] + seq_lengths[lengths / 2 - 1]) / 2.0
            else:
              options().median_full_length = seq_lengths[lengths / 2]              
            
          (min_length,max_length) = (int(options().median_full_length*(1-options().backbone_threshold)),int(options().median_full_length*(1+options().backbone_threshold)))
          frag_names = [name for name in sequences if len(sequences[name]) > max_length or len(sequences[name]) < min_length]
          if (len(frag_names) > 0):
              fragments = sequences.get_hard_sub_alignment(frag_names)        
              [sequences.pop(i) for i in fragments.keys()]        
        if (options().backbone_size is None):            
            options().backbone_size = min(1000,int(sequences.get_num_taxa()))
            _LOG.info("Backbone size set to: %d" %(options().backbone_size))
        if (options().backbone_size > len(sequences.keys())):
          options().backbone_size = len(sequences.keys())
        backbone_sequences = sequences.get_hard_sub_alignment(random.sample(sequences.keys(), options().backbone_size))        
        [sequences.pop(i) for i in backbone_sequences.keys()]
                
        _LOG.info("Writing backbone set. ")
        backbone = get_temp_file("backbone", "backbone", ".fas")
        _write_fasta(backbone_sequences, backbone)
         
        _LOG.info("Generating pasta backbone alignment and tree. ")
        pastaalignJob = PastaAlignJob()
        moleculeType = options().molecule
        if (options().molecule == 'amino'):
            moleculeType =  'protein'
        pastaalignJob.setup(backbone,options().backbone_size,self.options.outdir,moleculeType,options().cpu)
        pastaalignJob.run()
        pastaalignJob.read_results()
        
        options().placement_size = self.options.backbone_size
        options().alignment_file = open(self.options.outdir + "/pasta.fasta")
        options().tree_file = open(self.options.outdir + "/pasta.fasttree")
        _LOG.info("Backbone alignment written to %s.\nBackbone tree written to %s" % (options().alignment_file, options().tree_file))
        sequences.set_alignment(fragments)        
        if (len(sequences) == 0):
          _LOG.info("No query sequences to align.  Final alignment saved as %s" % self.get_output_filename("alignment.fasta"))   
          shutil.copyfile(self.options.outdir + "/pasta.fasta", self.get_output_filename("alignment.fasta"))
          sys.exit(0)
        else:
          query = get_temp_file("query", "backbone", ".fas")
          options().fragment_file = query          
          _write_fasta(sequences, query)               

    def check_options(self):
        options().info_file = "A_dummy_value"

        #Check to see if tree/alignment/fragment file provided, if not, generate it
        #from sequence file                
        if not options().tree_file is None and not options().alignment_file is None and not options().sequence_file is None:
            options().fragment_file = options().sequence_file        
        elif options().tree_file is None and options().alignment_file is None and not options().sequence_file is None:            
            self.generate_backbone()
        else:
            _LOG.error("Either specify the backbone alignment and tree and query sequences or only the query sequences.  Any other combination is invalid")
            exit(-1)
        sequences = MutableAlignment()
        sequences.read_file_object(open(self.options.alignment_file.name))            
        backbone_size = sequences.get_num_taxa()
        if options().backbone_size is None:
            options().backbone_size = backbone_size        
        assert options().backbone_size == backbone_size, ("Backbone parameter needs to match actual size of backbone; backbone parameter:%s backbone_size:%s" 
                %(options().backbone_size, backbone_size))                    
        if options().placement_size is None:
            options().placement_size = options().backbone_size
        return ExhaustiveAlgorithm.check_options(self)
        
    def merge_results(self):
        assert len(self.root_problem.get_children()) == 1, "Currently UPP works with only one placement subset."
        '''
        Merge alignment subset extended alignments to get one extended alignment
        for current placement subset.
        '''     
        pp = self.root_problem.get_children()[0]        
        _LOG.info("Merging sub-alignments for placement problem : %s." %(pp.label))
        ''' First assign fragments to the placement problem'''
        pp.fragments = pp.parent.fragments.get_soft_sub_alignment([])
        frags = []      
        for ap in pp.get_children():
            frags.extend(ap.fragments)        
        pp.fragments.seq_names.update(frags)   
        ''' Then Build an extended alignment by merging all hmmalign results''' 
        extendedAlignment = ExtendedAlignment(pp.fragments.seq_names)
        for ap in pp.children:
            assert isinstance(ap, SeppProblem)
            ''' Get all fragment chunk alignments for this alignment subset'''
            aligned_files = [fp.get_job_result_by_name('hmmalign') for 
                                fp in ap.children if 
                                fp.get_job_result_by_name('hmmalign') is not None]
            _LOG.info("Merging fragment chunks for subalignment : %s." %(ap.label))
            ap_alg = ap.read_extendend_alignment_and_relabel_columns\
                        (ap.jobs["hmmbuild"].infile , aligned_files)
            _LOG.info("Merging alignment subset into placement subset: %s." %(ap.label))
            extendedAlignment.merge_in(ap_alg,convert_to_string=False)
        
        extendedAlignment.from_bytearray_to_string()
        self.results = extendedAlignment

# Useful for multi-core merging if ever needed
#    def parallel_merge_results(self):
#        assert len(self.root_problem.get_children()) == 1, "Currently UPP works with only one placement subset."
#        '''
#        Merge alignment subset extended alignments to get one extended alignment
#        for current placement subset.
#        '''     
#        pp = self.root_problem.get_children()[0]        
#        _LOG.info("Merging sub-alignments for placement problem : %s." %(pp.label))       
#        ''' Then Build an extended alignment by merging all hmmalign results'''
#        manager = Manager() 
#        extendedAlignments = manager.list()        
#        for ap in pp.children:
#            assert isinstance(ap, SeppProblem)
#            ''' Get all fragment chunk alignments for this alignment subset'''
#            aligned_files = [fp.get_job_result_by_name('hmmalign') for 
#                                fp in ap.children if 
#                                fp.get_job_result_by_name('hmmalign') is not None]
#            _LOG.info("Merging fragment chunks for subalignment : %s." %(ap.label))
#            ap_alg = ap.read_extendend_alignment_and_relabel_columns\
#                        (ap.jobs["hmmbuild"].infile , aligned_files)
#            _LOG.info("Merging alignment subset into placement subset: %s." %(ap.label))
#            extendedAlignments.append(ap_alg) 
#            
#        while len(extendedAlignments)>1:     
#            a=range(0,len(extendedAlignments))    
#            #print [len(x) for x in extendedAlignments]
#            x = zip(a[0::2],a[1::2])
#            mapin = zip (x,[extendedAlignments]*len(x))         
#            _LOG.debug("One round of merging started. Currently have %d alignments left. " %len(extendedAlignments)) 
#            Pool(max(12,len(extendedAlignments))).map(mergetwo,mapin)
#            #print [len(x) if x is not None else "None" for x in extendedAlignments]
#            extendedAlignments = manager.list([x for x in extendedAlignments if x is not None])
#            extendedAlignments.reverse()            
#            _LOG.debug("One round of merging finished. Still have %d alignments left. " %len(extendedAlignments)) 
#        extendedAlignment = extendedAlignments[0] 
#        extendedAlignment.from_bytearray_to_string()
#        self.results = extendedAlignment
        

    def output_results(self):        
        extended_alignment = self.results        
        _LOG.info("Generating output. ")
        outfilename = self.get_output_filename("alignment.fasta")
        extended_alignment.write_to_path(outfilename)
        _LOG.info("Unmasked alignment written to %s" %outfilename)
        outfilename = self.get_output_filename("insertion_columns.txt")
        extended_alignment.write_insertion_column_indexes(outfilename)
        _LOG.info("The index of insertion columns written to %s" %outfilename)
        extended_alignment.remove_insertion_columns()
        outfilename = self.get_output_filename("alignment_masked.fasta")
        extended_alignment.write_to_path(outfilename)
        _LOG.info("Masked alignment written to %s" %outfilename)
        
    def check_and_set_sizes(self, total):
        assert (self.options.placement_size is None) or (
                self.options.placement_size >= total), \
                "currently UPP works with only one placement subset. Please leave placement subset size option blank."
        ExhaustiveAlgorithm.check_and_set_sizes(self, total)
        self.options.placement_size = total
    
    def _get_new_Join_Align_Job(self):        
        return UPPJoinAlignJobs()
    
    def modify_tree(self,a_tree):
        ''' Filter out taxa on long branches '''
        self.filtered_taxa=[]                              
        if self.options.long_branch_filter is not None:
            tr = a_tree.get_tree()
            elen = {}
            for e in tr.leaf_edge_iter():
                elen[e] = e.length
            elensort = sorted(elen.values())
            mid = elensort[len(elensort)/2]
            torem = []
            for k,v in elen.items():
                if v > mid * self.options.long_branch_filter:
                    self.filtered_taxa.append(k.head_node.taxon.label)
                    torem.append(k.head_node.taxon)
            tr.prune_taxa(torem)
            
    def create_fragment_files(self):
        alg_subset_count = len(list(self.root_problem.iter_leaves()))
        frag_chunk_count = lcm(alg_subset_count,self.options.cpu)//alg_subset_count
        _LOG.info("%d taxa pruned from backbone and added to fragments: %s" %(len(self.filtered_taxa), " , ".join(self.filtered_taxa)))        
        return self.read_and_divide_fragments(frag_chunk_count, extra_frags =\
                       self.root_problem.subalignment.get_soft_sub_alignment(\
                                                         self.filtered_taxa))
                
def augment_parser():
    sepp.config.set_main_config_path(os.path.expanduser("~/.sepp/upp.config"))
    parser = sepp.config.get_parser()    
    parser.description = "This script runs the UPP algorithm on set of sequences.  A backbone alignment and tree can be given as input.  If none is provided, a backbone will be automatically generated."
    
    decompGroup = parser.groups['decompGroup']                                 
    decompGroup.__dict__['description'] = ' '.join(["These options",
        "determine the alignment decomposition size, backbone size, and how to decompose the backbone set."])
        
    
    decompGroup.add_argument("-A", "--alignmentSize", type = int, 
                      dest = "alignment_size", metavar = "N", 
                      default = 10,
                      help = "max alignment subset size of N "
                             "[default: 10]")    
    decompGroup.add_argument("-M", "--median_full_length", type = int, 
                      dest = "median_full_length", metavar = "N", 
                      default = None,
                      help = "Consider all fragments that are 25%% longer or shorter than N to be excluded from the backbone.  If value is -1, then UPP will use the median of the sequences as the median full length "
                             "[default: None]")                                 
    decompGroup.add_argument("-T", "--backbone_threshold", type = float, 
                      dest = "backbone_threshold", metavar = "N", 
                      default = 0.25,
                      help = "Only consider sequences with lengths within (1-N)*M and (1+N)*M as full-length, where M is the median length of the full-length sequence given by the -M option and N is the percentage given by the -T option."
                             "[default: 0.25]")                              
    decompGroup.add_argument("-B", "--backboneSize", type = int,
                      dest = "backbone_size", metavar = "N", 
                      default = None,
                      help = "(Optional) size of backbone set.  "
                             "If no backbone tree and alignment is given, the sequence file will be randomly split into a backbone set (size set to N) and query set (remaining sequences), [default: min(1000,input size)]")    
    decompGroup.add_argument("-S", "--decomp_strategy", type = valid_decomp_strategy, 
                      dest = "decomp_strategy", metavar = "DECOMP",
                      default = "hierarchical", 
                      help = "decomposition strategy "
                             "[default: ensemble of HMMs (hierarchical)]")                              
                             
    inputGroup = parser.groups['inputGroup']                             
    inputGroup .add_argument("-s", "--sequence_file", type = argparse.FileType('r'),
                      dest = "sequence_file", metavar = "SEQ", 
                      default = None,
                      help = "Unaligned sequence file.  "
                             "If no backbone tree and alignment is given, the sequence file will be randomly split into a backbone set (size set to B) and query set (remaining sequences), [default: None]")                             
    inputGroup.add_argument("-c", "--config", 
                      dest = "config_file", metavar = "CONFIG",
                      type = argparse.FileType('r'), 
                      help = "A config file, including options used to run UPP. Options provided as command line arguments overwrite config file values for those options. "
                             "[default: %(default)s]")    
    inputGroup.add_argument("-t", "--tree", 
                      dest = "tree_file", metavar = "TREE",
                      type = argparse.FileType('r'), 
                      help = "Input tree file (newick format) "
                             "[default: %(default)s]")    
    inputGroup.add_argument("-a", "--alignment", 
                      dest = "alignment_file", metavar = "ALIGN",
                      type = argparse.FileType('r'), 
                      help = "Aligned fasta file "
                             "[default: %(default)s]")                                 
                             
    uppGroup = parser.add_argument_group("UPP Options".upper(), 
                         "These options set settings specific to UPP")                                 
    
    uppGroup.add_argument("-l", "--longbranchfilter", type = int, 
                      dest = "long_branch_filter", metavar = "N", 
                      default = None,
                      help = "Branches longer than N times the median branch length are filtered from backbone and added to fragments."
                             " [default: None (no filtering)]")
                             
    seppGroup = parser.add_argument_group("SEPP Options".upper(), 
                         "These options set settings specific to SEPP and are not used for UPP.")                                 
    seppGroup.add_argument("-P", "--placementSize", type = int, 
                      dest = "placement_size", metavar = "N",
                      default = None, 
                      help = "max placement subset size of N "
                             "[default: 10%% of the total number of taxa]")                              
    seppGroup.add_argument("-r", "--raxml", 
                      dest = "info_file", metavar = "RAXML",
                      type = argparse.FileType('r'), 
                      help = "RAxML_info file including model parameters, generated by RAxML."
                             "[default: %(default)s]")    
    seppGroup.add_argument("-f", "--fragment",
                      dest = "fragment_file", metavar = "FRAG",
                      type = argparse.FileType('r'), 
                      help = "fragment file "
                             "[default: %(default)s]")          
                             
                                                   
def main():
    augment_parser() 
    UPPExhaustiveAlgorithm().run()

if __name__ == '__main__':   
    main()
        
