'''
Created on Oct 10, 2012

@author: smirarab
'''
import sys,random,argparse
from argparse import ArgumentParser, Namespace
from sepp import get_logger
from sepp.alignment import MutableAlignment, ExtendedAlignment,_write_fasta
from sepp.exhaustive import JoinAlignJobs, ExhaustiveAlgorithm
from sepp.jobs import PplacerJob,MafftAlignJob,FastTreeJob,SateAlignJob
from sepp.filemgr import get_temp_file
from sepp.config import options
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
        backbone_sequences = sequences.get_hard_sub_alignment(random.sample(sequences.keys(), options().backbone_size))        
        [sequences.pop(i) for i in backbone_sequences.keys()]
        
        _LOG.info("Writing query and backbone set. ")
        query = get_temp_file("backbone", "query", ".fas")
        backbone = get_temp_file("backbone", "backbone", ".fas")
        _write_fasta(sequences, query)
        _write_fasta(backbone_sequences, backbone)
        
        _LOG.info("Generating mafft backbone alignment and tree. ")
        mafftJob = MafftAlignJob()
        mafft_alignment = "/projects/sate5/ultra_large/temp/upp_100_10_size_alignment.WkuZl3/mafft/backboneZ7DgWh.fasta"#get_temp_file("backbone", "mafft", ".fasta")
        mafft_tree = '/projects/sate5/ultra_large/temp/upp_100_10_size_alignment.WkuZl3/mafft/backboneZ7DgWh.fasttree' #get_temp_file("backbone", "mafft", ".fasttree")        
        #mafftJob.setup(backbone,options().backbone_size,mafft_alignment,options().cpu)
        #mafftJob.run()
        #mafftJob.read_results()
        fasttreeJob = FastTreeJob()
        fasttreeJob.setup(mafft_alignment,mafft_tree,options().molecule)
        fasttreeJob.run()
        fasttreeJob.read_results()
        
        backbone_alignment = self.options.outdir + "/backbone.fasta"
        backbone_tree = self.options.outdir + "/backbone.fasttree"
        
        outfilename = self.get_output_filename("alignment.fasta")
        extended_alignment.write_to_path(outfilename)               

    def check_options(self):
        options().info_file = "A_dummy_value"
        
        #Check to see if tree/alignment/fragment file provided, if not, generate it
        #from sequence file        
        if (options().tree_file is None and options().alignment_file is None and options().fragment_file is None and not options().sequence_file is None):
            self.generate_backbone()
        
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
        for ap in pp.get_children():
            pp.fragments.seq_names.extend(ap.fragments)   
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
    parser = sepp.config.get_parser()
    
    uppGroup = parser.add_argument_group("UPP Options".upper(), 
                         "These options set settings specific to UPP")                                 
    
    uppGroup.add_argument("-l", "--longbranchfilter", type = int, 
                      dest = "long_branch_filter", metavar = "N", 
                      default = None,
                      help = "Branches longer than N times the median branch length are filtered from backbone and added to fragments."
                             " [default: None (no filtering)]")
                             
    uppGroup.add_argument("-s", "--sequence_file", type = argparse.FileType('r'),
                      dest = "sequence_file", metavar = "SEQ", 
                      default = None,
                      help = "(Optional) unaligned sequence file"
                             "If no backbone tree and alignment is given, the sequence file will be randomly split into a backbone set (size set to P) and query set (remaining sequences), [default: None]")

    uppGroup.add_argument("-B", "--backboneSize", type = int,
                      dest = "backbone_size", metavar = "N", 
                      default = 100,
                      help = "(Optional) size of backbone set"
                             "If no backbone tree and alignment is given, the sequence file will be randomly split into a backbone set (size set to N) and query set (remaining sequences), [default: None]")    
                             

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
                                                   
if __name__ == '__main__':   
    augment_parser() 
    UPPExhaustiveAlgorithm().run()
