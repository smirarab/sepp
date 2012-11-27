'''
Created on Oct 10, 2012

@author: smirarab
'''
from sepp import get_logger
from sepp.exhaustive import JoinAlignJobs, ExhaustiveAlgorithm
from sepp.jobs import PplacerJob
from sepp.config import options
import sepp.config
from sepp.math_utils import lcm

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
        
        extendedAlignment = self.merge_subalignments()
        
        pj = pp.jobs["pplacer"]
        assert isinstance(pj, PplacerJob)
        pj.set_attribute("extended_alignment_object", extendedAlignment)                    

class UPPExhaustiveAlgorithm(ExhaustiveAlgorithm):
    '''
    This implements the exhaustive algorithm where all alignments subsets
    are searched for every fragment. This is for UPP, meaning that no placement
    is performed, and that there is always only one placement subset (currently).
    '''
    def __init__(self):
        ExhaustiveAlgorithm.__init__(self)     

    def check_options(self):
        options().info_file = "A_dummy_value"
        return ExhaustiveAlgorithm.check_options(self)


    def merge_results(self):
        ''' no .json files to merge'''
        pass

    def output_results(self):
        _LOG.info("Generating output. ")
        assert len(self.root_problem.get_children()) == 1, "Currently UPP works with only one placement subset."
        pp = self.root_problem.get_children()[0]
        extended_alignment = pp.jobs["pplacer"].get_attribute("extended_alignment_object")
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
        return self.read_and_divide_fragments(frag_chunk_count, extra_frags=self.filtered_taxa)
                
def augment_parser():
    parser = sepp.config.get_parser()
    uppGroup = parser.add_argument_group("UPP Options".upper(), 
                         "These options set settings specific to UPP")                                 
    
    uppGroup.add_argument("-l", "--longbranchfilter", type = int, 
                      dest = "long_branch_filter", metavar = "N", 
                      default = None,
                      help = "Branches longer than N times the median branch length are filtered from backbone and added to fragments."
                             " [default: None (no filtering)]")                            

if __name__ == '__main__':   
    augment_parser() 
    UPPExhaustiveAlgorithm().run()
