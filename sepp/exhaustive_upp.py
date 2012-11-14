'''
Created on Oct 10, 2012

@author: smirarab
'''
from sepp import get_logger
from sepp.exhaustive import JoinAlignJobs, ExhaustiveAlgorithm
from sepp.jobs import PplacerJob
from sepp.scheduler import JobPool

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
        assert isinstance(pj,PplacerJob)
        pj.set_keyword_attribute("extended_alignment_object", extendedAlignment)
        #extendedAlignment.write_to_path(pj.extended_alignment_file)                       

class UPPExhaustiveAlgorithm(ExhaustiveAlgorithm):
    '''
    This implements the exhaustive algorithm where all alignments subsets
    are searched for every fragment. This is for UPP, meaning that no placement
    is performed, and that there is always only one placement subset (currently).
    '''
    def __init__(self):
        ExhaustiveAlgorithm.__init__(self)     

    def merge_results(self):
        ''' no .json files to merge'''
        pass

    def output_results(self):
        assert len(self.root_problem.get_children()) == 1, "Currently UPP works with only one placement subset."
        pp = self.root_problem.get_children()[0]
        extended_alignment = pp.jobs["pplacer"].get_keyword_attribute("extended_alignment_object")
        outfilename = self.get_output_filename("alignment.fasta")
        extended_alignment.write_to_path(outfilename)
        masked = extended_alignment.get_insertion_masked_alignment()
        outfilename = self.get_output_filename("masked_alignment.fasta")
        masked.write_to_path(outfilename)
        
    def check_and_set_sizes(self, total):
        assert (self.options.placement_size is None) or (self.options.placement_size >= total), "currently UPP works with only one placement subset"
        ExhaustiveAlgorithm.check_and_set_sizes(self, total)
        self.options.placement_size = total
    
    def _get_new_Join_Align_Job(self):        
        return UPPJoinAlignJobs()

if __name__ == '__main__':
    UPPExhaustiveAlgorithm().run()