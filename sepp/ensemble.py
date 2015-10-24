'''
Created on May 31, 2015

@author: namphuon
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
from sepp.scheduler import JobPool,Join
from multiprocessing import Pool, Manager
from sepp.alignment import ExtendedAlignment

_LOG = get_logger(__name__)
 
class EnsembleJoinSearchJobs(Join):
    '''
    After all search jobs have finished on tips, we need return the distribution
    of the bitscores for the search.  This join accomplishes this
    '''
    def __init__(self):
        Join.__init__(self)
        
    def setup_with_root_problem(self, root_problem):
        self.root_problem = root_problem            
        for p in root_problem.iter_leaves():
            self.add_job(p.jobs["hmmsearch"])

    def perform(self):
        '''
        A dummy join that waits for all the search results to complete     
        '''
        print ""
                
                
    def __str__(self):
        return "join search jobs for all tips of ", self.root_problem

class EnsembleExhaustiveAlgorithm(ExhaustiveAlgorithm):
    '''
    This implements the exhaustive algorithm where all alignments subsets
    are searched for every fragment. This is for UPP, meaning that no placement
    is performed, and that there is always only one placement subset (currently).
    '''
    def __init__(self):
        ExhaustiveAlgorithm.__init__(self)    
        self.symfrac = False
        self.elim = None
        self.filters = True 
               
    def check_options(self):
        options().info_file = "A_dummy_value"

        if options().tree_file is None or options().alignment_file is None:
            _LOG.error("Specify the backbone alignment and tree and query sequences")
            exit(-1)
        sequences = MutableAlignment()
        sequences.read_file_object(open(self.options.alignment_file.name))  
        return ExhaustiveAlgorithm.check_options(self)
            
    def check_and_set_sizes(self, total):
        assert (self.options.placement_size is None) or (
                self.options.placement_size >= total), \
                "currently eHMMs works with only one placement subset. Please leave placement subset size option blank."        
        ExhaustiveAlgorithm.check_and_set_sizes(self, total)
        self.options.placement_size = total
        
        
    def merge_results(self):
        ''' merges search results'''
        if self.root_problem.annotations.has_key("fragments.distribution.done"):
            return
        sequence_scores = dict([(name, []) for name in self.root_problem.fragments.keys()])        
        for fragment_chunk_problem in self.root_problem.iter_leaves():
            align_problem = fragment_chunk_problem.get_parent()
            assert isinstance(align_problem, SeppProblem)
            '''For each subproblem start with an empty set of fragments, 
            and add to them as we encounter new best hits for that subproblem'''
            if align_problem.fragments is None: 
                align_problem.fragments = self.root_problem.fragments.get_soft_sub_alignment([])
            search_res = fragment_chunk_problem.get_job_result_by_name("hmmsearch")
            for key in search_res.keys():
                 sequence_scores[key].append([search_res[key][1],search_res[key][0]])

                    
        # TODO: is the following efficient enough? Do we need to make lists
        # and then turn them to sets?
        notScored = []
        for key,v in sequence_scores.iteritems():
            if len(v) == 0:
                notScored.append(key)
                    
        self.root_problem.annotations["fragments.distribution.done"] = 1

        ''' Make sure all fragments are in at least one subproblem. 
        TODO: what to do with those that are not?  For now, only output warning message'''
        #notScored = [k for k, v in max_evalues.iteritems() if v[1] is None]
        _LOG.warning("Fragments %s are not scored against any subset" %str(notScored))
        #assert len(notScored) == 0, "Fragments %s are not scored against any subset" %str(notScored)  
        self.results = sequence_scores      
        
    def connect_jobs(self):
        ''' a callback function called after hmmbuild jobs are finished'''
        def enq_job_searchfragment(result, search_job):
            search_job.hmmmodel = result 
            JobPool().enqueue_job(search_job)        
        assert isinstance(self.root_problem, SeppProblem)
        for placement_problem in self.root_problem.get_children():
            '''For each alignment subproblem, ...'''
            for alg_problem in placement_problem.children:
                assert isinstance(alg_problem, SeppProblem)                
                ''' create the build model job'''
                bj = alg_problem.jobs["hmmbuild"]            
                ''' create the search jobs'''
                for fc_problem in alg_problem.get_children():
                    sj = fc_problem.jobs["hmmsearch"]
                    ''' connect bulid and search jobs'''
                    bj.add_call_Back(lambda result, next_job = sj: enq_job_searchfragment(result, next_job))
        jsj = EnsembleJoinSearchJobs()
        jsj.setup_with_root_problem(self.root_problem)         
                
    def output_results(self):        
        search_results = self.results        
        _LOG.info("Generating csv of search results. ")
        outfilename = self.get_output_filename("scores.csv")
        not_matched = self.get_output_filename("unmatched.csv")
        f = open(outfilename, 'w')
        unmatched = open(not_matched, 'w')
        f.write("seq,bitscore,evalue\n")
        
        for key,value in search_results.items():
            if len(value) == 0:
                unmatched.write("%s " % key)
            else:
                for pair in value:
                    f.write("%s,%0.4f,%s\n" % (key,pair[0],"{:.3e}".format(pair[1])))
        f.close()
        unmatched.close()
                        
    def create_fragment_files(self):
        alg_subset_count = len(list(self.root_problem.iter_leaves()))
        frag_chunk_count = lcm(alg_subset_count,self.options.cpu)//alg_subset_count
        return self.read_and_divide_fragments(frag_chunk_count)

                
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
    EnsembleExhaustiveAlgorithm().run()

if __name__ == '__main__':   
    main()
        
