'''
Created on Oct 2, 2012

@author: smirarab
'''
from sepp.config import options
from abc import abstractmethod, ABCMeta
from sepp.scheduler import JobPool
from sepp.filemgr import directory_has_files_with_prefix, get_temp_file
from sepp.alignment import MutableAlignment
from sepp.tree import PhylogeneticTree
import dendropy
from sepp import get_logger
import sys
import os
from sepp.problem import SeppProblem

_LOG = get_logger(__name__)

class AbstractAlgorithm(object):
    '''
    This class provides the interface for an abstract algorithm. All different
    ways of implementing SEPP are considered to be different 'algorithms'.
    
    New algorithms should be implemented by subclassing AbastractAlgorithm and 
    implementing its abstract methods.    
    '''

    __metaclass__ = ABCMeta
    
    def __init__(self):
        '''
        Constructor
        '''
        self.root_problem = None
        self.results = None
        self.options = options() # for ease of access
        pass
    
    
    def check_options(self):
        '''
        This method should check the input values stored in config.option to make
        sure every necessary argument is provided, and that the provided values
        are all fine. 
        
        In the event of recognizing invalid or missing input, a an Exception
        (maybe an ArgumentError) should be raised. 
        
        By default expects tree_file, raxml_file, and fragment_file. Overwrite if required. 
        '''
        supply = []
        if (options().tree_file is None):
            supply = supply + ["tree file"]
        if (options().alignment_file is None):
            supply = supply + ["alignment file"]
        if (options().info_file is None):
            supply = supply + ["raxml file"];
        if (options().fragment_file is None):
            supply = supply + ["fragment file"]
        if (len(supply) != 0):
            raise ValueError ("Failed to supply: %s\nRun with -h option to see a list of options." % " , ".join(supply))
        
        self.check_outputprefix()
        pass

    
    @abstractmethod
    def build_subproblems(self):
        '''
        This method should read the config.options() and build a hierarchy of
        subproblems (see sepp.problem). This hierarchy will be used by 
        build_job_dag to create a DAG of jobs (see sepp.jobs and sepp.scheduler)        
        '''
        raise NotImplementedError()
    
    @abstractmethod
    def build_job_dag(self):
        '''
        Builds separate jobs for different tasks that need to be done. 
        Jobs are joined together to form a DAG using Joins (see sepp.scheduler)
        and call_back functions (see sepp.scheduler.Job). 

        Once the first level of jobs (those with no dependency) are enqueued
        (using enqueue_firstlevel_job) everything else should be automatically 
        enqueued when all their dependencies are satisfied. 
        For example, if the conceptual DAG is:
        A----> E---|--> F
                   |
        B--|-> D---|
        C--| 
        then enqueue_firstlevel_job should enqueue A B and C. When A finishes,
        it should enqueue E (either using a callback or a Join with only one
        dependency). B and C should be joined together using a Join, and that 
        Join needs to enqueue D. Also E and D need to be joined, and their join
        needs to enqueue F.          
        '''
        raise NotImplementedError()

    @abstractmethod
    def enqueue_firstlevel_job(self):
        '''
        This is called after the DAG is created (see buld_job_dag) to enqueue
        the first level of jobs (those with no dependency). These jobs should
        automatically enqueue the rest of the DAG upon completion.
        '''
        raise NotImplementedError()    

    @abstractmethod
    def merge_results(self):
        '''
        This method should read results and prepare all the final results
        for output. Resutls should be saved in the problem hierarchy, and 
        by setting ?? attribute        
        '''
        raise NotImplementedError()

    @abstractmethod
    def output_results(self):
        '''
        This method should output the final results. Final results can be found
        in the problem hierarchy        
        '''
        raise NotImplementedError()
            
    def run(self):
        
        '''check input arguments'''
        self.check_options()
        
        '''build the problem structure'''
        self.root_problem = self.build_subproblems()          
        
        '''build a DAG for running all jobs'''
        self.build_job_dag()
        
        '''Queu up first level jobs (i.e. those no dependency).
        Once these run, they should automatically enqueue the rest of the
        DAG through joins and callbacks '''                            
        self.enqueue_firstlevel_job()
        
        '''Wait for all jobs to finish'''
        if (not JobPool().wait_for_all_jobs()):
            _LOG.exception("There have been errors in executed jobs. Terminating.")
            sys.exit(1)    
        
        '''Merge results into final outputs'''
        self.merge_results()
        
        '''Output final results'''
        self.output_results()         
        

    ''' The following are a bunch of helper methods that will be needed in 
    most implementations of sepp'''
        
    def check_and_set_sizes(self, total):
        #If sizes are not set, then use 10% rule
        options = self.options
        if (options.alignment_size is None):
            options.alignment_size = int(total*.10)
        if (options.placement_size is None):
            options.placement_size = options.alignment_size
        if options.placement_size is not None and options.placement_size < options.alignment_size:
            raise ValueError("alignment_size (%d) cannot be larger than placement_size (%d). " %(options.alignment_size,options.placement_size))                    
        if options.placement_size > total or options.placement_size > total:
            raise ValueError("alignment_size (%d) and placement_size (%d) should be smaller than total number of taxa (%d). " %(options.alignment_size,options.placement_size,total))

    def check_outputprefix(self):
        if directory_has_files_with_prefix(self.options.outdir,self.options.output):
            raise ValueError("Output directory [%s] already contains files with prefix [%s]...\nTerminating to avoid loss of existing files." % (self.options.outdir,self.options.output))

    def get_output_filename(self, name):
        return os.path.join(self.options.outdir,"%s_%s" %(self.options.output,name))
    
    def read_alignment_and_tree(self):
        alignment = MutableAlignment()
        alignment.read_file_object(self.options.alignment_file)
        
        #fragments = MutableAlignment()
        #fragments.read_file_object(self.options.fragment_file);   
        
        tree = PhylogeneticTree( dendropy.Tree(stream=self.options.tree_file, 
                                               schema="newick", 
                                               preserve_underscores=True))        
        
        return (alignment, tree)

    def read_and_divide_fragments(self, chunks):
        alignment = MutableAlignment()
        alignment.read_file_object(self.options.fragment_file)
        
        chunksize = alignment.get_num_taxa()//chunks + 1;
        names = alignment.get_sequence_names()
        
        ret = []
        for i in xrange(0,chunks):
            temp_file = get_temp_file("fragment_chunk_%d" %i, "fragment_chunks", ".fasta")
            subset = names[i*chunksize:max((i+1)*chunksize,len(names))]
            subset_alg = alignment.get_soft_sub_alignment(subset)
            subset_alg.write_to_path(temp_file)
            ret.append(temp_file)            
        return ret
    
    def _create_root_problem(self, tree, alignment):
        ''' Create the root problem'''   
        self.root_problem = SeppProblem(tree.leaf_node_names())
        self.root_problem.label = "root"
        self.root_problem.subalignment = alignment
        self.root_problem.subtree = tree