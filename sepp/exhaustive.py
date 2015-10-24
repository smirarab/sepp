'''
Created on Oct 10, 2012

@author: smirarab
'''
from sepp.algorithm import AbstractAlgorithm
from sepp.config import options
from sepp.tree import PhylogeneticTree
from sepp.alignment import MutableAlignment, ExtendedAlignment, hamming_distance
from sepp.problem import SeppProblem
from dendropy.dataobject.tree import Tree
from sepp.jobs import HMMBuildJob, HMMSearchJob, HMMAlignJob, PplacerJob,\
    MergeJsonJob
from sepp.scheduler import JobPool, Join
from sepp import get_logger
from sepp.math_utils import lcm
import pdb
_LOG = get_logger(__name__)

class JoinSearchJobs(Join):
    '''
    After all search jobs have finished on tips, we need to figure out which 
    fragment goes  to which subset and start aligning fragments. 
    This join takes care of that step. 
    '''
    def __init__(self):
        Join.__init__(self)
        
    def setup_with_root_problem(self, root_problem):
        self.root_problem = root_problem            
        for p in root_problem.iter_leaves():
            self.add_job(p.jobs["hmmsearch"])         
    
    def figureout_fragment_subset(self):
        ''' Figure out which fragment should go to which subproblem'''
        if self.root_problem.annotations.has_key("fragments.distribution.done"):
            return
        max_evalues = dict([(name, (None, None)) for name in self.root_problem.fragments.keys()])
        for fragment_chunk_problem in self.root_problem.iter_leaves():
            align_problem = fragment_chunk_problem.get_parent()
            assert isinstance(align_problem, SeppProblem)
            '''For each subproblem start with an empty set of fragments, 
            and add to them as we encounter new best hits for that subproblem'''
            if align_problem.fragments is None: 
                align_problem.fragments = self.root_problem.fragments.get_soft_sub_alignment([])
            search_res = fragment_chunk_problem.get_job_result_by_name("hmmsearch")
            for key in search_res.keys():
                (best_value, prev_align_problem) = max_evalues[key]
                ''' If this is better than previous best hit, remove this
                fragment from the previous hit, and add it to this subproblem 
                '''
                if best_value is None or (best_value < search_res[key][1]):
                    max_evalues[key] = (search_res[key][1], align_problem)
                    
        # TODO: is the following efficient enough? Do we need to make lists
        # and then turn them to sets?
        notScored = []
        for key,v in max_evalues.iteritems():
            if v[1] is None:
                notScored.append(key)
            else:
                v[1].fragments.seq_names.add(key)
                    
        self.root_problem.annotations["fragments.distribution.done"] = 1

        ''' Make sure all fragments are in at least one subproblem. 
        TODO: what to do with those that are not?  For now, only output warning message'''
        #notScored = [k for k, v in max_evalues.iteritems() if v[1] is None]
        _LOG.warning("Fragments %s are not scored against any subset" %str(notScored))
        #assert len(notScored) == 0, "Fragments %s are not scored against any subset" %str(notScored)

    def perform(self):
        '''
        Distributes fragments to alignments subsets with best score, 
        and runs align jobs on those. Also, creates new chunks of fragments
        for better parallelism.     
        '''
  
        ''' Figure out which fragment should go to which subproblem'''
        self.figureout_fragment_subset()
                
        ''' For each alignment subproblem, 
        1) make sure its fragments are evenly distributed to fragment chunks. 
        2) Setup alignment jobs for its children and enqueue them'''
        alg_problems = [alg for p in self.root_problem.children for alg in p.children ]
        for alg_problem in alg_problems:
            assert isinstance(alg_problem, SeppProblem)
            chunks = len(alg_problem.get_children())
            fragment_chunks = alg_problem.fragments.divide_to_equal_chunks(chunks)

            ''' Now setup alignment jobs and enqueue them'''
            for (i,fragment_chunk_problem) in enumerate(alg_problem.children):        
                fragment_chunk_problem.fragments = fragment_chunks[i] 
                aj = fragment_chunk_problem.jobs['hmmalign']
                assert isinstance(aj,HMMAlignJob)            
                ''' First Complete setting up alignments''' 
                aj.hmmmodel = alg_problem.get_job_result_by_name('hmmbuild')
                aj.base_alignment = alg_problem.jobs["hmmbuild"].infile    

                if fragment_chunk_problem.fragments is None or fragment_chunk_problem.fragments.is_empty():
                    aj.fake_run = True
                else:
                    fragment_chunk_problem.fragments.write_to_path(aj.fragments)
                ''' Now the align job can be put on the queue '''
                JobPool().enqueue_job(aj)                
                
    def __str__(self):
        return "join search jobs for all tips of ", self.root_problem

class JoinAlignJobs(Join):
    '''
    After all alignments jobs for a placement subset have finished, 
    we need to build those extended alignments and start placing fragments. 
    This join takes care of that step. 
    '''
    def __init__(self):
        Join.__init__(self)
        
    def setup_with_placement_problem(self, placement_problem):
        self.placement_problem = placement_problem            
        for p in placement_problem.iter_leaves():
            self.add_job(p.jobs["hmmalign"])         
    
    def merge_subalignments(self):   
        '''
        Merge alignment subset extended alignments to get one extended alignment
        for current placement subset.
        '''     
        pp = self.placement_problem
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
            del ap_alg
        
        extendedAlignment.from_bytearray_to_string()
        return extendedAlignment
    
    def perform(self):            
        pp = self.placement_problem

        fullExtendedAlignment = self.merge_subalignments()
        
        #Split the backbone alignment and query sequences into separate files        
        queryExtendedAlignment = fullExtendedAlignment.get_fragments_readonly_alignment()
        baseAlignment = fullExtendedAlignment.get_base_readonly_alignment()

        pj = pp.jobs["pplacer"]
        assert isinstance(pj,PplacerJob)
        if (queryExtendedAlignment.is_empty()):
          pj.fake_run = True
        
        #Write out the extended alignments, split into query and full-length for pplacer
        queryExtendedAlignment.write_to_path(pj.extended_alignment_file)          
        baseAlignment.write_to_path(pj.backbone_alignment_file)
        
        #But keep the extended alignment on everything 
        pj.set_attribute("full_extended_alignment_object", fullExtendedAlignment)

        JobPool().enqueue_job(pj)

    def __str__(self):
        return "join align jobs for tips of ", self.placement_problem

class ExhaustiveAlgorithm(AbstractAlgorithm):
    '''
    This implements the exhaustive algorithm where all alignments subsets
    are searched for every fragment. 
    '''
    def __init__(self):
        AbstractAlgorithm.__init__(self)
        self.place_nomatch_fragments = False
        ''' Hardcoded E-Lim for hmmsearch ''' #TODO: what to do with this
        self.elim = 99999999
        self.filters = False
        self.symfrac = True
        self.strategy = options().exhaustive.strategy
        self.minsubsetsize = int(options().exhaustive.minsubsetsize)
        #Temp fix for now, 
        self.molecule = self.options.molecule
        self.distances = dict()

    def compute_distances(self, sequences):
        for seq1 in sequences.keys():
            for seq2 in sequences.keys():
                if ("".join([seq1,seq2]) not in self.distances):
                    self.distances["".join([seq1,seq2])] = hamming_distance(sequences[seq1],sequences[seq2])
                    self.distances["".join([seq2,seq1])] = self.distances["".join([seq1,seq2])]

    def merge_results(self):
        ''' TODO: implement this 
        '''        
        assert isinstance(self.root_problem,SeppProblem)
        
        '''Generate single extended alignment'''
        fullExtendedAlignment = self.root_problem.get_children()[0].jobs["pplacer"].get_attribute("full_extended_alignment_object")
        for pp in self.root_problem.get_children()[1:]:
            extended_alignment = pp.jobs["pplacer"].get_attribute("full_extended_alignment_object")
            fullExtendedAlignment.merge_in(extended_alignment,convert_to_string=True)
        self.results = fullExtendedAlignment
        
        mergeinput = []
        '''Append main tree to merge input'''
        mergeinput.append("%s;" %(self.root_problem.subtree.compose_newick(labels = True)))
        jsons = []
        for pp in self.root_problem.get_children():
            assert isinstance(pp,SeppProblem)
            if (pp.get_job_result_by_name("pplacer") is None):
              continue
            '''Append subset trees and json locations to merge input'''
            mergeinput.append("%s;\n%s" %(pp.subtree.compose_newick(labels = True),
                              pp.get_job_result_by_name("pplacer")))
        mergeinput.append("")
        mergeinput.append("")
        meregeinputstring = "\n".join(mergeinput)
        mergeJsonJob = MergeJsonJob()
        mergeJsonJob.setup(meregeinputstring, 
                           self.get_output_filename("placement.json"))
        mergeJsonJob.run()

    def output_results(self):
        ''' Merged json file is already saved in merge_results function and
            full extended alignment already created in merge_results function
        '''
        outfilename = self.get_output_filename("alignment.fasta")
        self.results .write_to_path(outfilename)
        self.results.remove_insertion_columns()
        outfilename = self.get_output_filename("alignment_masked.fasta")
        self.results.write_to_path(outfilename)

    def check_options(self, supply = []):
        AbstractAlgorithm.check_options(self,supply)
            

    def modify_tree(self,a_tree):
        pass

    def build_subproblems(self):
        (alignment, tree) = self.read_alignment_and_tree()
        
        if options().distance != 1:
            self.compute_distances(alignment)
        
        
        assert isinstance(tree, PhylogeneticTree)
        assert isinstance(alignment, MutableAlignment)

        tree.get_tree().resolve_polytomies()
        # Label edges with numbers so that we could assemble things back
        # at the end
        tree.lable_edges()        

        ''' Make sure size values are set, and are meaningful. '''
        self.check_and_set_sizes(alignment.get_num_taxa())        

        self._create_root_problem(tree, alignment)

        ''' Decompose the tree based on placement subsets'''
        placement_tree_map = PhylogeneticTree(Tree(tree.den_tree)).decompose_tree(
                                        self.options.placement_size, 
                                        strategy=self.strategy, 
                                        minSize = self.minsubsetsize,
                                        tree_map = {},pdistance = 1,distances = self.distances)
        assert len(placement_tree_map) > 0, ("Tree could not be decomposed"
                " given the following settings; strategy:%s minsubsetsize:%s placement_size:%s" 
                %(self.strategy, self.minsubsetsize, self.options.placement_size))                    
        _LOG.info("Breaking into %d placement subsets." %len(placement_tree_map))

        ''' For placement subsets create a placement subproblem, and decompose further'''
        for (p_key,p_tree) in placement_tree_map.items():
            assert isinstance(p_tree, PhylogeneticTree)
            placement_problem  = SeppProblem(p_tree.leaf_node_names(), self.root_problem)
            placement_problem.subtree = p_tree
            placement_problem.label = "P_%s" %str(p_key)
            _LOG.debug("Placement subset %s has %d nodes" %(placement_problem.label,len(p_tree.leaf_node_names())))
            ''' Further decompose to alignment subsets '''
            alignment_tree_map = PhylogeneticTree(Tree(p_tree.den_tree)).decompose_tree(
                                        self.options.alignment_size, 
                                        strategy=self.strategy, 
                                        minSize = self.minsubsetsize,
                                        tree_map = {}, decomp_strategy = self.options.decomp_strategy, pdistance = options().distance,distances = self.distances)
            assert len(alignment_tree_map) > 0, ("Tree could not be decomposed"
            " given the following settings; strategy:%s minsubsetsize:%s alignmet_size:%s" 
            %(self.strategy, self.minsubsetsize, self.options.alignment_size))
                        
            _LOG.debug("Placement subset %s has %d alignment subsets: %s" %(placement_problem.label,len(alignment_tree_map.keys()),str(sorted(alignment_tree_map.keys()))))
            _LOG.debug("Placement subset %s has %d taxa:" %(placement_problem.label,sum([len(a_tree.leaf_node_names()) for a_tree in alignment_tree_map.values()])))
            for (a_key, a_tree) in alignment_tree_map.items():
                assert isinstance(a_tree, PhylogeneticTree)  
                self.modify_tree(a_tree)
                alignment_problem  = SeppProblem(a_tree.leaf_node_names(), 
                                                  placement_problem)
                alignment_problem.subtree = a_tree
                alignment_problem.label = "A_%s_%s" %(str(p_key),str(a_key))                                                       
        
        ''' Divide fragments into chunks, to help achieve better parallelism'''
        fragment_chunk_files = self.create_fragment_files()                
        for alignment_problem in self.root_problem.iter_leaves():       
            for afc in xrange(0,len(fragment_chunk_files)):
                frag_chunk_problem  = SeppProblem(alignment_problem.taxa, 
                                              alignment_problem)
                frag_chunk_problem.subtree = alignment_problem.subtree
                frag_chunk_problem.label = alignment_problem.label.replace("A_", "FC_") + "_" +str(afc)
                frag_chunk_problem.fragments = fragment_chunk_files[afc]
                    
        _LOG.info("Breaking into %d alignment subsets." %(len(list(self.root_problem.iter_leaves()))))    
        _LOG.info("Breaking each alignment subset into %d fragment chunks." %len(fragment_chunk_files))
        _LOG.info("Subproblem structure: %s" %str(self.root_problem))
        return self.root_problem
    
    def create_fragment_files(self):
        alg_subset_count = len(list(self.root_problem.iter_leaves()))
        frag_chunk_count = lcm(alg_subset_count,self.options.cpu)//alg_subset_count
        return self.read_and_divide_fragments(frag_chunk_count)
         
    def _get_new_Join_Align_Job(self):
        return JoinAlignJobs()
    
    def build_jobs(self):        
        assert isinstance(self.root_problem, SeppProblem)
        for placement_problem in self.root_problem.get_children():
            ''' Create pplacer jobs'''
            pj = PplacerJob()
            placement_problem.add_job(pj.job_type,pj)
            pj.partial_setup_for_subproblem(placement_problem, self.options.info_file)
            
            '''For each alignment subproblem, ...'''
            for alg_problem in placement_problem.children:
                assert isinstance(alg_problem, SeppProblem)                
                ''' create the build model job'''
                bj = HMMBuildJob()
                #pdb.set_trace()
                bj.setup_for_subproblem(alg_problem,symfrac=self.symfrac, molecule=self.molecule,**vars(self.options.hmmbuild))
                alg_problem.add_job(bj.job_type, bj)                
                ''' create the search jobs'''
                for fc_problem in alg_problem.get_children():
                    sj = HMMSearchJob()
                    sj.partial_setup_for_subproblem(fc_problem.fragments, fc_problem, self.elim, self.filters)
                    fc_problem.add_job(sj.job_type, sj)                                    
                    ''' create the align job'''
                    aj = HMMAlignJob()
                    fc_problem.add_job(aj.job_type, aj)
                    aj.partial_setup_for_subproblem(fc_problem, molecule=self.molecule)

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
            '''Join all align jobs of a placement subset (enqueues placement job)'''
            jaj = self._get_new_Join_Align_Job()
            jaj.setup_with_placement_problem(placement_problem)                        
        ''' Join all search jobs together (enqueues align jobs)'''
        jsj = JoinSearchJobs()
        jsj.setup_with_root_problem(self.root_problem)        
        
    def enqueue_firstlevel_job(self):
        for p in self.root_problem.children:
            for ap in p.children:
                JobPool().enqueue_job(ap.jobs["hmmbuild"])

if __name__ == '__main__':
    ExhaustiveAlgorithm().run()
