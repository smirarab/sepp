'''
Collection of functions for metagenomic pipeline for taxonomic classification
Created on June 3, 2014

@author: namphuon
'''

class JoinBlastJobs(Join):
    '''
    After all blast search jobs have finished on markers, we need to figure out which 
    fragment goes to which marker. 
    This join takes care of that step. 
    '''
    def __init__(self):
        Join.__init__(self)
        
    def setup_with_root_problem(self, root_problem):
        self.root_problem = root_problem            
        for p in root_problem.iter_leaves():
            self.add_job(p.jobs["blastsearch"])         
    
    def figureout_fragment_marker(self):
        ''' Figure out which fragment should go to which marker'''
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

                if not fragment_chunk_problem.fragments.is_empty():
                    fragment_chunk_problem.fragments.write_to_path(aj.fragments)
                else:
                    aj.fake_run = True
                ''' Now the align job can be put on the queue '''
                JobPool().enqueue_job(aj)                
                
    def __str__(self):
        return "join search jobs for all tips of ", self.root_problem


def blast_fragments(config, input, output):
  '''Blast the fragments against all marker genes+16S sequences
  '''
  
  
def fix_direction(config,input,output):
  '''Fixes the direction of all the reads by searching the
     sequences against each HMM and finding out which one works the best
  '''  
def reverse_sequence(sequence):
  '''Reverse a sequence to be in the same direction as marker sequences  
  '''
  
def read_blast_results(input):
  '''Reads the results for blast
  '''
  
    #def read_results(self):
        #'''
        #Read the Sate log file and get the alignment and tree from file, copy to output directory
        #'''
        #assert os.path.exists('%s/sateout/satejob.out.txt' % sepp.filemgr.get_root_temp_dir())
        #assert os.stat('%s/sateout/satejob.out.txt' % sepp.filemgr.get_root_temp_dir())[stat.ST_SIZE] != 0        
        #outfile = open('%s/sateout/satejob.out.txt' % sepp.filemgr.get_root_temp_dir(), 'r');
        #alignment_pattern = re.compile('Writing resulting alignment to (.*)')
        #tree_pattern = re.compile('Writing resulting tree to (.*)')
        #tree_file = ''
        #alignment_file = ''        
        #for line in outfile:            
            #line = line.strip()            
            #result = alignment_pattern.findall(line)
            #if (len(result) != 0):
                #alignment_file = result[0]
            #result = tree_pattern.findall(line)
            #if (len(result) != 0):
                #tree_file = result[0]
        #shutil.copyfile(tree_file, "%s/sate.fasttree" % self.output)
        #shutil.copyfile(alignment_file, "%s/sate.fasta" % self.output)
        #return (tree_file,alignment_file)
