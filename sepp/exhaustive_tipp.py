from sepp.exhaustive import ExhaustiveAlgorithm
import sepp,os,stat,math
from sepp.config import options
import argparse,json
from sepp.algorithm import AbstractAlgorithm
from sepp.alignment import MutableAlignment, ExtendedAlignment,_write_fasta
from sepp.jobs import HMMBuildJob, HMMSearchJob, HMMAlignJob, PplacerJob,\
    MergeJsonJob,ExternalSeppJob
from sepp.problem import SeppProblem
from sepp.scheduler import JobPool, Join
from sepp.tree import PhylogeneticTree
from dendropy.dataobject.tree import Tree
import dendropy,pickle,pdb
from sepp import get_logger

_LOG = get_logger(__name__)

class TIPPJoinSearchJobs(Join):
    '''
    After all search jobs have finished on tips, we need to figure out which
    fragment goes to which subset and start aligning fragments.
    This join takes care of that step.
    '''
    def __init__(self, alignment_threshold):
        Join.__init__(self)
        self.alignment_threshold = alignment_threshold
        
    def setup_with_root_problem(self, root_problem):
        self.root_problem = root_problem
        for p in root_problem.iter_leaves():
            self.add_job(p.jobs["hmmsearch"])
        
    def figureout_fragment_subset(self):
        ''' Figure out which fragment should go to which subproblem'''
        # We need to keep and check the following flag because of checkpoining scenarios (join already done before!)
        if self.root_problem.annotations.has_key("fragments.distribution.done"):
            return
        bitscores = dict([(name, []) for name in self.root_problem.fragments.keys()])
        for fragment_chunk_problem in self.root_problem.iter_leaves():
            align_problem = fragment_chunk_problem.get_parent()
            assert isinstance(align_problem, SeppProblem)
            '''For each subproblem start with an empty set of fragments,
            and add to them as we encounter new best hits for that subproblem'''
            if align_problem.fragments is None:
                align_problem.fragments = MutableAlignment()
            search_res = fragment_chunk_problem.get_job_result_by_name("hmmsearch")
            for key in search_res.keys():
                ''' keep a list of all hits, and their bit scores'''
                bitscores[key].append( (search_res[key][1], align_problem) )
                
        for frag, tuplelist in bitscores.iteritems():
            ''' TODO: what to do with those that are not? For now, only output warning message'''
            #TODO:  Need to double check and fix the math
            if len(tuplelist) == 0:
                _LOG.warning("Fragment %s is not scored against any subset" %str(frag))
                continue
            ''' convert bit scores to probabilities '''            
            denum = sum(math.pow(2, min(x[0],1022)) for x in tuplelist)
            tuplelist = [((math.pow(2,min(x[0],1022))/denum*1000000),x[1]) for x in tuplelist]
            ''' Sort subsets by their probability'''
            tuplelist.sort(reverse=True)
            ''' Find enough subsets to reach the threshold '''
            selected = tuplelist[ 0 : max(1,
                reduce(lambda x, y: (x[0],None) if x[1] is None else
                                    (y[0],x[1]+y[1]) if x[1] < int(1000000 * self.alignment_threshold) else
                                    (y[0],None),
                       enumerate([x[0] for x in tuplelist]))[0]) ]
            
            ''' Renormalize the selected list to add up to 1'''
            renorm = 0
            for (prob,align_problem) in selected:	      
              renorm = renorm + prob/1000000
            renorm = 1/renorm
            
            _LOG.debug("Fragment %s assigned to %d subsets" %(frag,len(selected)))
            ''' Rename the fragment and assign it to the respective subsets'''
            for (prob,align_problem) in selected:
                postfix = prob*renorm if options().exhaustive.weight_placement_by_alignment.lower() == "true" else 1000000
                frag_rename = "%s_%s_%d" %(frag,align_problem.label,postfix)
                align_problem.fragments[frag_rename] = self.root_problem.fragments[frag]
        
        self.root_problem.annotations["fragments.distribution.done"] = 1

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

                if fragment_chunk_problem.fragments is not None and not fragment_chunk_problem.fragments.is_empty():
                    fragment_chunk_problem.fragments.write_to_path(aj.fragments)
                else:
                    aj.fake_run = True
                ''' Now the align job can be put on the queue '''
                JobPool().enqueue_job(aj)
                
    def __str__(self):
        return "join search jobs for all tips of ", self.root_problem

class TIPPJoinAlignJobs(Join):
    '''
    After all alignments jobs for a placement subset have finished,
    we need to build those extended alignments and start placing fragments.
    This join takes care of that step.
    '''
    def __init__(self, placer):
        Join.__init__(self)
        self.placer = placer
        
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
        for ap in pp.get_children():
            pp.fragments.seq_names |= set(ap.fragments)
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
        pj = pp.jobs["placer"]

        #Split the backbone alignment and query sequences into separate files        
        queryExtendedAlignment = fullExtendedAlignment.get_fragments_readonly_alignment()
        baseAlignment = fullExtendedAlignment.get_base_readonly_alignment()
    
        # Check for empty fragment files
        if (queryExtendedAlignment.is_empty()):
            pj.fake_run = True
                    
        elif self.placer == "pplacer":
            assert isinstance(pj,PplacerJob)
            
            #Write out the extended alignments, split into query and full-length for pplacer
            queryExtendedAlignment.write_to_path(pj.extended_alignment_file)
            baseAlignment.write_to_path(pj.backbone_alignment_file)
            
        elif self.placer == "epa":
            assert isinstance(pj,EPAJob)
                        
            #Write out the extended alignments in phylip for EPA
            fullExtendedAlignment.write_to_path(pj.extended_alignment_file, schema="PHYLIP")
                        
        #keep the extended alignment on everything
        #pj.set_attribute("full_extended_alignment_object", fullExtendedAlignment)
        
        #TODO:  Removed this, as it can cause unexpected lockups
        output = open(pj.full_extended_alignment_file, 'wb')
        pickle.dump(fullExtendedAlignment, output)
        output.close()                    
        
        # Enqueue the placement job        
        JobPool().enqueue_job(pj)
        
    def __str__(self):
        return "join align jobs for tips of ", self.placement_problem

class TIPPMergeJsonJob(ExternalSeppJob):
    def __init__(self, **kwargs):
        self.job_type = 'jsonmerger'
        ExternalSeppJob.__init__(self, self.job_type, **kwargs)
        self.out_file = None
        self.distribution = False
        self.taxonomy = None
        self.mapping = None
        self.threshold = None
        self.classification_file = None        
        self.elim = float(options().hmmsearch.elim)
        self.filters = True if options().hmmsearch.filters.upper() == "TRUE" else False if options().hmmsearch.filters.upper() == "FALSE" else None
        if self.filters is None:
            raise Exception("Expecting true/false for options().hmmsearch.filters")
        self.strategy = options().exhaustive.strategy
        self.minsubsetsize = int(options().exhaustive.minsubsetsize)
        self.alignment_threshold = float(options().alignment_threshold)
        #Temp fix for now,
        self.molecule = options().molecule
        self.placer = options().exhaustive.__dict__['placer'].lower()
        self.cutoff = 0;
        
    def setup(self, inString, output_file, **kwargs):
        self.stdindata = inString
        self.out_file = output_file
        self._kwargs = kwargs

    def setup_for_tipp(self, inString, output_file, taxonomy, mapping, threshold, classification_file, push_down, distribution=False,cutoff=0, **kwargs):
        self.stdindata = inString
        self.out_file = output_file
        self.taxonomy = taxonomy.name
        self.mapping = mapping.name
        self.distribution = distribution
        self.threshold = str(threshold)
        self.classification_file = classification_file
        self.push_down = push_down
        self._kwargs = kwargs
        self.cutoff = cutoff;
        
    def get_invocation(self):
        invoc = ["java", "-jar", self.path, "-", "-", self.out_file, "-r", "4"]        
        if self.taxonomy is not None:
            invoc.extend(["-t", self.taxonomy])
        if self.mapping is not None:
            invoc.extend(["-m", self.mapping])
        if self.threshold is not None:
            invoc.extend(["-p", self.threshold])
        if self.classification_file is not None:
            invoc.extend(["-c", self.classification_file])
        if self.distribution:
            invoc.extend(["-d"])            
        if not self.push_down:
            invoc.extend(["-u"])
        invoc.extend(["-C", str(self.cutoff)])
        #print " ".join(invoc)
        return invoc

    def characterize_input(self):
        return "input:pipe output:%s; Pipe:\n%s" %(self.out_file,self.stdindata)
 
    def read_results(self):
        '''
        Since the output file can be huge, we don't want to read it here, because
        it will need to get pickled and unpickled. Instead, we just send
        back the file name, and will let the caller figure out what to do with it.
        '''
        assert os.path.exists(self.out_file)
        assert os.stat(self.out_file)[stat.ST_SIZE] != 0
        return (self.out_file, self.stdoutdata)
        
class TIPPExhaustiveAlgorithm(ExhaustiveAlgorithm):
    '''
    This implements the exhaustive algorithm where all alignments subsets
    are searched for every fragment. This is for TIPP, meaning that we also
    perform classification based on a given taxonomy.
    '''
    def __init__(self):
        ExhaustiveAlgorithm.__init__(self)    
        self.alignment_threshold = self.options.alignment_threshold 
        self.placer = self.options.exhaustive.placer.lower()
        self.push_down = True if self.options.push_down == True else False
        _LOG.info("Will push fragments %s from their placement edge." %("down" if self.push_down else "up"))                               
    def _get_new_Join_Align_Job(self):
        return TIPPJoinAlignJobs(self.placer)

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
        jsj = TIPPJoinSearchJobs(self.alignment_threshold)
        jsj.setup_with_root_problem(self.root_problem)
        
    def merge_results(self):
        ''' TODO: implement this
        '''
        assert isinstance(self.root_problem,SeppProblem)
        
        '''Generate single extended alignment'''
        align_input = open(self.root_problem.get_children()[0].jobs["placer"].full_extended_alignment_file,'r')
        fullExtendedAlignment = pickle.load(align_input)
        align_input.close()
        #fullExtendedAlignment = self.root_problem.get_children()[0].jobs["placer"].get_attribute("full_extended_alignment_file")
        for pp in self.root_problem.get_children()[1:]:
            #Removed this because it can cause unexpected lockups
            #extended_alignment = pp.jobs["placer"].get_attribute("full_extended_alignment_object")
            align_input = open(pp.jobs["placer"].full_extended_alignment_file,'r')
            extended_alignment = pickle.load(align_input)
            align_input.close()
            
            #fullExtendedAlignment.merge_in(extended_alignment,convert_to_string=True)
            fullExtendedAlignment.merge_in(extended_alignment,convert_to_string=True)
        self.results = fullExtendedAlignment
        
        mergeinput = []
        '''Append main tree to merge input'''
        mergeinput.append("%s;" %(self.root_problem.subtree.compose_newick(labels = True)))
        jsons = []
        for pp in self.root_problem.get_children():
            assert isinstance(pp,SeppProblem)
            if (pp.get_job_result_by_name("placer") is None):
              continue
            '''Append subset trees and json locations to merge input'''
            mergeinput.append("%s;\n%s" %(pp.subtree.compose_newick(labels = True),
                              pp.get_job_result_by_name("placer")))
        mergeinput.append("")
        mergeinput.append("")
        meregeinputstring = "\n".join(mergeinput)
        mergeJsonJob = self.get_merge_job(meregeinputstring)
        mergeJsonJob.run()        
        
    def build_jobs(self):
        assert isinstance(self.root_problem, SeppProblem)
        for placement_problem in self.root_problem.get_children():
            ''' Create placer jobs'''
            if self.placer == "pplacer":
                pj = PplacerJob()
                pj.partial_setup_for_subproblem(placement_problem, self.options.info_file)
            elif self.placer == "epa":
                pj = EPAJob()
                pj.partial_setup_for_subproblem(placement_problem, self.molecule)
                
            placement_problem.add_job("placer",pj)
            
            '''For each alignment subproblem, ...'''
            for alg_problem in placement_problem.children:
                assert isinstance(alg_problem, SeppProblem)
                ''' create the build model job'''
                bj = HMMBuildJob()
                bj.setup_for_subproblem(alg_problem,molecule=self.molecule)
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

        

    def check_options(self, supply=[]):
        if (options().reference_pkg is not None):
            self.load_reference(os.path.join(options().reference.path, 'refpkg/%s.refpkg/' % options().reference_pkg))                  
        if (options().taxonomy_file is None):
            supply = supply + ["taxonomy file"]
        if (options().taxonomy_name_mapping_file is None):
            supply = supply + ["taxonomy name mapping file"]
        ExhaustiveAlgorithm.check_options(self, supply)
        
    def load_reference(self, reference_pkg):
        file = open(reference_pkg + 'CONTENTS.json')
        result=json.load(file)
        file.close()
        options().taxonomy_name_mapping_file = open(reference_pkg + result['files']['seq_info'])
        options().taxonomy_file = open(reference_pkg + result['files']['taxonomy'])
        options().alignment_file = open(reference_pkg + result['files']['aln_fasta'])
        options().tree_file = open(reference_pkg + result['files']['tree'])
        options().info_file = reference_pkg + result['files']['tree_stats']
        
        
    def read_alignment_and_tree(self):
        (alignment, tree) = AbstractAlgorithm.read_alignment_and_tree(self)
        # TODO: Check for rooted input
        #if not tree.is_rooted:
        #    raise Exception ("For TIPP, backbone tree should be correctly rooted according to the taxonomy.")
        return (alignment, tree)
                
    def get_merge_job(self,meregeinputstring):
        mergeJsonJob = TIPPMergeJsonJob()
        mergeJsonJob.setup_for_tipp(meregeinputstring, 
                           self.get_output_filename("placement.json"), 
                           self.options.taxonomy_file, 
                           self.options.taxonomy_name_mapping_file, 
                           self.options.placement_threshold,
                           self.get_output_filename("classification.txt"),
                           self.push_down,
                           self.options.distribution,
                           self.options.cutoff)
        return mergeJsonJob

    def get_alignment_decomposition_tree(self, p_tree):
        assert isinstance(p_tree, PhylogeneticTree)
        if self.options.alignment_decomposition_tree is None:
            return PhylogeneticTree(Tree(p_tree.den_tree))
        elif p_tree.count_leaves() != self.root_problem.subtree.count_leaves():
            raise ValueError("Alignment decomposition tree can be different from placement tree only if placement subset size is set to the number of taxa (i.e. entire tree)")
        else:
            _LOG.info("Reading alignment decomposition input tree: %s" %(self.options.alignment_decomposition_tree))        
            d_tree = PhylogeneticTree( dendropy.Tree(stream=self.options.alignment_decomposition_tree, 
                                               schema="newick", 
                                               preserve_underscores=True,
                                               taxon_set=self.root_problem.subtree.get_tree().taxon_set))               
            return d_tree
            
def augment_parser():
    sepp.config.set_main_config_path(os.path.expanduser("~/.sepp/tipp.config"))
    #default_settings['DEF_P'] = (100 , "Number of taxa (i.e. no decomposition)")
    parser = sepp.config.get_parser()
    tippGroup = parser.add_argument_group("TIPP Options".upper(), 
                         "These arguments set settings specific to TIPP")

    tippGroup.add_argument("-R", "--reference_pkg", type = str, 
                      dest = "reference_pkg", metavar = "N", 
                      default = None,
                      help = "Use a pre-computed reference package [default: None]")                         
    
    tippGroup.add_argument("-at", "--alignmentThreshold", type = float, 
                      dest = "alignment_threshold", metavar = "N", 
                      default = 0.95,
                      help = "Enough alignment subsets are selected to reach a commulative probability of N. "
                             "This should be a number between 0 and 1 [default: 0.95]")                            
    tippGroup.add_argument("-D", "--dist", 
                      dest = "distribution", action='store_true', 
                      default = False,
                      help = "Treat fragments as distribution")    
                             
    tippGroup.add_argument("-pt", "--placementThreshold", type = float, 
                      dest = "placement_threshold", metavar = "N", 
                      default = 0.95,
                      help = "Enough placements are selected to reach a commulative probability of N. "
                             "This should be a number between 0 and 1 [default: 0.95]")                            
    
    tippGroup.add_argument("-PD", "--push_down", type = bool,
                      dest = "push_down", metavar = "N",
                      default = True,
                      help = "Whether to classify based on children below or above insertion point.  [default: True]")


    tippGroup.add_argument("-tx", "--taxonomy", type = argparse.FileType('r'), 
                      dest = "taxonomy_file", metavar = "TAXONOMY", 
                      help = "A file describing the taxonomy. This is a comma-separated text file that has the following fields: "
                             "taxon_id,parent_id,taxon_name,rank. " 
                             "If there are other columns, they are ignored. The first line is also ignored.")
    
    tippGroup.add_argument("-txm", "--taxonomyNameMapping", type = argparse.FileType('r'), 
                      dest = "taxonomy_name_mapping_file", metavar = "MAPPING", 
                      help = "A comma-separated text file mapping alignment sequence names to taxonomic ids. "
                      "Formats (each line): "
                             "sequence_name,taxon_id. " 
                             "If there are other columns, they are ignored. The first line is also ignored.")

    tippGroup.add_argument("-adt", "--alignmentDecompositionTree", type = argparse.FileType('r'), 
                      dest = "alignment_decomposition_tree", metavar = "TREE", default = None,
                      help = "A newick tree file used for decomposing taxa into alignment subsets. " 
                             "[default: the backbone tree]")
                             
    tippGroup.add_argument("-C", "--cutoff", type = float, 
                      dest = "cutoff", metavar = "N", 
                      default = 0.0,
                      help = "Placement probability requirement to count toward the distribution. "
                             "This should be a number between 0 and 1 [default: 0.0]")    
                             
        
def main():
    augment_parser() 
    TIPPExhaustiveAlgorithm().run()

if __name__ == '__main__':   
    main()
    
    
