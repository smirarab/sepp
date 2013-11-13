from sepp.exhaustive import ExhaustiveAlgorithm
import sepp
from sepp.config import options, default_settings
import argparse
from sepp.algorithm import AbstractAlgorithm
from sepp.jobs import MergeJsonJob
from sepp.tree import PhylogeneticTree
from dendropy.dataobject.tree import Tree
import dendropy
from sepp import get_logger

_LOG = get_logger(__name__)

class TIPPExhaustiveAlgorithm(ExhaustiveAlgorithm):
    '''
    This implements the exhaustive algorithm where all alignments subsets
    are searched for every fragment. This is for TIPP, meaning that we also
    perform classification based on a given taxonomy.
    '''
    def __init__(self):
        ExhaustiveAlgorithm.__init__(self)    
        self.alignment_threshold = self.options.alignment_threshold 
        self.push_down = True if self.options.tipp.push_down.lower() == "true" else False
        _LOG.info("Will push fragments %s from their placement edge." %("down" if self.push_down else "up"))

    def check_options(self, supply=[]):
        if (options().taxonomy_file is None):
            supply = supply + ["taxonomy file"]
        if (options().taxonomy_name_mapping_file is None):
            supply = supply + ["taxonomy name mapping file"]
        ExhaustiveAlgorithm.check_options(self, supply)
        
    def read_alignment_and_tree(self):
        (alignment, tree) = AbstractAlgorithm.read_alignment_and_tree(self)
        # TODO: Check for rooted input
        #if not tree.is_rooted:
        #    raise Exception ("For TIPP, backbone tree should be correctly rooted according to the taxonomy.")
        return (alignment, tree)
                
    def get_merge_job(self,meregeinputstring):
        mergeJsonJob = MergeJsonJob()
        mergeJsonJob.setup_for_tipp(meregeinputstring, 
                           self.get_output_filename("placement.json"), 
                           self.options.taxonomy_file, 
                           self.options.taxonomy_name_mapping_file, 
                           self.options.placement_threshold,
                           self.get_output_filename("classification.txt"),
                           self.push_down)
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
    default_settings['DEF_P'] = (100 , "Number of taxa (i.e. no decomposition)")
    parser = sepp.config.get_parser()
    uppGroup = parser.add_argument_group("TIPP Options".upper(), 
                         "These arguments set settings specific to TIPP")                                 
    
    uppGroup.add_argument("-at", "--alignmentThreshold", type = float, 
                      dest = "alignment_threshold", metavar = "N", 
                      default = 0.95,
                      help = "Enough alignment subsets are selected to reach a commulative probability of N. "
                             "This should be a number between 0 and 1 [default: 0.95]")                            

    uppGroup.add_argument("-pt", "--plaecementThreshold", type = float, 
                      dest = "placement_threshold", metavar = "N", 
                      default = 0.95,
                      help = "Enough placements are selected to reach a commulative probability of N. "
                             "This should be a number between 0 and 1 [default: 0.95]")                            
    
    uppGroup.add_argument("-tx", "--taxonomy", type = argparse.FileType('r'), 
                      dest = "taxonomy_file", metavar = "TAXONOMY", 
                      help = "A file describing the taxonomy. This is a comma-separated text file that has the following fields: "
                             "taxon_id,parent_id,taxon_name,rank. " 
                             "If there are other columns, they are ignored. The first line is also ignored.")
    
    uppGroup.add_argument("-txm", "--taxonomyNameMapping", type = argparse.FileType('r'), 
                      dest = "taxonomy_name_mapping_file", metavar = "MAPPING", 
                      help = "A comma-separated text file mapping alignment sequence names to taxonomic ids. "
                      "Formats (each line): "
                             "sequence_name,taxon_id. " 
                             "If there are other columns, they are ignored. The first line is also ignored.")

    uppGroup.add_argument("-adt", "--alignmentDecompositionTree", type = argparse.FileType('r'), 
                      dest = "alignment_decomposition_tree", metavar = "TREE", default = None,
                      help = "A newick tree file used for decomposing taxa into alignment subsets. " 
                             "[default: the backbone tree]")
        
def main():
    augment_parser() 
    TIPPExhaustiveAlgorithm().run()

if __name__ == '__main__':   
    main()