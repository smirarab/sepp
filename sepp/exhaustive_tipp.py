from sepp.exhaustive import ExhaustiveAlgorithm
import sepp
from sepp.config import options
import argparse
from sepp.algorithm import AbstractAlgorithm
from sepp.jobs import MergeJsonJob


class TIPPExhaustiveAlgorithm(ExhaustiveAlgorithm):
    '''
    This implements the exhaustive algorithm where all alignments subsets
    are searched for every fragment. This is for TIPP, meaning that we also
    perform classification based on a given taxonomy.
    '''
    def __init__(self):
        ExhaustiveAlgorithm.__init__(self)    
        self.alignment_threshold = self.options.alignment_threshold 

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
                           self.get_output_filename("classification.txt"))
        return mergeJsonJob
        
def augment_parser():
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
                      help = "A file containing the mapping between the alignment entries and taxonomic ids. "
                      "This is a comma-separated text file that has the following fields: "
                             "sequence_name,taxon_id. " 
                             "If there are other columns, they are ignored. The first line is also ignored.")
        


if __name__ == '__main__':   
    augment_parser() 
    TIPPExhaustiveAlgorithm().run()