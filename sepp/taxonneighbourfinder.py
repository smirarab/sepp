'''
Created on Jun 7, 2011

@author: smirarab
'''
from sepp.utilities import mytreeclac, is_taxon_a_full_length_seq

class TaxonNeighbourFinder():

    def __init__(self, taxonSet):
        '''
        Constructor
        '''
        self.taxonSet = taxonSet
        self.subsets = None
 
    def findTaxonSubsets(self, readsTaxonList):
        raise NotImplementedError('Abstract TaxonNeighbourFinder class cannot find subsets.')         

    def writeSubsetsToFile(self, fileName):
        f = open(fileName,'w')
        for s in self.subsets.keys():
            f.write("%s %s\n" %(s, self.subsets.get(s)));
        f.close();
        
class PPlacerResultsNeighbour(TaxonNeighbourFinder):
    
    def __init__(self,taxonSet, pplacerTogTree):
        TaxonNeighbourFinder.__init__(self, taxonSet)
        self.togTree = pplacerTogTree
        self.filterToCount = 100
        
    def findTaxonSubsets(self, readsTaxonList):
        for taxon in readsTaxonList:
            subset = mytreeclac.branchOut(self.togTree, taxon, self.filterToCount, is_taxon_a_full_length_seq)
            self.subsets[taxon] = subset                    