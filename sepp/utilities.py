'''
Created on Jun 7, 2011

@author: smirarab
'''
from dendropy import treecalc
from sepp import sortByValue
from Bio import AlignIO
from sepp.alignment import Alignment

class MyTreeCalc():
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''

    def branchOut(self,tree,centerTaxon,subsetSize,**kwargs):
        dist = {}                
        pdm = treecalc.PatristicDistanceMatrix(tree)
        for i , s in enumerate(tree.taxon_set): #@UnusedVariable
            if kwargs.has_key("filterTaxon"):
                if not kwargs["filterTaxon"](s):
                    continue;
            dist [s.label] = pdm(centerTaxon, s);
        incircle = sortByValue(dist)[0:subsetSize]
        return [node[0] for node in incircle]       

def convert_fasta_to_stockholm(in_path, out_path):
        in_fmt = 'fasta'
        out_fmt = 'stockholm'
        # open files
        in_hndl = open (in_path, 'rb')
        out_hndl = open (out_path, 'wb')
        # read and write
        in_aligns = [x for x in AlignIO.parse (in_hndl, in_fmt)]
        assert (in_aligns), '''Alignments at %s not valid''' % (in_path, in_fmt)
        AlignIO.write (in_aligns, out_hndl, out_fmt)

def convert_stockholm_to_fasta(in_path, out_path):
        out_fmt = 'fasta'
        in_fmt = 'stockholm'
        # open files
        in_hndl = open (in_path, 'rb')
        out_hndl = open (out_path, 'wb')
        # read and write
        in_aligns = [x for x in AlignIO.parse (in_hndl, in_fmt)]
        assert (in_aligns), '''Alignments at %s not valid''' % (in_path, in_fmt)
        AlignIO.write (in_aligns, out_hndl, out_fmt)

def read_sto_alignment(fn):
    try:
        ioA = AlignIO.read(fn, 'stockholm')
        alignment = Alignment()
        for row in ioA:
            alignment[row.id]=row.seq.data.upper().replace(".","-")
    except Exception as inst:
        print type(inst)   
        print inst.args    
        print inst
        raise RuntimeError("Could not read STO alignment")
        
    return alignment        

mytreeclac = MyTreeCalc()    