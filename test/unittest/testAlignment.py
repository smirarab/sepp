'''
Created on Sep 19, 2012

@author: smirarab
'''
import unittest
import sepp
sepp._DEBUG = True
from sepp.alignment import MutableAlignment, ReadonlySubalignment,\
    ExtendedAlignment
from sepp.problem import SeppProblem        



class Test(unittest.TestCase):

    def setUp(self):
        pass
    
    def tearDown(self):
        pass


    def testAlignmentReadFasta(self):
        print "====== starting testAlignmentReadFasta ==========" 
        alg = MutableAlignment()
        alg.read_filepath("data/mock/pyrg/sate.fasta")
        
        print "Maing alignment is:\n\n", alg
        
        assert len(alg) == 65, "MutableAlignment length is %s" %len(alg)
        
        assert all([not alg.is_all_gap(i) for i in xrange(0,alg.get_length())])            
            

    def testReadOnlySubAlignment(self):
        print "======= starting testReadOnlySubAlignment =========" 
        alg = MutableAlignment()
        alg.read_filepath("data/mock/pyrg/sate.fasta")
        
        subset  = alg.keys()[9:12]
        readonly_subalignment = ReadonlySubalignment(subset, alg)
        
        print "subalignment is:\n\n", readonly_subalignment
        
        assert len(readonly_subalignment) == 3, len(readonly_subalignment) 
        
        assert readonly_subalignment.keys() == readonly_subalignment.get_sequence_names() == subset, "Subalignment keys not matching given keys %s vs %s" %(readonly_subalignment.keys() , subset)
        
        for (k, s) in readonly_subalignment.items():
            assert k in subset, "%s not found in subset but returned by subalignment" %k
            assert s == alg[k], "sequence associated with %k not matching parent alignment" %k 
        
        try:
            readonly_subalignment[2] = "ACGT"
            assert False, "Readony alignment is successfully modified. "
        except TypeError:
            pass
        
        assert readonly_subalignment.get_length() == alg.get_length(), "alignment length should not change"
        
        assert readonly_subalignment.is_aligned() == True
        
        assert readonly_subalignment.is_all_gap(2) == True, "Site 2 should be all gaps"
        assert readonly_subalignment.is_all_gap(150) == False, "Site 100 should not be all gaps"        
        
        readonly_subalignment.write_to_path("data/mock/pyrg/sate.sub.fasta")
        
        mutable_subalignment = readonly_subalignment.get_mutable_alignment()
        mutable_subalignment.delete_all_gap()
        
        assert all([not mutable_subalignment.is_all_gap(i) for i in xrange(0,mutable_subalignment.get_length())])
        
        print "======= finishing testReadOnlySubAlignment =========" 

    def testExtendedAlignment(self):
        print "======= starting testExtendedAlignment ========="

        subset = ["SFIF","SFII","SCFC","SGHD","SDCC","SBGE","SFBB","SDI","SCGB","SJGF","SGBI","SCJA","SGAD","SHEB","SFHB","SDJI","SHED","SJJJ","SBBE","SCCH","SDJB","SDAC","SHEH","SFDC","SFEI","SHHB","SC","SIAB","SDDI","SBCB","SJB","SEBD","SFGD","SHA","SIDA","SGHI","SGIB","SBFJ","SFIE","SCJF","SJHJ","SJBG","SEJI","SFFF","SJ","SIII","SJHH","SEIH","SBDC","SHDJ","SJDD","SGDB","SIHA","SIBB","SECC","SCAD","SGBB","SGIF","SJHC","SFCD","SEAA","SEFF","SDFG","SDJE","SCFG","SFH","SCJ","SDDD","SEGD","SCIH","SDAG","SCJE","SFAJ","SIDJ","SE","SHBC","SJFF","SCHD","SBHA","SEDF","SFAF","SEDD","SDHD","SGJD","SIBH","SGDF","SIFA","SJGA","SIJB","SFI","SGA","SBFC","SBJA","SFFC","SFDH","SFEE","SBDF","SGBJ","SDHE","SJIB","SHHI","SIDE","SJII"]
         
        alg = MutableAlignment()
        alg.read_filepath("data/simulated/test.fasta")
        alg.delete_all_gap()
        tlen = alg.get_length()                    
        
        frg = MutableAlignment()
        frg.read_filepath("data/simulated/test.fas")
        #print frg.get_num_taxa()
        
        pp = SeppProblem(alg.keys())
        pp.fragments = frg
        pp.subalignment = alg
        
        cp1 = SeppProblem(subset, pp)
        cp2 = SeppProblem(list(set(alg.keys()) -set(subset)), pp)
        cp1.fragments = ReadonlySubalignment([k for k in frg.keys() if int(k[-1]) >= 9], frg)
        cp2.fragments = ReadonlySubalignment([k for k in frg.keys() if int(k[-1]) <= 1], frg)
        
        cp1labels = cp1.write_subalignment_without_allgap_columns("data/tmp/cp1.fasta")
        cp2labels = cp2.write_subalignment_without_allgap_columns("data/tmp/cp2.fasta")
        tmp = MutableAlignment().read_filepath("data/tmp/cp1.fasta")
        assert all([not tmp.is_all_gap(pos) for pos in xrange(0,tmp.get_length())])        
        tmp = MutableAlignment().read_filepath("data/tmp/cp2.fasta")
        assert all([not tmp.is_all_gap(pos) for pos in xrange(0,tmp.get_length())])
        
        cp1.fragments.write_to_path("data/tmp/cp1.frags.fas")
        cp2.fragments.write_to_path("data/tmp/cp2.frags.fas")
        
        '''We have done the hmmalign before. don't worry about that right now'''
        
        ext1 = ExtendedAlignment(cp1.fragments)
        ext1.build_extended_alignment("data/tmp/cp1.fasta", "data/tmp/cp1.extended.sto")
        ext1.relabel_original_columns(cp1labels)
        ext2 = ExtendedAlignment(cp2.fragments)
        ext2.build_extended_alignment("data/tmp/cp2.fasta", "data/tmp/cp2.extended.sto")
        ext2.relabel_original_columns(cp2labels)
        
        extmerger = ExtendedAlignment([])
        extmerger.merge_in(ext1)
        mixed = extmerger.merge_in(ext2)
                        
        extmerger.write_to_path("data/tmp/extended.merged.fasta")        

        assert extmerger.is_aligned(), "Merged alignment is not aligned"
        in1 = len([x for x in ext1._col_labels if x<0])
        in2 = len([x for x in ext2._col_labels if x<0])
        print "Merged:%d. Insertion1:%d Insertion2:%d BaseLen:%d" %(extmerger.get_length(),in1 , in2 , tlen)
        assert ( in1 + in2 + tlen - mixed) == extmerger.get_length(), "Lengths don't match up after merging. Merged:%d. Insertion1:%d Insertion2:%d BaseLen:%d Mixed-insertion: %d"  %(extmerger.get_length(),in1, in2 , tlen, mixed)
        assert ( in1 + in2 - mixed) == len(list(extmerger.iter_insertion_columns())), "Columns are not correctly labeled after merging. Merged insertion count:%d. Insertion1:%d Insertion2:%d Mixed-insertion: %d"  %(len(list(extmerger.iter_insertion_columns())),in1 , in1, mixed)
         
        
        tmp = extmerger.get_base_readonly_alignment().get_mutable_alignment()
        tmp.delete_all_gap()
        assert tmp.is_aligned(), "merged alignment should be aligned!"
        assert tmp.get_length() == tlen, "merged alignment has wrong length"
        assert all([alg[k] == s for (k,s) in tmp.items()]), "merged alignment should match original alignment"

        
        print "======= finished testExtendedAlignment ========="
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testAlignmentReadFasta']
    unittest.main()