'''
Created on Sep 19, 2012

@author: smirarab
'''
import unittest
from sepp.alignment import MutableAlignment, ReadonlySubalignment,\
    ExtendedAlignment
from sepp.problem import SeppProblem
from sepp.filemgr import get_data_path
from tempfile import mkstemp
from os import remove


class Test(unittest.TestCase):
    fp_dummy1 = None
    fp_dummy2 = None
    fp_dummy3 = None
    fp_dummy4 = None
    fp_dummy5 = None

    def setUp(self):
        _, self.fp_dummy1 = mkstemp()
        _, self.fp_dummy2 = mkstemp()
        _, self.fp_dummy3 = mkstemp()
        _, self.fp_dummy4 = mkstemp()
        _, self.fp_dummy5 = mkstemp()

    def tearDown(self):
        remove(self.fp_dummy1)
        remove(self.fp_dummy2)
        remove(self.fp_dummy3)
        remove(self.fp_dummy4)
        remove(self.fp_dummy5)

    def testAlignmentReadFasta(self):
        alg = MutableAlignment()
        alg.read_filepath(get_data_path("mock/pyrg/sate.fasta"))

        assert len(alg) == 65, "MutableAlignment length is %s" % len(alg)

        assert all([not alg.is_all_gap(i) for i in range(0, alg.get_length())])

    def testReadOnlySubAlignment(self):
        alg = MutableAlignment()
        alg.read_filepath(get_data_path("mock/pyrg/sate.fasta"))

        subset = ['NC_008701_720717_722309', 'NC_013156_149033_150643',
                  'NC_013887_802739_801129']
        readonly_subalignment = ReadonlySubalignment(subset, alg)

        assert len(readonly_subalignment) == 3, len(readonly_subalignment)

        assert set(readonly_subalignment.keys()) == set(
            readonly_subalignment.get_sequence_names()) == set(subset), \
            "Subalignment keys not matching given keys %s vs %s" % (
            list(readonly_subalignment.keys()), subset)

        for (k, s) in list(readonly_subalignment.items()):
            assert k in subset, \
                "%s not found in subset but returned by subalignment" % k
            assert s == alg[k], \
                "sequence associated with %s not matching parent alignment" % k

        try:
            readonly_subalignment[2] = "ACGT"
            assert False, "Readony alignment is successfully modified. "
        except TypeError:
            pass

        assert readonly_subalignment.get_length() == alg.get_length(), \
            "alignment length should not change"

        assert readonly_subalignment.is_aligned() is True

        assert readonly_subalignment.is_all_gap(2) is True, \
            "Site 2 should be all gaps"
        assert readonly_subalignment.is_all_gap(150) is False, \
            "Site 100 should not be all gaps"

        readonly_subalignment.write_to_path(
            self.fp_dummy1)  # "mock/pyrg/sate.sub.fasta"

        mutable_subalignment = readonly_subalignment.get_mutable_alignment()
        mutable_subalignment.delete_all_gap()

        assert all([not mutable_subalignment.is_all_gap(i)
                    for i in range(0, mutable_subalignment.get_length())])

    def testExtendedAlignment(self):
        subset = [
            "SFIF", "SFII", "SCFC", "SGHD", "SDCC", "SBGE", "SFBB", "SDI",
            "SCGB", "SJGF", "SGBI", "SCJA", "SGAD", "SHEB", "SFHB", "SDJI",
            "SHED", "SJJJ", "SBBE", "SCCH", "SDJB", "SDAC", "SHEH", "SFDC",
            "SFEI", "SHHB", "SC", "SIAB", "SDDI", "SBCB", "SJB", "SEBD",
            "SFGD", "SHA", "SIDA", "SGHI", "SGIB", "SBFJ", "SFIE", "SCJF",
            "SJHJ", "SJBG", "SEJI", "SFFF", "SJ", "SIII", "SJHH", "SEIH",
            "SBDC", "SHDJ", "SJDD", "SGDB", "SIHA", "SIBB", "SECC", "SCAD",
            "SGBB", "SGIF", "SJHC", "SFCD", "SEAA", "SEFF", "SDFG", "SDJE",
            "SCFG", "SFH", "SCJ", "SDDD", "SEGD", "SCIH", "SDAG", "SCJE",
            "SFAJ", "SIDJ", "SE", "SHBC", "SJFF", "SCHD", "SBHA", "SEDF",
            "SFAF", "SEDD", "SDHD", "SGJD", "SIBH", "SGDF", "SIFA", "SJGA",
            "SIJB", "SFI", "SGA", "SBFC", "SBJA", "SFFC", "SFDH", "SFEE",
            "SBDF", "SGBJ", "SDHE", "SJIB", "SHHI", "SIDE", "SJII"]

        alg = MutableAlignment()
        alg.read_filepath(get_data_path("simulated/test.fasta"))
        alg.delete_all_gap()
        tlen = alg.get_length()

        frg = MutableAlignment()
        frg.read_filepath(get_data_path("simulated/test.fas"))
        # print frg.get_num_taxa()

        pp = SeppProblem(list(alg.keys()))
        pp.fragments = frg
        pp.subalignment = alg

        cp1 = SeppProblem(subset, pp)
        cp2 = SeppProblem(list(set(alg.keys()) - set(subset)), pp)
        cp1.fragments = ReadonlySubalignment(
            [k for k in list(frg.keys()) if int(k[-1]) >= 9], frg)
        cp2.fragments = ReadonlySubalignment(
            [k for k in list(frg.keys()) if int(k[-1]) <= 1], frg)

        cp1labels = cp1.write_subalignment_without_allgap_columns(
            self.fp_dummy1)  # tmp/cp1.fasta
        cp2labels = cp2.write_subalignment_without_allgap_columns(
            self.fp_dummy2)  # tmp/cp2.fasta
        tmp = MutableAlignment().read_filepath(self.fp_dummy1)
        assert all([not tmp.is_all_gap(pos)
                    for pos in range(0, tmp.get_length())])
        tmp = MutableAlignment().read_filepath(self.fp_dummy2)
        assert all([not tmp.is_all_gap(pos)
                    for pos in range(0, tmp.get_length())])

        cp1.fragments.write_to_path(self.fp_dummy3)  # tmp/cp1.frags.fas
        cp2.fragments.write_to_path(self.fp_dummy4)  # tmp/cp2.frags.fas

        '''We have done the hmmalign before.
           don't worry about that right now'''

        ext1 = ExtendedAlignment(cp1.fragments)
        ext1.build_extended_alignment(
            self.fp_dummy1,
            get_data_path("tmp/cp1.extended.sto"))
        ext1.relabel_original_columns(cp1labels)
        ext2 = ExtendedAlignment(cp2.fragments)
        ext2.build_extended_alignment(
            self.fp_dummy2,
            get_data_path("tmp/cp2.extended.sto"))
        ext2.relabel_original_columns(cp2labels)

        extmerger = ExtendedAlignment([])
        extmerger.merge_in(ext1)
        mixed = extmerger.merge_in(ext2)

        extmerger.write_to_path(self.fp_dummy5)  # tmp/extended.merged.fasta

        assert extmerger.is_aligned(), "Merged alignment is not aligned"
        in1 = len([x for x in ext1._col_labels if x < 0])
        in2 = len([x for x in ext2._col_labels if x < 0])
        assert (in1 + in2 + tlen - mixed) == extmerger.get_length(), \
            ("Lengths don't match up after merging. Merged:%d. Insertion1:%d "
             "Insertion2:%d BaseLen:%d Mixed-insertion: %d") % (
                extmerger.get_length(), in1, in2, tlen, mixed)
        assert (in1 + in2 - mixed) == len(extmerger.get_insertion_columns()), \
            ("Columns are not correctly labeled after merging. Merged "
             "insertion count:%d. Insertion1:%d Insertion2:%d Mixed-insertion:"
             " %d") % (
             len(list(extmerger.iter_insertion_columns())), in1, in1, mixed)

        tmp = extmerger.get_base_readonly_alignment().get_mutable_alignment()
        tmp.delete_all_gap()
        assert tmp.is_aligned(), "merged alignment should be aligned!"
        assert tmp.get_length() == tlen, "merged alignment has wrong length"
        assert all([alg[k] == s for (k, s) in list(tmp.items())]), \
            "merged alignment should match original alignment"


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testAlignmentReadFasta']
    unittest.main()
