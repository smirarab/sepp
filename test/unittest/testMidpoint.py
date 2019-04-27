'''
Created on June 14, 2018

@author: Stefan.M.Janssen@gmail.com
'''
import sys
import tempfile
import shutil
import unittest
import sepp
from sepp.filemgr import get_data_path
import platform
from argparse import Namespace
import sepp.scheduler

from sepp.exhaustive import ExhaustiveAlgorithm
# sepp._DEBUG = True
# sepp.reset_loggers()


class Test(unittest.TestCase):
    x = None

    def resetSepp(self):
        sepp.scheduler._jobPool = None
        self.x = ExhaustiveAlgorithm()
        self.x.options.alignment_file = open(
            get_data_path(
                "q2-fragment-insertion/reference_alignment_tiny.fasta"), "r")
        self.x.options.info_file = open(
            get_data_path(
                "q2-fragment-insertion/RAxML_info-reference-gg-raxml-bl.info"),
            "r")
        self.x.options.tree_file = open(
            get_data_path(
                "q2-fragment-insertion/reference_phylogeny_tiny.nwk"), "r")
        self.x.options.outdir = tempfile.mkdtemp()

        suff_bit = "-64" if sys.maxsize > 2**32 else "-32"
        if platform.system() == 'Darwin':
            suff_bit = ""
        for prog in ['hmmalign', 'hmmbuild', 'hmmsearch', 'pplacer']:
            setattr(self.x.options, prog, Namespace(
                path=get_data_path("../../../tools/bundled/%s/%s%s" % (
                    platform.system(), prog, suff_bit))))

    def setUp(self):
        # ensure necessary settings are made
        sys.argv = [sys.argv[0], "-c", get_data_path("configs/test2.config")]
        self.resetSepp()

    def tearDown(self):
        shutil.rmtree(self.x.options.outdir, ignore_errors=True)
        self.resetSepp()

    def test_diamMid(self):
        self.x.options.fragment_file = open(
            get_data_path(
                "q2-fragment-insertion/input_fragments.fasta"), "r")
        self.x.options.maxDiam = 0.1
        self.x.options.fragmentChunkSize = 1
        self.x.options.decomp_strategy = "midpoint"
        self.x.run()
        self.assertTrue(self.x.results is not None)

    def test_diamCent(self):
        self.x.options.fragment_file = open(
            get_data_path(
                "q2-fragment-insertion/input_fragments.fasta"), "r")
        self.x.options.maxDiam = 0.1
        self.x.options.fragmentChunkSize = 1
        self.x.options.decomp_strategy = "centroid"
        self.x.run()
        self.assertTrue(self.x.results is not None)


if __name__ == "__main__":
    unittest.main()
