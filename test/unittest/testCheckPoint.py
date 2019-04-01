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
from sepp.config import set_checkpoint
import sepp.scheduler

from sepp.exhaustive import ExhaustiveAlgorithm

# sepp._DEBUG = True
# sepp.reset_loggers()


class Test(unittest.TestCase):
    x = None

    def resetSepp(self):
        # ensure necessary settings are made
        sepp.scheduler._jobPool = None
        sepp.scheduler._parser = None
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
        self.x.options.fragment_file = open(
            get_data_path(
                "q2-fragment-insertion/input_fragments.fasta"), "r")
        self.x.options.outdir = tempfile.mkdtemp()
        self.x.options.placement_size = 20

    def setUp(self):
        sys.argv = [sys.argv[0]]
        self.resetSepp()
        self._checkpointfile = tempfile.mktemp()

    def tearDown(self):
        shutil.rmtree(self.x.options.outdir, ignore_errors=True)
        self.resetSepp()

    def test_make_checkpoints(self):
        self.x.options.checkpoint = set_checkpoint(self._checkpointfile)
        self.x.options.checkpoint_interval = 5
        self.x.run()
        self.assertTrue(self.x.results is not None)


if __name__ == "__main__":
    unittest.main()
