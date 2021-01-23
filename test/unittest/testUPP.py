"""
Created on Nov 4, 2020

@author: Siavash Mirarab
"""
import sys
import tempfile
import shutil
import unittest
import sepp
from sepp.exhaustive_upp import UPPExhaustiveAlgorithm
from sepp.filemgr import get_data_path
import platform
from argparse import Namespace

# sepp._DEBUG = True
# sepp.reset_loggers()


class Test(unittest.TestCase):
    x = None

    def setUp(self):
        # ensure necessary settings are made
        sepp.scheduler._jobPool = None
        sys.argv = [sys.argv[0], "-c", get_data_path("configs/test3.config")]
        self.x = UPPExhaustiveAlgorithm()
        self.x.options.backbone_size = 407
        self.x.options.alignment_size = 50
        self.x.options.placement_size = None
        self.x.options.backtranslation_sequence_file = None
        self.x.options.long_branch_filter = 100000
        self.x.options.molecule = 'amino'
        self.x.molecule = 'amino'
        self.x.options.alignment_file = open(
            get_data_path(
                "upp_frag/backbone_pasta.fasta"), "r")
        self.x.options.sequence_file = open(
            get_data_path(
                "upp_frag/query.fas"),
            "r")
        self.x.options.tree_file = open(
            get_data_path(
                "upp_frag/backbone_pasta.fasttree"), "r")
        self.x.options.outdir = tempfile.mkdtemp()

        suff_bit = "-64" if sys.maxsize > 2**32 else "-32"
        if platform.system() == 'Darwin':
            suff_bit = ""
        for prog in ['hmmalign', 'hmmbuild', 'hmmsearch']:
            setattr(self.x.options, prog, Namespace(
                path=get_data_path("../../../tools/bundled/%s/%s%s" % (
                    platform.system(), prog, suff_bit))))

    def tearDown(self):
        self.x.options.alignment_file.close()
        self.x.options.sequence_file.close()
        self.x.options.tree_file.close()
        sepp.scheduler._jobPool = None
        shutil.rmtree(self.x.options.outdir, ignore_errors=True)

    def test_id_collision_working(self):

        self.x.run()
        self.assertTrue(self.x.results is not None)
        assert(len(self.x.results) == 490)
        assert(300 < len(self.x.results['SEQ396']) < 600)
        assert(len(self.x.results['SEQ554'].replace('-', '')) == 57)


if __name__ == "__main__":
    unittest.main()
