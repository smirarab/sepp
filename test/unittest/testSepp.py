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

from sepp.exhaustive import ExhaustiveAlgorithm
# sepp._DEBUG = True
# sepp.reset_loggers()


class Test(unittest.TestCase):
    x = None

    def setUp(self):
        # ensure necessary settings are made
        sepp.scheduler._jobPool = None
        sys.argv = [sys.argv[0], "-c", get_data_path("configs/test2.config")]
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

    def tearDown(self):
        self.x.options.alignment_file.close()
        self.x.options.info_file.close()
        self.x.options.tree_file.close()
        self.x.options.fragment_file.close()
        sepp.scheduler._jobPool = None
        shutil.rmtree(self.x.options.outdir, ignore_errors=True)

    def test_id_collision_working(self):
        self.x.options.fragment_file = open(
            get_data_path(
                "q2-fragment-insertion/input_fragments.fasta"), "r")
        self.x.run()
        self.assertTrue(self.x.results is not None)

    def test_id_collision_collision(self):
        self.x.options.fragment_file = open(
            get_data_path(
                "q2-fragment-insertion/input_fragments_collide.fasta"), "r")
        with self.assertRaisesRegex(
                ValueError,
                ' whose names overlap with names in your reference'):
            self.x.run()
        self.assertTrue(self.x.results is None)

    def test_seqnames_whitespaces(self):
        self.x.options.fragment_file = open(
            get_data_path(
                "q2-fragment-insertion/input_fragments_spaces.fasta"), "r")
        with self.assertRaisesRegex(
                ValueError,
                "contain either whitespaces: "):
            self.x.run()
        self.assertTrue(self.x.results is None)

    def test_fake_jobs(self):
        self.x.options.fragment_file = open(
            get_data_path(
                "q2-fragment-insertion/input_fragments_tiny.fasta"), "r")
        self.x.run()
        self.assertTrue(self.x.results is not None)

    def test_notpiped_jobs(self):
        sepp.config.options().hmmsearch.piped = "False"
        self.x.options.fragment_file = open(
            get_data_path(
                "q2-fragment-insertion/input_fragments_tiny.fasta"), "r")
        self.x.run()
        self.assertTrue(self.x.results is not None)


if __name__ == "__main__":
    unittest.main()
