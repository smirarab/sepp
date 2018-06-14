'''
Created on June 14, 2018

@author: Stefan.M.Janssen@gmail.com
'''
import tempfile
import shutil
import unittest
import sepp

from sepp.exhaustive import ExhaustiveAlgorithm
sepp._DEBUG = True


class Test(unittest.TestCase):
    x = None

    def setUp(self):
        self.x = ExhaustiveAlgorithm()
        self.x.options.alignment_file = open(
            "data/q2-fragment-insertion/reference_alignment_tiny.fasta", "r")
        self.x.options.info_file = open(
            "data/q2-fragment-insertion/RAxML_info-reference-gg-raxml-bl.info",
            "r")
        self.x.options.tree_file = open(
            "data/q2-fragment-insertion/reference_phylogeny_tiny.nwk", "r")
        self.x.options.outdir = tempfile.mkdtemp()

    def tearDown(self):
        self.x.options.alignment_file.close()
        self.x.options.info_file.close()
        self.x.options.tree_file.close()
        self.x.options.fragment_file.close()
        shutil.rmtree(self.x.options.outdir, ignore_errors=True)

    def test_id_collision_working(self):
        self.x.options.fragment_file = open(
            "data/q2-fragment-insertion/input_fragments.fasta", "r")
        self.x.run()
        self.assertTrue(self.x.results is not None)

    def test_id_collision_collision(self):
        self.x.options.fragment_file = open(
            "data/q2-fragment-insertion/input_fragments_collide.fasta", "r")
        with self.assertRaisesRegex(
                ValueError,
                ' whose names overlap with names in your reference'):
            self.x.run()
        self.assertTrue(self.x.results is None)


if __name__ == "__main__":
    unittest.main()
