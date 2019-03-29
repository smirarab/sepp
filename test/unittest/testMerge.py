'''
Created on Nov 20, 2012

@author: smirarab
'''
import unittest
from sepp.jobs import MergeJsonJob
import sys
import os.path
from sepp.filemgr import get_data_path
import sepp.config
import logging
from argparse import Namespace


logging.disable(logging.INFO)

class Test(unittest.TestCase):
    def testMerge(self):
        sys.argv = [sys.argv[0]]
        stdindata = open(get_data_path("tmp/tempmerge")).read()

        # write path of seppJsonMerger.jar into configuration
        setattr(sepp.config.options(), 'jsonmerger', Namespace(
            path=get_data_path("../../../tools/merge/seppJsonMerger.jar")))
        mergeJsonJob = MergeJsonJob()
        mergeJsonJob.setup(stdindata.replace(
            'data/tmp/pplacer.extended',
            os.path.dirname(get_data_path("tmp/tempmerge")) +
            '/pplacer.extended'), get_data_path("tmp/mergedfile"))
        mergeJsonJob.run()


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testMerge']
    unittest.main()
