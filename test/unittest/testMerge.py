'''
Created on Nov 20, 2012

@author: smirarab
'''
import unittest
from sepp.jobs import MergeJsonJob
import sys


class Test(unittest.TestCase):
    def testMerge(self):
        sys.argv = [sys.argv[0]]
        stdindata = open("data/tmp/tempmerge").read()
        mergeJsonJob = MergeJsonJob()
        mergeJsonJob.setup(stdindata, "data/tmp/mergedfile")
        mergeJsonJob.run()


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testMerge']
    unittest.main()
