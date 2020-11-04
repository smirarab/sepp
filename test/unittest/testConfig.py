'''
Created on Oct 2, 2012

@author: smirarab
'''
import unittest
import sys

from sepp.config import options
from sepp import config
import os
import io
try:
    filetypes = (io.IOBase, file)
except NameError:
    filetypes = io.IOBase
from sepp.filemgr import get_data_path
from tempfile import mkstemp
from sepp.scheduler import JobPool
from multiprocessing import cpu_count


class Test(unittest.TestCase):
    fp_config = None

    def setUp(self):
        _, self.fp_config = mkstemp()

    def tearDown(self):
        os.remove(self.fp_config)

    def testConfigFile(self):
        # Just to make different test cases independent of each other.
        config._options_singelton = None
        # Diasable main config path for this test
        config.main_config_path = self.fp_config

        sys.argv = [
            sys.argv[0], "-A", "2",
            "-c", get_data_path("configs/test.config"),
            "--outdir", "dir_form_commandline"]

        assert options().alignment_size == 2, \
            "Commandline option -A not read properly"

        assert isinstance(options().config_file, filetypes) and \
            options().config_file.name.endswith("data/configs/test.config"), \
            "Commandline option -c not read properly"

        assert (options().pplacer is not None and
                options().pplacer.path == "pplacer"), \
            "config file options not read properly"

        assert options().placement_size == 10, \
            "Config file option placementSize not read properly"

        assert options().outdir.endswith("dir_form_commandline"), \
            "Config file value outdir is not properly overwritten:%s " % \
            options().outdir

        assert options().tempdir is not None, \
            "Default value not properly set for tempfile attribute"

    def testConfigFileMissingFile(self):
        # Just to make different test cases independent of each other.
        config._options_singelton = None
        # Diasable main config path for this test
        config.main_config_path = self.fp_config

        sys.argv = [sys.argv[0],
                    "-c", get_data_path("configs/test2.config"),
                    "-f", get_data_path("simulated/test.fas"),
                    "-a", get_data_path("simulated/test.small.fas")]
        assert isinstance(options().config_file, filetypes) and \
            options().config_file.name.endswith(
                "data/configs/test2.config"), \
            "Commandline option -c not read properly"

        assert isinstance(options().alignment_file, filetypes) and\
            options().alignment_file.name.endswith(
                "data/simulated/test.small.fas"), \
            "Config file option alignment not read properly"

        assert isinstance(options().fragment_file, filetypes) and\
            options().fragment_file.name.endswith(
                "data/simulated/test.fas"), \
            "Command-line option -f alignment not read properly"

    def testMainConfigFile(self):
        # Just to make different test cases independent of each other.
        config._options_singelton = None

        sys.argv = [sys.argv[0], "-c", get_data_path("configs/test2.config")]
        # set pplacer filepath to a file shipped with the code base
        options().pplacer.path = get_data_path(
            "../../../tools/bundled/Darwin/pplacer")

        assert (options().pplacer is not None and os.path.exists(
                options().pplacer.path)), \
            ("main config file options not read properly,"
             "or nonexistent binaries: pplacer = %s" %
             options().pplacer.path)

        options().hmmalign.path = get_data_path(
            "../../../tools/bundled/Darwin/hmmalign")
        assert (options().hmmalign is not None and os.path.exists(
                options().hmmalign.path)), \
            ("main config file options not read properly, or nonexistent "
             "binaries: hmmalign = %s" % options().hmmalign.path)

        options().hmmsearch.path = get_data_path(
            "../../../tools/bundled/Darwin/hmmsearch")
        assert (options().hmmsearch is not None and os.path.exists(
                options().hmmsearch.path)), \
            ("main config file options not read properly, or nonexistent "
             "binaries: hmmsearch = %s" % options().hmmsearch.path)

    def testCpuCount(self):
        # Just to make different test cases independent of each other.
        config._options_singelton = None
        # Disable main config path for this test
        config.main_config_path = self.fp_config
        JobPool().terminate()
        JobPool().__init__(7)
        sys.argv = [sys.argv[0], "-x", "7"]

        assert options().cpu == 7, "Commandline option -x not read properly"

        # clean up after test:
        # 1) the JobPool CPU counts needs to be reset to the default
        # 2) the command line arguments must be restored
        JobPool().terminate()
        JobPool().__init__(cpu_count())
        sys.argv = [sys.argv[0], "-x", str(cpu_count())]
        config._options_singelton = None
        options()

    def testLog(self):

        import logging
        import sepp.jobs

        sdb = sepp._DEBUG

        sepp._DEBUG = True
        sepp.reset_loggers()
        sepp.jobs._LOG.debug("test debugging works")
        assert(sepp.jobs._LOG.getEffectiveLevel() == logging.DEBUG)

        sepp._DEBUG = False
        sepp.reset_loggers()
        sepp.jobs._LOG.debug("test debugging is disabled")
        assert(sepp.jobs._LOG.getEffectiveLevel() == logging.INFO)

        sepp._DEBUG = sdb
        sepp.reset_loggers()


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testConfigFile']
    unittest.main()
