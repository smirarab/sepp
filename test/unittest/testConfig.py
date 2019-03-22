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


class Test(unittest.TestCase):
    def testConfigFile(self):
        # Just to make different test cases independent of each other.
        config._options_singelton = None
        back = config.main_config_path
        # Diasable main config path for this test
        config.main_config_path = os.path.expanduser(
            "~/.sepp/main.config.notexistentfile")
        sys.argv = [sys.argv[0], "-A", "2", "-c", "data/configs/test.config",
                    "--outdir", "dir_form_commandline"]

        assert options().alignment_size == 2,
        "Commandline option -A not read properly"

        assert isinstance(options().config_file, filetypes)
        and options().config_file.name == "data/configs/test.config",
        "Commandline option -c not read properly"

        assert (options().pplacer is not None
                and options().pplacer.path == "pplacer"),
        "config file options not read properly"

        assert options().placement_size == 10,
        "Config file option placementSize not read properly"

        assert options().outdir.endswith("dir_form_commandline"),
        "Config file value outdir is not properly overwritten:%s " %
        options().outdir

        assert options().tempdir is not None,
        "Default value not properly set for tempfile attribute"

        config.main_config_path = back

    def testConfigFileMissingFile(self):
        # Just to make different test cases independent of each other.
        config._options_singelton = None
        back = config.main_config_path
        # Diasable main config path for this test
        config.main_config_path = os.path.expanduser(
            "~/.sepp/main.config.notexistentfile")

        sys.argv = [sys.argv[0], "-c", "data/configs/test2.config", "-f",
                    "data/simulated/test.fas"]
        assert isinstance(options().config_file, filetypes)
        and options().config_file.name == "data/configs/test2.config",
        "Commandline option -c not read properly"

        assert isinstance(options().alignment_file, filetypes)
        and options().alignment_file.name == "data/simulated/test.small.fas",
        "Config file option alignment not read properly"

        assert isinstance(options().fragment_file, filetypes)
        and options().fragment_file.name == "data/simulated/test.fas",
        "Command-line option -f alignment not read properly"

        config.main_config_path = back

    def testMainConfigFile(self):
        # Just to make different test cases independent of each other.
        config._options_singelton = None

        sys.argv = [sys.argv[0]]

        assert (options().pplacer is not None and os.path.exists(
            options().pplacer.path)),
        ("main config file options not read properly, or nonexistent binaries:"
         " pplacer = %s" % options().pplacer.path)

        assert (options().hmmalign is not None and os.path.exists(
            options().hmmalign.path)),
        ("main config file options not read properly, or nonexistent "
         "binaries: hmmalign = %s" % options().pplacer.path)

        assert (options().hmmsearch is not None and os.path.exists(
            options().hmmalign.path)),
        ("main config file options not read properly, or nonexistent bina"
         "ries: hmmsearch = %s" % options().pplacer.path)

    def testCpuCount(self):
        # Just to make different test cases independent of each other.
        config._options_singelton = None
        back = config.main_config_path
        # Diasable main config path for this test
        config.main_config_path = os.path.expanduser(
            "~/.sepp/main.config.notexistentfile")
        sys.argv = [sys.argv[0], "-x", "7"]

        assert options().cpu == 7, "Commandline option -x not read properly"

        config.main_config_path = back


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testConfigFile']
    unittest.main()
