import os
import platform
import sys
import site

import shutil
from setuptools import Command, Distribution
import argparse
from sepp import version


def get_tools_dir(where):
    platform_name = platform.system()
    if where is None:
        where = os.path.join("bundled", platform_name)
    path = os.path.join(os.getcwd(), "tools", where)
    if not os.path.exists(path):
        raise OSError("SEPP does not bundle tools for '%s' at this time!" %
                      platform_name)
    return path


def get_tool_name(tool, bits):
    if platform.system() == "Darwin" or not bits:  # MAC doesn't have 32/64
        return tool
    is_64bits = sys.maxsize > 2**32
    return "%s-%s" % (tool, "64" if is_64bits else "32")


class ConfigSepp(Command):
    """setuptools Command"""
    def initopts(self, contained=None):
        self.contained = contained
        self.configfile = None
        self.basepath = None
        self.version = version

    def initpath(self, name):
        if self.contained:
            self.configfile = os.path.expanduser(
                os.path.abspath(os.path.join(".sepp", name)))
            self.basepath = os.path.dirname(self.configfile)
        else:
            self.configfile = os.path.expanduser("~/.sepp/%s" % name)
            self.basepath = os.path.expanduser("~/.sepp")
        fp_home_path = 'home.path'
        with open(fp_home_path, 'w') as fo:
            fo.write(self.basepath)
            fo.close()
        
        # copy created home.path file to site-packages directory
        target_dir = site.getsitepackages()[0]
        if not self.contained:
            target_dir = site.getusersitepackages()
        shutil.copy(fp_home_path, os.path.join(target_dir, fp_home_path))

    def get_tools_dest(self):
        return os.path.join(self.basepath, "bundled-v%s" % self.version)

    def copy_tool_to_lib(self, tool, where=None, bits=True):
        shutil.copy2(
            os.path.join(get_tools_dir(where), get_tool_name(tool, bits)),
            os.path.join(self.get_tools_dest(), tool))

    def initialize_options(self, contained=None):
        """init options"""
        self.initopts(contained)

    def finalize_options(self):
        """finalize options"""
        self.initpath("main.config")
        print("\nCreating main sepp config file at %s and tools at %s" % (
            self.configfile, self.basepath))

    def run(self):
        def get_tool_name(tool, bits):
            # MAC doesn't have 32/64
            if platform.system() == "Darwin" or not bits:
                return tool
            is_64bits = sys.maxsize > 2**32
            return "%s-%s" % (tool, "64" if is_64bits else "32")

        # Create the default config file
        if not os.path.exists(self.basepath):
            os.mkdir(self.basepath)
        if not os.path.exists(self.get_tools_dest()):
            os.mkdir(self.get_tools_dest())
        c = open("default.main.config")
        d = open(self.configfile, "w")
        for l1 in c:
            l1 = l1.replace("~", self.get_tools_dest())
            d.write(l1)
        d.close()

        # Copy tools to a bundled directory inside .sepp
        self.copy_tool_to_lib("guppy")
        self.copy_tool_to_lib("pplacer")
        self.copy_tool_to_lib("hmmalign")
        self.copy_tool_to_lib("hmmsearch")
        self.copy_tool_to_lib("hmmbuild")
        # TODO: should we compile and build merge.jar?
        self.copy_tool_to_lib("seppJsonMerger.jar", where="merge", bits=False)


class ConfigUPP(ConfigSepp):
    """setuptools Command"""
    def initialize_options(self, contained=None):
        """init options"""
        self.initopts(contained=contained)

    def finalize_options(self):
        """finalize options"""
        self.initpath("upp.config")
        print("\nCreating main UPP config file at %s and tools at %s" % (
            self.configfile, self.basepath))

    def run(self):
        # Create the default config file
        c = open("default.main.config")
        d = open(self.configfile, "w")
        for l2 in c:
            l2 = l2.replace("~", self.get_tools_dest())
            d.write(l2)
        d.write(('\n[pasta]\n'
                 'path=run_pasta.py\nuser_options= --max-mem-mb=2000'))
        d.close()


def config_sepp(cmd=ConfigSepp, name="SEPP"):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--contained",
        help=("Whether %s should be installed in a self-contained "
              "manner or on user's home") % name,
        action='store_true')
    args = parser.parse_args()

    dist = Distribution()
    conf = cmd(dist)
    conf.initialize_options(args.contained)
    conf.finalize_options()
    conf.run()


def config_upp():
    config_sepp(cmd=ConfigUPP, name="UPP")
