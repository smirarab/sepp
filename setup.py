#!/usr/bin/env python

# ##########################################################################
#    Copyright 2012 Siavash Mirarab, Nam Nguyen, and Tandy Warnow.
#    This file is part of SEPP.
#
#    SEPP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SEPP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with SEPP.  If not, see <http://www.gnu.org/licenses/>.
# ##########################################################################

import os
import platform
import sys

# from distutils.core import setup
from distribute_setup import use_setuptools
import shutil
from setuptools import find_packages
from distutils.core import setup, Command

use_setuptools(version="0.6.24")
version = "4.5.1"


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
    description = "Configures Sepp for the current user"
    user_options = [(
        'contained', 'c',
        ("Whether SEPP should be installed in a self-contained "
         "manner or on user's home"))]

    def initopts(self):
        self.contained = None
        self.configfile = None
        self.basepath = None

    def initpath(self, name):
        if self.contained:
            self.configfile = os.path.expanduser(
                os.path.abspath(os.path.join(".sepp", name)))
            self.basepath = os.path.dirname(self.configfile)
        else:
            self.configfile = os.path.expanduser("~/.sepp/%s" % name)
            self.basepath = os.path.expanduser("~/.sepp")
        with open('home.path', 'w') as fo:
            fo.write(self.basepath)
            fo.close()

    def get_tools_dest(self):
        return os.path.join(self.basepath, "bundled-v%s" % version)

    def copy_tool_to_lib(self, tool, where=None, bits=True):
        shutil.copy2(
            os.path.join(get_tools_dir(where), get_tool_name(tool, bits)),
            os.path.join(self.get_tools_dest(), tool))

    def initialize_options(self):
        """init options"""
        self.initopts()

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
    description = "Configures UPP for the current user"
    user_options = [('contained', 'c',
                     ("Whether SEPP should be installed in a self-contained "
                      "manner or on user's home"))]

    def initialize_options(self):
        """init options"""
        self.initopts()

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


setup(name="sepp",
      version=version,
      description="SATe enabled phylogenetic placement.",
      packages=find_packages(),

      url="https://github.com/smirarab/sepp",
      author="Siavash Mirarab and Nam Nguyen",
      author_email="smirarab@gmail.com, namphuon@cs.utah.edu",

      license="General Public License (GPL)",
      install_requires=["dendropy >= 4.0.0"],
      provides=["sepp"],
      scripts=["run_sepp.py", 'run_upp.py', "split_sequences.py"],
      cmdclass={"config": ConfigSepp, "upp": ConfigUPP},
      data_files=[('', ['home.path'])],

      classifiers=["Environment :: Console",
                   "Intended Audience :: Developers",
                   "Intended Audience :: Science/Research",
                   ("License :: OSI Approved :: GNU General Public "
                    "License (GPL)"),
                   "Natural Language :: English",
                   "Operating System :: OS Independent",
                   "Programming Language :: Python",
                   "Topic :: Scientific/Engineering :: Bio-Informatics"])
