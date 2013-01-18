#!/usr/bin/env python

###########################################################################
##    Copyright 2012 Siavash Mirarab, Nam Nguyen, and Tandy Warnow.
##    This file is part of SEPP.
##
##    SEPP is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    SEPP is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with SEPP.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################

import os, platform, sys

#from distutils.core import setup
from distribute_setup import use_setuptools 
import shutil
use_setuptools(version="0.6.24")

from setuptools import setup, find_packages

def get_tools_dir():    
    platform_name = platform.system()
    path = os.path.join(os.getcwd(),"tools","bundled",platform_name)
    if not os.path.exists(path):
        raise OSError("SEPP does not bundle tools for '%s' at this time!" % platform_name)
    return path

def get_tool_name(tool):
    is_64bits = sys.maxsize > 2**32
    if platform.system() == "Darwin":#MAC doesn't have 32/64
        return tool
    return "%s-%s" %(tool,"64" if is_64bits else "32")

def copy_tool_to_lib(tool):    
    shutil.copy2(os.path.join(get_tools_dir(),get_tool_name(tool)), 
                os.path.join(os.getcwd(),"src","sepp","lib",tool))

if not os.path.exists(os.path.join(os.getcwd(),"src","sepp","lib")):
    os.mkdir(os.path.join(os.getcwd(),"src","sepp","lib"))

copy_tool_to_lib("guppy")
copy_tool_to_lib("pplacer")
copy_tool_to_lib("hmmalign")
copy_tool_to_lib("hmmsearch")
copy_tool_to_lib("hmmbuild")
#TODO: should we compile and build merge.jar?
shutil.copy2(os.path.join(os.getcwd(),"tools","merge","seppJsonMerger.jar"),
             os.path.join(os.getcwd(),"src","sepp","lib","seppJsonMerger.jar"))

a = setup(name = "sepp",
      version = "2.1",
      description = "SATe enabled phylogenetic placement.",
      packages = find_packages("src"),
      package_data = { "sepp" : ["lib/*"]},
      package_dir = {'':'src'},

      url = "http://www.cs.utexas.edu/~phylo/software/sepp", 
      author = "Siavash Mirarab and Nam Nguyen",
      author_email = "smirarab@gmail.com, namphuon@cs.utexas.edu",

      license="General Public License (GPL)",
      install_requires = ["dendropy >= 3.4"],
      provides = ["sepp"],
      scripts = ["run_sepp.py"],

      classifiers = ["Environment :: Console",
                     "Intended Audience :: Developers",
                     "Intended Audience :: Science/Research",
                     "License :: OSI Approved :: GNU General Public License (GPL)",
                     "Natural Language :: English",
                     "Operating System :: OS Independent",
                     "Programming Language :: Python",
                     "Topic :: Scientific/Engineering :: Bio-Informatics"])


if not os.path.exists(os.path.expanduser("~/.sepp")):
    os.mkdir(os.path.expanduser("~/.sepp"))

import sepp
lib=os.path.dirname(sepp.__file__)
c = open("default.main.config")
d = open(os.path.expanduser("~/.sepp/main.config"),"w")
for l in c:
    l = l.replace("~",lib)
    d.write(l)