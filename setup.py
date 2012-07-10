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

#from distutils.core import setup
from distribute_setup import use_setuptools
use_setuptools(version="0.6.24")

from setuptools import setup, find_packages



setup(name = "sepp",
      version = "1.0",
      description = "SATe enabled phylogenetic placement.",
      packages = find_packages(),
      package_data = {'sepp' : ["*.py", "lib/merge.jar"] },

      url = "http://www.cs.utexas.edu/~phylo/software/sepp", 
      author = "Siavash Mirarab and Nam Nguyen",
      author_email = "smirarab@gmail.com, namphuon@cs.utexas.edu",

      license="General Public License (GPL)",
      install_requires = ["dendropy >= 3.4", "numpy >= 1.6","biopython >= 1.58"],
      provides = ["sepp"],
      scripts = ["sepp/scripts/run_sepp.py"],

      classifiers = ["Environment :: Console",
                     "Intended Audience :: Developers",
                     "Intended Audience :: Science/Research",
                     "License :: OSI Approved :: GNU General Public License (GPL)",
                     "Natural Language :: English",
                     "Operating System :: OS Independent",
                     "Programming Language :: Python",
                     "Topic :: Scientific/Engineering :: Bio-Informatics"])
