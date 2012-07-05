#!/usr/bin/env python

###########################################################################
##    Copyright 2010 Rahul Suri and Tandy Warnow.
##    This file is part of ReUP.
##
##    ReUP is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    ReUP is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with ReUP.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################

from distutils.core import setup

setup(name = "sepp",
      version = "1.0",
      description = "SATe enabled phylogenetic placement.",
      packages = ["sepp", "satelib"],
      package_data = {'sepp' : ["sepp/*.py", "sepp/tools/merge.*"] },

      url = "http://www.cs.utexas.edu/~phylo/software/sepp", 
      author = "Siavash Mirarab and Nam Nguyen",
      author_email = "smirarab@gmail.com, namphuon@cs.utexas.edu",

      license="General Public License (GPL)",
      requires = ["dendropy (>= 3.4)", "bio (>= 1.0)"],
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

