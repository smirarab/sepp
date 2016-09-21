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
from distutils.core import setup, Command
from distutils.command.install import install

version = "3.1"

class ConfigSepp(Command):
    """setuptools Command"""
    description = "Configures Sepp for the current user"
    user_options = []
    
    def initialize_options(self):
        """init options"""
        self.configfile = os.path.expanduser("~/.sepp/main.config")
        pass

    def finalize_options(self):
        """finalize options"""
        pass
    
    def run(self):        
        print "\nCreating main sepp config file at %s " %(self.configfile)            
        def get_tools_dir(where):    
            platform_name = platform.system()
            if where is None:
                where = os.path.join("bundled",platform_name)
            path = os.path.join(os.getcwd(),"tools",where)
            if not os.path.exists(path):
                raise OSError("SEPP does not bundle tools for '%s' at this time!" % platform_name)
            return path
    
        def get_tool_name(tool,bits):
            if platform.system() == "Darwin" or not bits:#MAC doesn't have 32/64
                return tool
            is_64bits = sys.maxsize > 2**32
            return "%s-%s" %(tool,"64" if is_64bits else "32")
        
        def get_tools_dest():
            return os.path.join(os.path.dirname(self.configfile),"bundled-v%s"%version) 
        
        def copy_tool_to_lib(tool,where=None,bits=True):    
            shutil.copy2(os.path.join(get_tools_dir(where),get_tool_name(tool,bits)), 
                        os.path.join(get_tools_dest(),tool))
                    
        # Create the default config file
        if not os.path.exists(os.path.expanduser("~/.sepp")):
            os.mkdir(os.path.expanduser("~/.sepp"))
        if not os.path.exists(get_tools_dest()):
            os.mkdir(get_tools_dest())
        c = open("default.main.config")
        d = open(self.configfile,"w")
        for l in c:
            l = l.replace("~",get_tools_dest())
            d.write(l)
        d.close()
    
        # Copy tools to a bundled directory inside .sepp
        copy_tool_to_lib("guppy")
        copy_tool_to_lib("pplacer")
        copy_tool_to_lib("hmmalign")
        copy_tool_to_lib("hmmsearch")
        copy_tool_to_lib("hmmbuild")
        #TODO: should we compile and build merge.jar?
        copy_tool_to_lib("seppJsonMerger.jar",where="merge",bits=False)
        
class ConfigUPP(Command):
    """setuptools Command"""
    description = "Configures UPP for the current user"
    user_options = []
    
    def initialize_options(self):
        """init options"""
        self.configfile = os.path.expanduser("~/.sepp/upp.config")
        pass

    def finalize_options(self):
        """finalize options"""
        pass
    
    def run(self):                
        print "\nCreating main upp config file at %s " %(self.configfile)
        def get_tools_dir(where):    
            platform_name = platform.system()
            if where is None:
                where = os.path.join("bundled",platform_name)
            path = os.path.join(os.getcwd(),"tools",where)
            if not os.path.exists(path):
                raise OSError("SEPP does not bundle tools for '%s' at this time!" % platform_name)
            return path
    
        def get_tool_name(tool,bits):
            if platform.system() == "Darwin" or not bits:#MAC doesn't have 32/64
                return tool
            is_64bits = sys.maxsize > 2**32
            return "%s-%s" %(tool,"64" if is_64bits else "32")
        
        def get_tools_dest():
            return os.path.join(os.path.dirname(self.configfile),"bundled-v%s"%version) 
        
        def copy_tool_to_lib(tool,where=None,bits=True):    
            shutil.copy2(os.path.join(get_tools_dir(where),get_tool_name(tool,bits)), 
                        os.path.join(get_tools_dest(),tool))
                    
        # Create the default config file
        c = open("default.main.config")
        d = open(self.configfile,"w")
        for l in c:
            l = l.replace("~",get_tools_dest())
            d.write(l)
        d.write('\n[pasta]\npath=run_pasta.py\n')            
        d.close()        


class ConfigTIPP(Command):
    """setuptools Command"""
    description = "Configures TIPP for the current user"
    user_options = []
    
    def initialize_options(self):
        """init options"""
        self.configfile = os.path.expanduser("~/.sepp/tipp.config")
        pass

    def finalize_options(self):
        """finalize options"""
        pass
    
    def run(self):                
        print "\nCreating main tipp config file at %s " %(self.configfile)
        def get_tools_dir(where):    
            platform_name = platform.system()
            if where is None:
                where = os.path.join("bundled",platform_name)
            path = os.path.join(os.getcwd(),"tools",where)
            if not os.path.exists(path):
                raise OSError("SEPP does not bundle tools for '%s' at this time!" % platform_name)
            return path
    
        def get_tool_name(tool,bits):
            if platform.system() == "Darwin" or not bits:#MAC doesn't have 32/64
                return tool
            is_64bits = sys.maxsize > 2**32
            return "%s-%s" %(tool,"64" if is_64bits else "32")
        
        def get_tools_dest():
            return os.path.join(os.path.dirname(self.configfile),"bundled-v%s"%version) 
        
        def copy_tool_to_lib(tool,where=None,bits=True):    
            shutil.copy2(os.path.join(get_tools_dir(where),get_tool_name(tool,bits)), 
                        os.path.join(get_tools_dest(),tool))
                    
        # Create the default config file
        c = open("default.main.config")
        d = open(self.configfile,"w")
        for l in c:
            l = l.replace("~",get_tools_dest())
            if (l.find('seppJsonMerger.jar') != -1):
              l=l.replace('seppJsonMerger.jar','tippJsonMerger.jar')
            d.write(l)
        if not os.getenv('SATE') is None:
            d.write('\n[sate]\npath=%s' % os.getenv('SATE'))
        if os.getenv('BLAST') is None:
            print "\nWarning! BLAST variable is not defined.  If you plan to run TIPP for abundance profiling, then have BLAST pointed to blastn executable.  You can also change your config to point to blastn by including the following line in your config:\n[blast]\npath=/location/of/blast_directory/blastn\n"            
        if os.getenv('REFERENCE') is None:
            print "\nWarning! REFERENCE variable is not defined.  If you plan to run TIPP for abundance profiling, then have REFERENCE pointed to Reference directory.  You can also change your config to point to the Reference directory by including the following line in your config:\n[reference]\npath=/location/of/reference_directory/\n"            
        d.write('\n[blast]\npath=%s\n' % os.getenv('BLAST'))
        d.write('\n[reference]\npath=%s\n' % os.getenv('REFERENCE'))
        d.write('\n[tipp]\npushdown = true\n')                    
        d.close()        

    
        # Copy tools to a bundled directory inside .sepp
        copy_tool_to_lib("tippJsonMerger.jar",where="merge",bits=False)
    
setup(name = "sepp",
      version = version,
      description = "SATe enabled phylogenetic placement.",
      packages = find_packages(),

      url = "https://github.com/smirarab/sepp", 
      author = "Siavash Mirarab and Nam Nguyen",
      author_email = "smirarab@gmail.com, namphuon@cs.utah.edu",

      license="General Public License (GPL)",
      install_requires = ["dendropy >= 3.4, <= 3.99"],
      provides = ["sepp"],
      scripts = ["run_sepp.py",'run_tipp.py','run_upp.py','run_abundance.py',"split_sequences.py","run_tipp_tool.py"],
      cmdclass = {"config": ConfigSepp,"tipp": ConfigTIPP,"upp":ConfigUPP},
      
      classifiers = ["Environment :: Console",
                     "Intended Audience :: Developers",
                     "Intended Audience :: Science/Research",
                     "License :: OSI Approved :: GNU General Public License (GPL)",
                     "Natural Language :: English",
                     "Operating System :: OS Independent",
                     "Programming Language :: Python",
                     "Topic :: Scientific/Engineering :: Bio-Informatics"])
