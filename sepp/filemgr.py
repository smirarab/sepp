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


from sepp import get_logger, is_temp_kept
import tempfile
import re
import sepp

"""
File and path management.
"""

import os
import shutil

_LOG = get_logger(__name__)

def open_with_intermediates(filepath, mode):
    """Opens a `filepath` in the `mode`, but will create intermediate directories if they are absent."""
    d = os.path.dirname(filepath)
    if d:
        if not os.path.exists(d):
            os.makedirs(d)
        elif not os.path.isdir(d):
            raise IOError('The file "%s" cannot be created because "%s" exists but is not a directory' % (filepath, d))
    return open(filepath, mode)

def get_temp_file(prefix, relative_path, suffix = ""):
    '''
    Returns a temporary directory given a relative_path (relative to root temp
    directory) and a suffix. Use tempfile_for_subproblem when possible. 
    '''
    d = os.path.join(get_root_temp_dir(), relative_path)
    if not os.path.exists(d):
        os.makedirs(d)
    return tempfile.mktemp(prefix = prefix,
                           suffix = suffix, 
                           dir = d)

def tempdir_for_subproblem(subproblem):
    '''
    Returns the path to a temporary directory for a given subproblem.
    '''
    largs = [x.label for x in reversed(subproblem.get_path_to_root())]
    return os.path.join(*largs)
    

def tempfile_for_subproblem(prefix, subproblem, suffix = ""):
    '''
    Returns the path to a temporary file for a given subproblem, and a prefix.
    '''
    from sepp.problem import SeppProblem
    assert isinstance(subproblem, SeppProblem)
    return get_temp_file(prefix, tempdir_for_subproblem(subproblem), suffix)
    

def is_my_temp_file(path):
    '''
    is this file part of my temp directory structure?
    '''
    return os.path.commonprefix([get_root_temp_dir(),path]) == get_root_temp_dir() 
                             
def remove_temp(path):
    '''
    Tries to remove a temp directory. Does so only if a temp file is created
    inside the tmp directory structure of current run.  
    '''
    if not is_my_temp_file(path):
        _LOG.warning("Temp File %s is not under temporary structure created by this job. It will not be removed. ")
        return  
    if not is_temp_kept() and os.path.exists(path):
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)
            
def get_default_temp_dir():
    ''' Returns a default temp directory. This is used internally, in case
    a temp directory is not provided by the user. Use get_root_temp_dir to find
    the current root temp directory. 
    '''
    tempfile.gettempdir()
    return os.path.join(tempfile.tempdir,"sepp")           

_root_temp_dir = None

def get_root_temp_dir():
    '''
    Return the root temp directory for the current run. Creates the temp file
    the first time it is called in any run. 
    '''
    global _root_temp_dir
    if _root_temp_dir is None:
        _root_temp_dir = tempfile.mkdtemp(prefix = sepp.config.options().output +".", 
                                          dir = sepp.config.options().tempdir)
        _LOG.info("Root temp directory built: %s" % _root_temp_dir)        
    return _root_temp_dir

def set_root_temp_dir(tmp):
    global _root_temp_dir
    _root_temp_dir = tmp

def check_or_make_dir_path(path):
    '''
    Make sure a path to a directory exists, and return the normalized path.
    If the path is to a file, or if it cannot be made, return None. 
    '''
    if not os.path.exists(path):
      #Race condition exists here
      try:
        os.mkdir(path)
      except:
        return os.path.normpath(os.path.abspath(path))
    if os.path.exists(path) and os.path.isdir(path):
        return os.path.normpath(os.path.abspath(path))        
    return None

def directory_has_files_with_prefix(d,prefix):
    '''
    Check to see if a directory has files with a given prefix
    '''
    files = os.listdir(d)
    for f in files:
        if (re.search('^' + prefix + '.*', f) is not None):
            return True
    return False
