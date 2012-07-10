#!/usr/bin/env python
from sepp import get_logger

# This file is part of SATe

# SATe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Mark Holder, University of Kansas

"""
File and path management.
"""

import os
import re
import tempfile
import shutil
from threading import Lock
_LOG = get_logger(__name__)

_ILLEGAL_FILENAME_PATTERN = re.compile(r'[^-_a-zA-Z0-9.]')
def get_safe_filename(filename):
    return "".join(_ILLEGAL_FILENAME_PATTERN.split(filename))

def open_with_intermediates(filepath, mode):
    """Opens a `filepath` in the `mode`, but will create intermediate directories if they are absent."""
    d = os.path.dirname(filepath)
    if d:
        if not os.path.exists(d):
            os.makedirs(d)
        elif not os.path.isdir(d):
            raise IOError('The file "%s" cannot be created because "%s" exists but is not a directory' % (filepath, d))
    return open(filepath, mode)

class TempFS(object):
    '''A wrapper for creating temporary directories that will safeguard against
    removing directories that were not created by SATe (there is no evidence
    that this has happened, but rmtree is a dangerous call and it would be a
    horrible bug to have).

    Note that this class does not protect against incorrect order of deletion.
        If you delete the top level tempdirectory then all of the subdirectories
        will be deleted.
    '''

    def __init__(self):
        self._directories_created = set()
        self._top_level_temp = None
        self._top_level_temp_real = None
        self._directories_created_lock = Lock()

    def _is_already_created(self, real_path):
        self._directories_created_lock.acquire()
        b = real_path in self._directories_created
        self._directories_created_lock.release()
        return b


    def create_subdir(self, dir):
        '''Creates a directory `dir`

        `dir` must be a subdirectory of self.top_level_temp and `dir`
        cannot exist.

        The canonical file path is returned.

        Client code is responsible for deleting directory by a call to
        `remove_dir` using the same `TempFS` instance.

        '''
        rp = os.path.realpath(dir)
        if os.path.exists(rp):
            raise OSError("Path exists: '%s'" % rp)
        if not self._is_in_top_level_temp(rp):
            raise OSError("Subdirectory is not under the top level temp dir: '%s'" % rp)
        self._directories_created_lock.acquire()
        try:
            already_in = rp in self._directories_created
            if already_in:
                raise OSError("Subdirectory is flagged as having already been created: '%s'" % rp)
            else:
                os.makedirs(rp)
                self._directories_created.add(rp)
                return rp
        finally:
            self._directories_created_lock.release()


    def create_top_level_temp(self, parent, prefix='temp'):
        '''Creates (and stores the path to) the top-level temporary
        directory under `parent`

        `parent` must already exist
        `prefix` is passed to tempfile.mkdtemp

        Client code is responsible for deleting directory by a call to
        `remove_dir` using the same `TempFS` instance.

        The canonical file path is returned.

        '''
        assert(self._top_level_temp is None)
        assert(self._top_level_temp_real is None)
        r_parent = os.path.realpath(parent)
        if not os.path.exists(r_parent):
            raise OSError("Path does not exist: '%s'" % r_parent)
        if not os.path.isdir(r_parent):
            raise OSError("Path is not a directory: '%s'" % r_parent)
        self._directories_created_lock.acquire()
        try:
            self._top_level_temp = tempfile.mkdtemp(prefix=prefix, dir=r_parent)
            self._top_level_temp_real = os.path.realpath(self._top_level_temp)
            self._directories_created.add(self._top_level_temp_real)
            return self._top_level_temp_real
        finally:
            self._directories_created_lock.release()


    def create_temp_subdir(self, parent, prefix='temp'):
        '''Creates (and stores the path to) a temporary directory under `parent`

        `parent` should be "below" the top_level_temp.
        `prefix` is passed to tempfile.mkdtemp

        Client code is responsible for deleting directory by a call to
        `remove_dir` using the same `TempFS` instance.

        The canonical file path is returned.

        '''
        r_parent = os.path.realpath(parent)
        if not os.path.exists(r_parent):
            raise OSError("Path does not exist: '%s'" % r_parent)
        if not os.path.isdir(r_parent):
            raise OSError("Path is not a directory: '%s'" % r_parent)
        if not self._is_in_top_level_temp(r_parent):
            raise OSError("Subdirectory is not under the top level temp dir: '%s'" % r_parent)
        self._directories_created_lock.acquire()
        try:
            d = tempfile.mkdtemp(prefix=prefix, dir=r_parent)
            rp = os.path.realpath(d)
            self._directories_created.add(rp)
            return rp
        finally:
            self._directories_created_lock.release()

    def remove_dir(self, real_path):
        '''Recursively removes `real_path` from the filesystem if it is
        listed as one of the directories created by this TempFS object (or raises
        a ValueError if it is not listed).
        '''
        self._directories_created_lock.acquire()
        try:
            if real_path in self._directories_created:
                self._directories_created.remove(real_path)
                if (real_path == self._top_level_temp) or (real_path == self._top_level_temp_real):
                    self._top_level_temp = None
                    self._top_level_temp_real = None
            else:
                raise ValueError("'%s' is not registered as a temporary directory that was created by this process!" % real_path)
        finally:
            self._directories_created_lock.release()
        _LOG.debug("Removing temp dir: '%s'" % real_path)
        # Because we raise an exception if real_path is not in _directories_created,
        #   we only get down here if real_path was in self._directories_created
        #   thus, this call should only delete a directory created by this
        #   TempFS instance.
        if os.path.exists(real_path):
            try:
                shutil.rmtree(real_path)
                return True
            except:
                return False
        return False

    def get_remaining_directories(self):
        '''Returns a copy of the set of directories that have been created but
        not deleted.
        '''
        self._directories_created_lock.acquire()
        c = set(self._directories_created)
        self._directories_created_lock.release()
        return c


    def get_top_level_temp(self):
        return self._top_level_temp_real
    top_level_temp = property(get_top_level_temp)



    def _is_in_top_level_temp(self, real_path):
        if self._top_level_temp_real is None:
            raise ValueError("_top_level_temp has not been set, yet!")
        in_common = os.path.commonprefix([self._top_level_temp_real, real_path])
        return in_common == self._top_level_temp_real

temp_fs = TempFS()