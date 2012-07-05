#!/usr/bin/env python

"""Interface to external tools (alignment and tree inference programs)
"""

# This file is part of SATe, adopted by Siavash Mirarab for rpass

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

# Jiaye Yu and Mark Holder, University of Kansas

import os


from scheduler import jobq, start_worker, DispatchableJob
from sepp import get_logger
from sepp.settings import settings

_LOG = get_logger(__name__)

def is_file_checker(p):
    if not p:
        return False, "Expecting the path to an executable, got an empty string"
    if not os.path.isfile(p):
        return False, "%s is not a valid file" % str(p)
    return True, p

def is_executable_checker(p):
    r, msg = is_file_checker(p)
    if not r:
        return (r, msg)
    try:
        import stat
    except:
        return True, p
    if (os.stat(p)[stat.ST_MODE] & (stat.S_IXUSR|stat.S_IXGRP|stat.S_IXOTH)) == 0:
        return False, "%s does not appear to be an executable file" % p
    else:
        return True, p


class ExternalTool (object):
    def __init__(self, name):
        self.name = name
        args = settings.get_setting(self.name)
        self.exe = args['path']
        if not os.path.exists(self.exe):
            raise ValueError('The path "%s" does not exist' % self.exe)
        self.delete_temps = args.get('delete_temps', True)

    @staticmethod
    def exists(self):
        return is_executable_checker(self.exe)[0]

    def make_temp_workdir(self, tmp_dir_par):
        pref = 'temp' + self.name
        scratch_dir = self.temp_fs.create_temp_subdir(parent=tmp_dir_par, prefix=pref)
        return scratch_dir

    def run(self, *args, **kwargs):
        start_worker(1)
        job = self.create_job(*args, **kwargs)
        jobq.put(job)
        return job.get_results()

    def create_job(self, *args, **kwargs):
        raise NotImplementedError('Abstract ExternalTool class cannot spawn jobs.')

#class ShortReadAligner(ExternalTool):
#    def __init__(self, name, temp_fs, **kwargs):
#        ExternalTool.__init__(self, name, temp_fs, **kwargs)
#        self.user_opts = kwargs.get('args', ' ').split()
#
#    def _prepare_input(self, alignment, **kwargs):
#        """Wraps up the writing of raw fasta, creation of temp dir, ... for common aligners.
#        Returns directory, input filename, output filename."""
#        tdp = kwargs.get('tmp_dir_par')
#        if not tdp:
#            raise AssertionError('The tmp_dir_par must be specified when calling create_job or _prepare_input')
#        scratch_dir = self.make_temp_workdir(tmp_dir_par=tdp)
#        seqfn = os.path.join(scratch_dir, "input.fasta")
#        alignment.write_unaligned_fasta(seqfn)
#        alignedfn = os.path.join(scratch_dir, 'input.aligned')
#        return scratch_dir, seqfn, alignedfn
#
#    def create_job(self, *args, **kwargs):
#        raise NotImplementedError('Abstract Aligner class cannot spawn jobs.')
#
#    def _finish_standard_job(self, alignedfn, datatype, invoc, scratch_dir, job_id, delete_temps):
#        dirs_to_delete = []
#        if delete_temps:
#            dirs_to_delete = [scratch_dir]
#        # create a results processor to read the alignment file
#        rpc = lambda : read_internal_alignment(alignedfn, 
#                                               datatype=datatype,
#                                               dirs_to_delete=dirs_to_delete,
#                                               temp_fs=self.temp_fs)
#        job = DispatchableJob(invoc, result_processor=rpc,  cwd=scratch_dir, context_str=job_id)
#        return job
