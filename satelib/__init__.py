'''Functions for configuring a runtime logger for the satelib module.
'''
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

# Jiaye Yu and Mark Holder, University of Kansas

# This file is copied to SEPP and is used as an external library

PROGRAM_NAME = "SATe"
PROGRAM_AUTHOR = ["Jiaye Yu", "Mark Holder", "Jeet Sukumaran"]
PROGRAM_LICENSE = "GNU General Public License, version 3"
PROGRAM_VERSION = "1.3.0"
PROGRAM_YEAR = "2009-2011"
PROGRAM_DESCRIPTION = "Simultaneous Alignment and Tree Estimation"
PROGRAM_WEBSITE = "http://phylo.bio.ku.edu/software/sate/sate.html"
PROGRAM_INSTITUTE = "Department of Ecology and Evolutionary Biology, Univesity of Kansas"
PROGRAM_LONG_DESCRIPTION = """
SATe performs iterative realignment and tree inference.

Minimally you must provide a sequence file (with the '--input' option) a starting tree is optional.

The command line allows you to alter the behavior of the algorithm (termination criteria, when the algorithm switches to "Blind" acceptance of new alignments, how the tree is decomposed to find subproblems to be used, and the external tools to use).

Options can also be passed in as configuration files.

With the format:
####################################################
[commandline]
option-name = value

[sate]
option-name = value
####################################################

If you tell sate to keep its temporary files (-k option), then the configuration for the run will be stored in a file "last_used.cfg" in the "temporary" directory.

If configuration files are read in the order they occur as arguments (with values in later files replacing previously read values). Options specified in the command line are read last. Thus these values "overwrite" any settings from the configuration files.
"""


import os
import glob
import logging
import re
import time
import sys


class Messenger(object):
    """
    Wraps reporting of messages, progress, warnings and errors to users.
    Singleton (instantiated below).
    """

    def __init__(self):
        self.err_log_streams = [sys.stderr]
        self.run_log_streams = [sys.stdout]

    def _format_msg(self, msg, msg_type):
        return "SATe %s: %s\n" % (msg_type, msg)

    def _write_to_streams(self, streams, msg, flush=True):
        for s in streams:
            s.write(msg)
            if flush:
                s.flush()

    def _write_to_err_streams(self, msg, flush=True):
        self._write_to_streams(self.err_log_streams, msg, flush=flush)

    def _write_to_out_streams(self, msg, flush=True):
        self._write_to_streams(self.run_log_streams, msg, flush=flush)

    def send_error(self, msg):
        msg = self._format_msg(msg, "ERROR")
        self._write_to_err_streams(msg)

    def send_warning(self, msg):
        msg = self._format_msg(msg, "WARNING")
        self._write_to_err_streams(msg)

    def send_info(self, msg):
        msg = self._format_msg(msg, "INFO")
        self._write_to_out_streams(msg)

##############################################
## Instantiation Of Singleton
## more idiomatic way might be to implement
## this functionality as a module instead of
## a class
MESSENGER = Messenger()

