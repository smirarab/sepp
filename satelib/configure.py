#!/usr/bin/env python

"""Accessor for runtime configuration settings.

In general one should be able to simply call get_configuration()

"""
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

import os
import platform
import sys

from satelib import get_logger
from satelib.settings import SateUserSettings
from satelib.tools import get_external_tool_classes

_LOG = get_logger(__name__)

def get_invoke_run_sate_command():
    """Used by GUI's this tries to return the python invocation for the run_sate
    program or script.
    """
    if sate_is_frozen():
        if platform.system() == 'Windows':
            return ['run_sate.exe']
        elif platform.system() == 'Darwin':
            return [sys.executable, os.path.join(sate_home_dir(), 'run_sate.py')]
        else:
            raise OSError('SATe is not frozen for %s' % platform.system())
    else:
        return [sys.executable, 'run_sate.py']

def sate_is_frozen():
    """Will return True if SATe is frozen.
    """
    import imp
    return (
        hasattr(sys, "frozen")          # new py2exe
        or hasattr(sys, "importers")    # old py2exe
        or imp.is_frozen("__main__")    # tools/freeze
    )

def sate_home_dir():
    """Attempts to return the directory that is the parent of the binary directories.
    """
    if sate_is_frozen():
        retpath = os.path.join( os.path.dirname(os.path.dirname(os.path.abspath(sys.executable))), 'Resources') if platform.system() == "Darwin" else os.path.dirname(sys.executable)
        return retpath
    else:
        # configure.py is expected to be in satelib directory of SATe main directory
        return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

_DEFAULT_TOOLS_PATH = None

def init_sate(sate_home=None):
    """Sets the _DEFAULT_TOOLS_PATH based on `sate_home`.
    Configuration objects created after this call will have paths to tools
    initially set based on this `sate_home`

    If `sate_home` is None, then
    """
    global _DEFAULT_TOOLS_PATH
    base_dir = sate_home_dir()
    if platform.system() == 'Windows':
        bin_dir = os.path.join(base_dir, 'bin_win')
    elif platform.system() == 'Darwin':
        bin_dir = os.path.join(base_dir, 'bin_mac')
    elif platform.system() == 'Linux':
        bin_dir = os.path.join(base_dir, 'bin_linux')
    else:
        raise OSError('SATe does not support %s at this time!' % platform.system())
    if _DEFAULT_TOOLS_PATH is None:
        _DEFAULT_TOOLS_PATH = {}
    for i in get_external_tool_classes():
        tool_name = i.section_name.split()[0]
        _DEFAULT_TOOLS_PATH[tool_name] = os.path.join(bin_dir, tool_name)
        if platform.system() == 'Windows' and tool_name != 'opal':
            _DEFAULT_TOOLS_PATH[tool_name] += '.exe'
    _DEFAULT_TOOLS_PATH['opal'] += '.jar'

def set_configuration_from_defaults(cfg):
    "Uses _DEFAULT_TOOLS_PATH to add paths to external tools to `cfg`"
    global _DEFAULT_TOOLS_PATH
    if _DEFAULT_TOOLS_PATH is None:
        init_sate()
    for name, path in _DEFAULT_TOOLS_PATH.items():
        x = getattr(cfg, name)
        x.path = path

def get_configuration(configfile=None):
    """Returns an instance of SateUserSettings that reflects the current
    defaults based on:
        1. paths inferred from sate_home_dir()
        2. any settings in `configfile` (these will overwrite the settings
            based on the defaults from sate_home).
    """
    cfg = SateUserSettings()
    set_configuration_from_defaults(cfg)
    if configfile is not None:
        if os.path.isfile(configfile):
            cfg.read_config_filepath(configfile)
        else:
            ### TODO: wrap up in messaging system
            sys.stderr.write("The specified configuration file %s cannot be found, the default settings are used instead.\n" % configfile)
    return cfg

def get_input_source_directory(config):
    """
    Given a configuration object, returns the directory of the input file(s).
    """
    options = config.commandline
    if options.multilocus:
        # multilocus dataset: assume directory is given as input source
        return os.path.abspath(options.input)
    else:
        # single locus dataset: return directory nanme
        return os.path.dirname(os.path.abspath(options.input))


