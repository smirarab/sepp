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

from operator import itemgetter    
import logging
import os

__all__ = ["alignment", "shortreadalignment", "taxonneighbourfinder", "tools", "problem"]

def sortByValue(d,reverse=False):
    return sorted(d.iteritems(), key=itemgetter(1), reverse=reverse)


_DEBUG = os.environ.has_key("SEPP_DEBUG") and os.environ["SEPP_DEBUG"].lower() == "true"
#print "Debug mode is %s." %("on" if _DEBUG else "off")

_INSTALL_PATH = __path__[0]

def is_temp_kept():
    return _DEBUG

def get_setup_path():
    return _INSTALL_PATH

def get_logging_level():    
    return logging.DEBUG

def get_logger(name="sepp"):
    logger_set = False
    logger = logging.getLogger(name)
    if not logger_set:
        level = get_logging_level()
        logging_formatter = logging.Formatter("[%(asctime)s] %(filename)s (line %(lineno)d): %(levelname) 8s: %(message)s")
        logging_formatter.datefmt = '%H:%M:%S'
        logger.setLevel(level)
        ch = logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(logging_formatter)
        logger.addHandler(ch)
    return logger

TIMING_LOG = logging.getLogger("TIMING_LOG")
TIMING_LOG.setLevel(logging.CRITICAL)

def set_timing_log_filepath(fp):
    global TIMING_LOG
    if not fp:
        TIMING_LOG.setLevel(logging.CRITICAL)
    else:
        TIMING_LOG.setLevel(logging.DEBUG)
        h = logging.FileHandler(fp)
        f = logging.Formatter("[%(asctime)s] : %(message)s")
        f.datefmt = '%D %H:%M:%S'
        h.setLevel(logging.DEBUG)
        h.setFormatter(f)
        TIMING_LOG.addHandler(h)

def log_exception(logger):
    '''Logs the exception trace to the logObj as an error'''
    import traceback, cStringIO
    s = cStringIO.StringIO()
    traceback.print_exc(None, s)
    logger.debug(s.getvalue())
