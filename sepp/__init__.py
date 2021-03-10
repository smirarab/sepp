###########################################################################
#    Copyright 2012 Siavash Mirarab, Nam Nguyen, and Tandy Warnow.
#    This file is part of SEPP.
#
#    SEPP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SEPP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with SEPP.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################

from operator import itemgetter
import logging
import os


__all__ = ['algorithm', 'alignment', 'backtranslate',
           'checkpointing', 'config', 'decompose_tree', 'ensemble',
           'exhaustive', 'exhaustive_upp', 'filemgr',
           'jobs', 'math_utils', 'problem', 'scheduler',
           'scratch', 'tree', 'get_logger', 'is_temp_kept', 'version']

version = "4.5.1"
_DEBUG = ("SEPP_DEBUG" in os.environ) and \
    (os.environ["SEPP_DEBUG"].lower() == "true")

_INSTALL_PATH = __path__[0]


def is_temp_kept():
    return _DEBUG


def get_setup_path():
    return _INSTALL_PATH


def get_logging_level():
    return logging.DEBUG if _DEBUG else logging.INFO


__set_loggers = set()


def get_logger(name="sepp"):
    logger = logging.getLogger(name)
    if name not in __set_loggers:
        level = get_logging_level()
        logging_formatter = logging.Formatter(
            ("[%(asctime)s] %(filename)s (line %(lineno)d):"
             " %(levelname) 8s: %(message)s"))
        logging_formatter.datefmt = '%H:%M:%S'
        logger.setLevel(level)
        ch = logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(logging_formatter)
        logger.addHandler(ch)
        __set_loggers.add(name)
    return logger


def reset_loggers():
    global __set_loggers
    __set_loggers = set()
    import pkgutil
    import sepp
    package = sepp
    for modl, name, _ in pkgutil.iter_modules(package.__path__):
        logger = (getattr(getattr(sepp, name, None), "_LOG", None))
        print("--- *", name, logger)
        if logger:
            setattr(getattr(sepp, name, None), "_LOG", get_logger(
                "sepp.%s" % name))


def log_exception(logger):
    """Logs the exception trace to the logObj as an error"""
    import traceback
    import io
    s = io.StringIO()
    traceback.print_exc(None, s)
    logger.debug(s.getvalue())


os.sys.setrecursionlimit(1000000)


def sort_by_value(d, reverse=False):
    return sorted(iter(d.items()), key=itemgetter(1), reverse=reverse)
