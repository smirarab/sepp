from operator import itemgetter    
import logging
import os

__all__ = ["alignment", "shortreadalignment", "taxonneighbourfinder", "tools"]

def sortByValue(d,reverse=False):
    return sorted(d.iteritems(), key=itemgetter(1), reverse=reverse)


_LOGGING_LEVEL_ENVAR = "SATE_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR = "SATE_LOGGING_FORMAT"

_INSTALL_PATH = __path__[0]

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
