__author__ = 'Michael'
import sys,random,argparse,os,shutil
from argparse import ArgumentParser, Namespace
from sepp import get_logger
from sepp.algorithm import AbstractAlgorithm
from sepp.alignment import MutableAlignment, ExtendedAlignment,_write_fasta
from sepp.exhaustive import JoinAlignJobs, ExhaustiveAlgorithm
from sepp.jobs import PastaAlignJob
from sepp.filemgr import get_temp_file
from sepp.config import options,valid_decomp_strategy
import sepp.config
from sepp.math_utils import lcm
from sepp.problem import SeppProblem
from sepp.scheduler import JobPool
from multiprocessing import Pool, Manager
from sepp.alignment import ExtendedAlignment

_LOG = get_logger(__name__)

class HierarchicalAlgorithm(AbstractAlgorithm):
    '''
    Class to implement a hierarchical search as compared to the Exhaustive Algorithm
    '''
