"""
This module provides options(), which returns a argparse.Namespace object
containing all configurations required by sepp. config.options() reads and
saves the configurations the first time it is called. Subsequent calls simply
return the already-saved configurations. A Typical usage is:
some_config_attribute = sepp.config.options().some_config_attribute

All command-line options are found directly inside the return value of
options() as attributes. Commandline options can be specified either directly
using the command-line, or under [commandline] section in a config file,
passed to SEPP using -c commandline option. commandline values overwrite
config file values.

Information stored in other sections of the config file are available as
nested arparse.Namespaes attributes inside the results of options(),
with the config file header used as the attribute name. For example,
imagine the config file has:
[pplacer]
path = /some/path
In this case, config.options().pplacer will be a arparse.Namespace and
config.options().pplacer.path will be "/some/path".

A "main" configuration file under {home}/.sepp/main.config is used
to store some basic configurations such as the location of extra programs, etc.
This main config file is read first, so that user provided config file can
overwrite its values.

In addition, each client of this module (e.g. a new algorithm)) can add new
commandline options by getting the parser object using get_parser() and then
adding extra options. This has to happen *before* the first call to options()
module. For an example see exhaustive_upp.
"""

from argparse import ArgumentParser, Namespace
from sepp.filemgr import get_default_temp_dir, check_or_make_dir_path
import sys

try:
    import configparser
except ImportError:
    import ConfigParser as configparser
from sepp import version, get_logger
import argparse
import os
import os.path
from multiprocessing import cpu_count
from sepp import scheduler
import random
import importlib.metadata

_LOG = get_logger(__name__)

# obtain filepath of original sepp source directory, which is
# stored in a file "home.path", which is located in the site-package
# installation
fp_home = os.path.join(
    importlib.metadata.distribution("sepp")._path, "home.path")
root_p = open(fp_home).readlines()[0].strip()
print("root_p='%s'" % root_p)
main_config_path = os.path.join(root_p, "main.config")


def set_main_config_path(filename):
    global main_config_path
    main_config_path = filename


def _read_config_file(filename, opts, expand=None):
    _LOG.debug("Reading config %s" % filename)
    config_defaults = []
    cparser = configparser.ConfigParser()
    cparser.optionxform = str
    cparser.read_file(filename)

    if cparser.has_section('commandline'):
        for (k, v) in cparser.items('commandline'):
            config_defaults.append("--%s" % k)
            config_defaults.append(v)

    for section in cparser.sections():
        if section == "commandline":
            continue
        if getattr(opts, section, None):
            section_name_space = getattr(opts, section)
        else:
            section_name_space = Namespace()
        for (k, v) in cparser.items(section):
            if expand and k == "path":
                v = os.path.join(expand, v)
            section_name_space.__setattr__(k, v)
        opts.__setattr__(section, section_name_space)

    return config_defaults


def valid_dir_path(path):
    ret = check_or_make_dir_path(path)
    if ret is None:
        raise argparse.ArgumentTypeError(
            "%s is not a valid directory path." % path)
    return ret


def valid_molecule(molecule):
    ret = molecule in ['dna', 'rna', 'amino']
    if ret is False:
        raise argparse.ArgumentTypeError(
            ("%s is not a valid molecule type.  Must be 'dna', 'rna', or "
             "'amino'.") % molecule)
    return molecule


def valid_decomp_strategy(strategy):
    ret = strategy in ['hierarchical', 'normal', 'midpoint', 'centroid']
    if ret is False:
        raise argparse.ArgumentTypeError(
            ("%s is not a valid strategy.  Must be 'normal', 'hierarchical',"
             " 'centroid', or 'midpoint'.") % strategy)
    return strategy


def valid_file_prefix(prefix):
    if os.path.dirname(prefix) != "":
        raise argparse.ArgumentTypeError(
            "%s is not a valid output prefix (includes a directory)." % prefix)
    return prefix


def set_cpu(cpus):
    c = int(cpus)
    scheduler.default_cpus = c
    return c


def set_checkpoint(checkpoint):
    import sepp.checkpointing
    return sepp.checkpointing.CheckPointManager(checkpoint)


_parser = None


def _init_parser():
    global _parser
    _parser = ArgumentParser(
        description=(
            "This script runs the SEPP algorithm on an input "
            "tree, alignment, fragment file, and RAxML info file."),
        conflict_handler='resolve')

    _parser.add_argument("-v", "--version", action='version',
                         version="%(prog)s " + version)

    decompGroup = _parser.add_argument_group(
        "Decomposition Options".upper(),
        ' '.join([
            "These options determine the alignment decomposition size and",
            "taxon insertion size.  If None is given, then the default",
            "is to align/place at 10% of total taxa.  "
            "The alignment decomosition size must be",
            "less than the taxon insertion size."]))
    _parser.groups = dict()
    _parser.groups['decompGroup'] = decompGroup

    decompGroup.add_argument(
        "-A", "--alignmentSize", type=int,
        dest="alignment_size", metavar="N", default=None,
        help=("max alignment subset size of N "
              "[default: 10%% of the total number of taxa or the placement"
              " subset size if given]"))
    decompGroup.add_argument(
        "-P", "--placementSize", type=int,
        dest="placement_size", metavar="N", default=None,
        help=("max placement subset size of N "
              "[default: 10%% of the total number of taxa or the alignment "
              "length (whichever bigger)]"))
    decompGroup.add_argument(
        "-F", "--fragmentChunkSize", type=int,
        dest="max_chunk_size", metavar="N", default=20000,
        help=("maximum fragment chunk size of N. Helps controlling memory.  "
              "[default: 20000]"))
    decompGroup.add_argument(
        "-D", "--distance", type=float,
        dest="distance", metavar="DISTANCE",
        default=1,
        help=("minimum p-distance before stopping the decomposition"
              "[default: 1]"))
    # uym2 added #
    decompGroup.add_argument(
        "-M", "--diameter", type=float,
        dest="maxDiam", metavar="DIAMETER",
        default=None,
        help=("maximum tree diameter before stopping the decomposition"
              "[default: None]"))

    decompGroup.add_argument(
        "-S", "--decomp_strategy", type=valid_decomp_strategy,
        dest="decomp_strategy", metavar="DECOMP",
        default="normal",
        # default = "midpoint",
        help="decomposition strategy "
             "[default: using tree branch length]")
    # "[default: only include smallest subsets]")

    outputGroup = _parser.add_argument_group(
        "Output Options".upper(), "These options control output.")
    _parser.groups['outputGroup'] = outputGroup

    outputGroup.add_argument(
        "-p", "--tempdir",
        dest="tempdir", metavar="DIR",
        type=valid_dir_path,
        default=get_default_temp_dir(),
        help=("Tempfile files will be written to DIR. Full-path required. "
              "[default: %(default)s]"))
    outputGroup.add_argument(
        "-rt", "--remtemp",
        dest="remtemp",
        action="store_true",
        help=("Remove tempfile directory.  "
              "[default: disabled]"))
    outputGroup.set_defaults(remtemp=False)
    outputGroup.add_argument(
        "-o", "--output",
        dest="output", metavar="OUTPUT",
        default="output",
        type=valid_file_prefix,
        help="output files with prefix OUTPUT. [default: %(default)s]")
    outputGroup.add_argument(
        "-d", "--outdir",
        dest="outdir", metavar="OUTPUT_DIR",
        default=os.path.curdir,
        type=valid_dir_path,
        help=("output to OUTPUT_DIR directory. full-path required. "
              "[default: %(default)s]"))

    inputGroup = _parser.add_argument_group(
        "Input Options".upper(),
        ' '.join([
            "These options control input. To run SEPP the following is "
            "required. A backbone tree (in newick format), a RAxML_info file "
            "(this is the file generated by RAxML during estimation of the "
            "backbone tree. Pplacer uses this info file to set model "
            "parameters), a backbone alignment file (in fasta format), and a"
            " fasta file including fragments.  The input sequences are assumed"
            " to be DNA unless specified otherwise."]))
    _parser.groups['inputGroup'] = inputGroup

    inputGroup.add_argument(
        "-c", "--config",
        dest="config_file", metavar="CONFIG",
        type=argparse.FileType('r'),
        help=("A config file, including options used to run SEPP. Options"
              " provided as command line arguments overwrite config file "
              "values for those options. "
              "[default: %(default)s]"))
    inputGroup.add_argument(
        "-t", "--tree",
        dest="tree_file", metavar="TREE",
        type=argparse.FileType('r'),
        help="Input tree file (newick format) [default: %(default)s]")
    inputGroup.add_argument(
        "-r", "--raxml",
        dest="info_file", metavar="RAXML",
        type=argparse.FileType('r'),
        help=("RAxML_info file including model parameters, generated by RAxML."
              "[default: %(default)s]"))
    inputGroup.add_argument(
        "-a", "--alignment",
        dest="alignment_file", metavar="ALIGN",
        type=argparse.FileType('r'),
        help="Aligned fasta file [default: %(default)s]")
    inputGroup.add_argument(
        "-f", "--fragment",
        dest="fragment_file", metavar="FRAG",
        type=argparse.FileType('r'),
        help="fragment file [default: %(default)s]")
    inputGroup.add_argument(
        "-m", "--molecule",
        dest="molecule", metavar="MOLECULE",
        type=valid_molecule,
        default="dna",
        help=("Molecule type of sequences. Can be amino, dna, or rna "
              "[default: %(default)s]"))
    inputGroup.add_argument(
        "--ignore-overlap",
        dest="ignore_overlap",
        default=False, action='store_true',
        help=("When a query sequence has the same name as a "
              "backbone sequence, ignore the query "
              "sequences and keep the backbone sequence"
              "[default: %(default)s]"))

    otherGroup = _parser.add_argument_group(
        "Other options".upper(), "These options control how SEPP is run")
    _parser.groups['otherGroup'] = otherGroup
    otherGroup.add_argument(
        "-x", "--cpu", type=set_cpu,
        dest="cpu", metavar="N",
        default=set_cpu(cpu_count()),
        help=("Use N cpus "
              "[default: number of cpus available on the machine]"))
    otherGroup.add_argument(
        "-cp", "--checkpoint", type=set_checkpoint,
        dest="checkpoint", metavar="CHCK_FILE",
        default=set_checkpoint(None),
        help="checkpoint file [default: no checkpointing]")
    otherGroup.add_argument(
        "-cpi", "--interval", type=int,
        dest="checkpoint_interval", metavar="N",
        default=3600,
        help=("Interval (in seconds) between checkpoint writes. Has effect "
              "only with -cp provided. [default: 3600]"))
    otherGroup.add_argument(
        "-seed", "--randomseed", type=int,
        dest="seed", metavar="N",
        default=297834,
        help="random seed number. [default: 297834]")
    # inputGroup.add_argument("-p", "--package",
    #                  dest="package", metavar="PKG",
    #                  help="package directory"
    #                         "[default: %(default)s]")
    #

    return _parser


def get_parser():
    global _parser
    if _parser is None:
        _parser = _init_parser()
    return _parser


def _parse_options():
    parser = get_parser()

    opts = Namespace()

    ''' First read the main configuration file '''
    _LOG.debug("Main configuration file at %s" % main_config_path)
    if not os.path.exists(main_config_path):
        _LOG.warning(
            "Main configuration file was not found at: %s\n" %
            main_config_path + "Proceeding without the main configuration...")
        main_cmd_defaults = []
    else:
        with open(main_config_path, 'r') as cfile:
            main_cmd_defaults = _read_config_file(
                cfile, opts)

    input_args = main_cmd_defaults + (sys.argv[1:])

    ''' Then read the commandline options '''
    opts = parser.parse_args(input_args, namespace=opts)

    ''' If there is a user-specified config file, read that '''
    if opts.config_file is not None:
        config_cmd_defaults = _read_config_file(opts.config_file, opts)

        input_args = main_cmd_defaults + config_cmd_defaults + (sys.argv[1:])

        def error_callback(message):
            newmessage = message.replace(
                "arguments:",
                "arguments (potentially from the config file):").replace(
                "--", "")
            ArgumentParser.error(parser, newmessage)

        parser.error = error_callback

        _LOG.debug(str(input_args))
        ''' Read commandline options again to overwrite config file values'''
        opts = parser.parse_args(input_args, namespace=opts)
    random.seed(opts.seed)
    _LOG.info("Seed number: %d" % opts.seed)
    return opts


_options_singelton = None


def options():
    """
    Returns the configurations read from main configuration file,
    commandline and the user input configuration file.
    """
    global _options_singelton
    if _options_singelton is None:
        _options_singelton = _parse_options()
    return _options_singelton
