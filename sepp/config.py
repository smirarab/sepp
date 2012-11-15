'''
This module provides options(), which returns a argparse.Namespace object
containing all configurations required by sepp. config.options() reads and saves 
the configurations the first time it is called. Subsequent calls simply return
the already-saved configurations. A Typical usage is:
some_config_attribute = sepp.config.options().some_config_attribute  

All command-line options are found directly inside the return value of options()
as attributes. Commandline options can be specified either directly using the 
command-line, or under [commandline] section in a config file, passed to SEPP 
using -c commandline option. commandline values overwrite config file values. 
 
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
'''                

from argparse import ArgumentParser, Namespace
from sepp.filemgr import get_default_temp_dir, check_or_make_dir_path
import sys
import ConfigParser
from sepp import version, get_logger
import argparse
import os
from sepp.scheduler import JobPool
from multiprocessing import cpu_count

_LOG = get_logger(__name__)
   
main_config_path = os.path.expanduser("~/.sepp/main.config")   

def _read_config_file(filename, opts):
    config_defaults = []
    configparser = ConfigParser.ConfigParser()    
    configparser.optionxform = str
    configparser.readfp(filename)

    if configparser.has_section('commandline'):
        for (k,v) in configparser.items('commandline'):    
            config_defaults.append("--%s" %k)
            config_defaults.append(v)
    
    for section in configparser.sections():
        if section == "commandline":
            continue
        section_name_space = Namespace()
        for (k,v) in configparser.items(section):    
            section_name_space.__setattr__(k,v)        
        opts.__setattr__(section, section_name_space)
    
    return config_defaults

def valid_dir_path(path):
    ret = check_or_make_dir_path(path)
    if ret is None:
        raise argparse.ArgumentTypeError("%s is not a valid directory path." %path)
    return ret

def valid_file_prefix(prefix):
    if os.path.dirname(prefix) != "":
        raise argparse.ArgumentTypeError("%s is not a valid output prefix (includes a directory)." %prefix)
    return prefix

def set_cpu(cpus):
    c = int(cpus)
    JobPool(c)
    return c
    
def _get_parser():
    parser = ArgumentParser(description= 
                            "This script runs the SEPP algorithm on an input "
                            "tree, alignment, fragment file, and RAxML info file.")    
    
    parser.add_argument("-v", "--version", action='version', version= "%(prog)s " + version)

    decompGroup = parser.add_argument_group("Decomposition Options".upper(), 
                         ' '.join(["These options determine the alignment decomposition size and", 
                                 "taxon insertion size.  If None is given, then the default",
                                 "is to align/place at 10% of total taxa.  The alignment decomosition size must be",
                                 "less than the taxon insertion size."]))                                 
    
    decompGroup.add_argument("-A", "--alignmentSize", type = int, 
                      dest = "alignment_size", metavar = "N", 
                      default = None,
                      help = "max alignment subset size of N "
                             "[default: 10%% of the total number of taxa]")    
    decompGroup.add_argument("-P", "--placementSize", type = int, 
                      dest = "placement_size", metavar = "N",
                      default = None, 
                      help = "max placement subset size of N "
                             "[default: 10%% of the total number of taxa]")                              
        
    outputGroup = parser.add_argument_group( "Output Options".upper(), 
                         "These options control output.")                                     
    outputGroup.add_argument("-p", "--tempdir", 
                      dest = "tempdir", metavar = "DIR",
                      type=valid_dir_path,
                      default = get_default_temp_dir(),                       
                      help = "Tempfile files will be written to DIR. Full-path required. "
                             "[default: %(default)s]")    
    outputGroup.add_argument("-o", "--output", 
                      dest = "output", metavar = "OUTPUT",
                      default = "output", 
                      type= valid_file_prefix,
                      help = "output files with prefix OUTPUT. "
                             "[default: %(default)s]")
    outputGroup.add_argument("-d", "--outdir", 
                      dest = "outdir", metavar = "OUTPUT_DIR", 
                      default = os.path.curdir, 
                      type = valid_dir_path,
                      help = "output to OUTPUT_DIR directory. full-path required. "
                             "[default: %(default)s]")                       
                             
    inputGroup = parser.add_argument_group ("Input Options".upper(), 
                         ' '.join(["These options control input. To run SEPP the following is required." 
                                 "A backbone tree (in newick format), a RAxML_info file (this is the file generated by RAxML during estimation of the backbone tree. " 
                                 "Pplacer uses this info file to set model parameters),"
                                 "a backbone alignment file (in fasta format), and a fasta file including fragments."]))                                     
    inputGroup.add_argument("-c", "--configparser", 
                      dest = "config_file", metavar = "CONFIG",
                      type = argparse.FileType('r'), 
                      help = "A configparser file, including options used to run SEPP. Options provided as command line arguments overwrite configparser file values for those options. "
                             "[default: %(default)s]")    
    inputGroup.add_argument("-t", "--tree", 
                      dest = "tree_file", metavar = "TREE",
                      type = argparse.FileType('r'), 
                      help = "Input tree file (newick format) "
                             "[default: %(default)s]")    
    inputGroup.add_argument("-r", "--raxml", 
                      dest = "info_file", metavar = "RAXML",
                      type = argparse.FileType('r'), 
                      help = "RAxML_info file including model parameters, generated by RAxML."
                             "[default: %(default)s]")    
    inputGroup.add_argument("-a", "--alignment", 
                      dest = "alignment_file", metavar = "ALIGN",
                      type = argparse.FileType('r'), 
                      help = "Aligned fasta file "
                             "[default: %(default)s]")    
    inputGroup.add_argument("-f", "--fragment",
                      dest = "fragment_file", metavar = "FRAG",
                      type = argparse.FileType('r'), 
                      help = "fragment file "
                             "[default: %(default)s]")          
        
    otherGroup = parser.add_argument_group( "Other options".upper(), 
                         "These options control how SEPP is run")                                     
    otherGroup.add_argument("-x", "--cpu", type = set_cpu, 
                      dest = "cpu", metavar = "N", 
                      default = cpu_count(),
                      help = "Use N cpus "
                             "[default: number of cpus available on the machine]")
                                                          
    #inputGroup.add_argument("-p", "--package", 
    #                  dest = "package", metavar = "PKG", 
    #                  help = "package directory"
    #                         "[default: %(default)s]")                                                          
    #                         

    
    return parser

def _parse_options ():
    
    parser = _get_parser()                                 
        
    opts = Namespace()
    
    ''' First read the main configuration file '''
    if not os.path.exists(main_config_path):
        _LOG.warn("Main configuration file was not found at: %s\n" %main_config_path + 
                  "Proceeding without the main configuration..." )
        main_cmd_defaults = []
    else:
        main_cmd_defaults = _read_config_file(open(main_config_path,'r'),opts)
    
    input_args = main_cmd_defaults + (sys.argv[1:])        
    
    ''' Then read the commandline options '''
    opts = parser.parse_args(input_args, namespace = opts)
    
    ''' If there is a user-specified config file, read that '''    
    if opts.config_file is not None:
        
        config_cmd_defaults = _read_config_file(opts.config_file,opts)
              
        input_args = main_cmd_defaults + config_cmd_defaults + (sys.argv[1:])
        
        def error_callback(message):
            newmessage = message.replace("arguments:",
                            "arguments (potentially from the config file):"
                            ).replace("--","")
            ArgumentParser.error(parser, newmessage)
            
        parser.error = error_callback 
        
        ''' Read commandline options again to overwrite config file values'''    
        opts = parser.parse_args(input_args, namespace = opts )    

    return opts

_options_singelton = None        
     
def options():
    '''
    Returns the configurations read from main configuration file, 
    commandline and the user input configuration file.    
    '''
    global _options_singelton
    if _options_singelton is None:
        _options_singelton = _parse_options()        
    return _options_singelton    
                                    