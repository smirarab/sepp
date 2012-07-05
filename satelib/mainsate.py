#! /usr/bin/env python

"""Main script of SATe in command-line mode
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
import re
import sys
import signal
import time
import glob
import optparse
import satelib

from satelib import PROGRAM_NAME, PROGRAM_VERSION, PROGRAM_LONG_DESCRIPTION, get_logger, set_timing_log_filepath, TIMING_LOG, MESSENGER, ensure_unique_filename
from satelib.alignment import Alignment, SequenceDataset, MultiLocusDataset, concatenate_alignments, disassemble_alignment
from satelib.configure import get_configuration, get_input_source_directory
from satelib.tree import PhylogeneticTree
from satelib.tools import *
from satelib.sate import *
from satelib.treeholder import read_and_encode_splits
from satelib.sate import SateJob, SateTeam
from satelib.scheduler import start_worker, jobq
from satelib.utility import IndentedHelpFormatterWithNL
from satelib.filemgr import open_with_intermediates
from satelib import filemgr

_RunningJobs = None

_LOG = get_logger(__name__)

def killed_handler(n, frame):
    global _RunningJobs
    if _RunningJobs:
        MESSENGER.send_warning("signal killed_handler called. Killing running jobs...\n")
        j = _RunningJobs
        j.kill()
    else:
        MESSENGER.send_warning("signal killed_handler called with no jobs running. Exiting.\n")
    sys.exit()

def read_input_sequences(seq_filename_list,
        datatype,
        missing=None):
    """
    Return a MultiLocusDataset object or raises an `Exception`

        - `seq_filename_list` should a be a list of paths to FASTA-formatted sequences
        - `datatype` should  be "DNA" or "PROTEIN"
        - `missing` should be "AMBIGUOUS" to "ABSENT" indicate whether these
          missing data symbols should be treated as "any residue" or "absent"

    """
    multilocus_dataset = MultiLocusDataset()
    datatype = datatype.upper()
    if datatype not in ["DNA", "PROTEIN"]:
        raise Exception("Expecting the datatype to by 'DNA' or 'PROTEIN', but found: %s\n" % datatype)

    for seq_fn in seq_filename_list:
        MESSENGER.send_info("Reading input sequences from '%s'..." % seq_fn)
        sd = SequenceDataset()
        try:
            fileobj = open(seq_fn, 'rU')
            sd.read(fileobj, file_format='FASTA', datatype=datatype)
            _LOG.debug("sd.datatype = %s" % sd.datatype)
        except Exception, x:
            raise Exception("Error reading file:\n%s\n" % str(x))

        try:
            if not sd.sequences_are_valid(remap_missing=False):
                m = missing.upper() if missing is not None else "ABSENT"
                if not m in ["AMBIGUOUS", "ABSENT"]:
                    if m:
                        msg = 'The value "%s" for --missing was not understood' % m
                    raise Exception('The missing data symbol ? was encountered.\nExpecting the "missing" command-line option to be either\n "Absent" to delete ? symbols, or\n "Ambiguous" to map them to "any residue".\n%s' % msg)
                assert(sd.alphabet)
                map_missing_to = (m == "AMBIGUOUS" and sd.alphabet.any_residue.symbol or None)
                if not sd.sequences_are_valid(remap_missing=True, map_missing_to=map_missing_to):
                    raise Exception("Input sequences could not be prepared for SATe.  Please report this error\n")
        except Exception, x:
            raise Exception('Error in processing file "%s":\n%s\n' % (seq_fn, str(x)))
        multilocus_dataset.append(sd)
    return multilocus_dataset

def finish_sate_execution(sate_team,
                          user_config,
                          temporaries_dir,
                          multilocus_dataset,
                          sate_products):
    global _RunningJobs


    options = user_config.commandline

    user_config.save_to_filepath(os.path.join(temporaries_dir, 'last_used.cfg'))
    if options.timesfile:
        f = open_with_intermediates(options.timesfile, 'a')
        f.close()
        set_timing_log_filepath(options.timesfile)
    ############################################################################
    # We must read the incoming tree in before we call the get_sequences_for_sate
    #   function that relabels that taxa in the dataset
    ######
    tree_file = options.treefile
    if tree_file:
        if not os.path.exists(tree_file):
            raise Exception('The tree file "%s" does not exist' % tree_file)
        tree_f = open(tree_file, 'rU')
        MESSENGER.send_info('Reading starting trees from "%s"...' % tree_file)
        tree_list = read_and_encode_splits(multilocus_dataset[0].dataset, tree_f)
        tree_f.close()
        if len(tree_list) > 1:
            MESSENGER.send_warning('%d starting trees found in "%s". The first tree will be used.' % (len(tree_list), tree_file))
        starting_tree = tree_list[0]
        score = None

    ############################################################################
    # This will relabel the taxa if they have problematic names
    #####

    alignments = []
    partitions = None
    whole_alignment = None
    for sd in multilocus_dataset:
        alignment = sd.relabel_for_sate()
        alignments.append(alignment)
        _LOG.debug("after relabel_for_sate, alignment.datatype = %s.  sd.datatype = %s" % (alignment.datatype, sd.datatype))
        _LOG.debug("alignment.get_sequence_names() = %s" % str(alignment.get_sequence_names()))

    options.aligned = all( [i.is_aligned() for i in alignments] )

    ############################################################################
    # Launch threads to do work
    #####
    sate_config = user_config.get("sate")
    start_worker(sate_config.num_cpus)

    ############################################################################
    # Be prepared to kill any long running jobs
    #####
    prev_signals = []
    for sig in [signal.SIGTERM, signal.SIGABRT, signal.SIGINT]: # signal.SIGABRT, signal.SIGBUS, signal.SIGINT, signal.SIGKILL, signal.SIGSTOP]:
        prev_handler = signal.signal(sig, killed_handler)
        prev_signals.append((sig, prev_handler))
    
    ######################### Align ####3
    try:
        if not options.aligned:
            MESSENGER.send_info("Performing initial alignment of the entire data matrix...")
            init_aln_dir = os.path.join(temporaries_dir, 'init_aln')
            init_aln_dir = sate_team.temp_fs.create_subdir(init_aln_dir)
            delete_aln_temps = not (options.keeptemp and options.keepalignmenttemps)
            new_alignments= []
            for alignment in alignments:
                job = sate_team.aligner.create_job(alignment,
                                                   tmp_dir_par=init_aln_dir,
                                                   context_str="initalign",
                                                   delete_temps=delete_aln_temps)
                _RunningJobs = job
                jobq.put(job)
                new_alignment = job.get_results()
                _RunningJobs = None
                new_alignments.append(new_alignment)
            alignments = new_alignments
            if delete_aln_temps:
                sate_team.temp_fs.remove_dir(init_aln_dir)
        else:
            MESSENGER.send_info("Input sequences assumed to be aligned (based on sequence lengths).")

        whole_alignment, partitions = concatenate_alignments(alignments)

        if tree_file:
            # getting the newick string here will allow us to get a string that is in terms of the correct taxon labels
            starting_tree_str = starting_tree.compose_newick()
        else:
            MESSENGER.send_info("Performing initial tree search to get starting tree...")


            init_tree_dir = os.path.join(temporaries_dir, 'init_tree')
            init_tree_dir = sate_team.temp_fs.create_subdir(init_tree_dir)
            delete_tree_temps = not options.keeptemp
            job = sate_team.tree_estimator.create_job(whole_alignment,
                                                    tmp_dir_par=init_tree_dir,
                                                    num_cpus=sate_config.num_cpus,
                                                    partitions=partitions,
                                                    context_str="inittree",
                                                    delete_temps=delete_tree_temps)
            _RunningJobs = job
            jobq.put(job)
            score, starting_tree_str = job.get_results()
            _RunningJobs = None
            if delete_tree_temps:
                sate_team.temp_fs.remove_dir(init_tree_dir)

        sate_config_dict = sate_config.dict()

        if options.keeptemp:
            sate_config_dict['keep_iteration_temporaries'] = True
            if options.keepalignmenttemps:
                sate_config_dict['keep_realignment_temporaries'] = True

        job = SateJob(alignments=alignments,
                        sate_team=sate_team,
                        multilocus_dataset=multilocus_dataset,
                        name=options.job,
                        status_messages=MESSENGER.send_info,
                        **sate_config_dict)
        job.partitions = partitions

        job.score = score
        job.best_score = score
        job.tree_str = starting_tree_str

        _RunningJobs = job
        MESSENGER.send_info("Starting SATe algorithm on initial tree...")
        job.run(tmp_dir_par=temporaries_dir)
        _RunningJobs = None

        partitions = job.partitions
        if job.alignment is None:
            new_alignments = disassemble_alignment(whole_alignment, partitions)
        else:
            new_alignments = disassemble_alignment(job.alignment, partitions)


        for p in partitions:
            i = partitions.index(p)
            multilocus_dataset[i].set_alignment(new_alignments[i])
            multilocus_dataset[i].restore_taxon_names()

        if not options.multilocus:
            alignment_stream = sate_products.alignment_streams[0]
            MESSENGER.send_info("Writing final alignment to %s" % alignment_stream.name)
            sd.dataset.write(alignment_stream, schema="FASTA")
        else:
            ### TODO : URGENT !!!!! ###
            raise NotImplementedError("Still need to figure out multilocus output")
            #aln_dirname = options.output
            #if aln_dirname is None:
            #    aln_dirname = seqdir
            #seqfns = glob.glob(os.path.join(os.path.abspath(seqdir), '*.fas'))
            #seqfns.extend(glob.glob( os.path.join(os.path.abspath(seqdir), '*.fasta')))

            #for seqfn in seqfns:
            #    i = seqfns.index(seqfn)
            #    aln_filename = os.path.join(aln_dirname, os.path.basename(seqfn)+"_%s.aln" % options.job)
            #    MESSENGER.send_info("Writing final alignment to %s" % aln_filename)
            #    aln_file = open_with_intermediates( aln_filename, 'w')
            #    sd.dataset.write(aln_file, schema="FASTA")
            #    aln_file.close()

        MESSENGER.send_info("Writing final tree to %s" % sate_products.tree_stream.name)
        tree_str = job.tree.compose_newick()
        sate_products.tree_stream.write("%s;\n" % tree_str)

        MESSENGER.send_info("Writing final likelihood score to %s" % sate_products.score_stream.name)
        sate_products.score_stream.write("%s\n" % job.score)
    finally:
        for el in prev_signals:
            sig, prev_handler = el
            if prev_handler is None:
                signal.signal(sig, signal.SIG_DFL)
            else:
                signal.signal(sig, prev_handler)

def run_sate_from_config(seq_filename_list,
        user_config,
        sate_products):
    """
    Returns `None` if no temporary directory is left over from the run
    or return the path to the temporary directory created for the
    scratch files.
    """

    multilocus_dataset = read_input_sequences(seq_filename_list,
            datatype=user_config.commandline.datatype,
            missing=user_config.commandline.missing)

    cmdline_options = user_config.commandline

    ############################################################################
    # Create the safe directory for temporaries
    # The general form of the directory is
    #   ${options.temporaries}/${options.job}/temp${RANDOM}
    ######
    par_dir = cmdline_options.temporaries
    if par_dir is None:
        par_dir = os.path.join(os.path.expanduser('~'), '.sate')
    cmdline_options.job = coerce_string_to_nice_outfilename(cmdline_options.job, "Job", "satejob")
    subdir = cmdline_options.job
    par_dir = os.path.abspath(os.path.join(par_dir, subdir))
    if not os.path.exists(par_dir):
        os.makedirs(par_dir) # this parent directory will not be deleted, so we don't store it in the sate_team.temp_fs

    sate_team = SateTeam(config=user_config)

    delete_dir = not cmdline_options.keeptemp

    temporaries_dir = sate_team.temp_fs.create_top_level_temp(parent=par_dir, prefix='temp')
    assert(os.path.exists(temporaries_dir))
    try:
        MESSENGER.send_info("Directory for temporary files created at %s" % temporaries_dir)
        finish_sate_execution(sate_team=sate_team,
                              user_config=user_config,
                              temporaries_dir=temporaries_dir,
                              multilocus_dataset=multilocus_dataset,
                              sate_products=sate_products)
    finally:
        if delete_dir:
            sate_team.temp_fs.remove_dir(temporaries_dir)
    if delete_dir:
        return None
    else:
        return temporaries_dir

def coerce_string_to_nice_outfilename(p, reason, default):
    illegal_filename_pattern = re.compile(r'[^-_a-zA-Z0-9.]')
    j = "".join(illegal_filename_pattern.split(p))
    if not j:
        j = default
    if j != p:
        MESSENGER.send_warning('%s name changed from "%s" to "%s" (a safer name for filepath)' % (reason, p, j))
    return j

def main(argv=sys.argv):
    '''Returns (`True`, dir) on successful execution or raises an exception.

    Where `dir` is either None or the undeleted directory of temporary files.

    Note that if `argv` is sys.argv then the first element will be skipped, but
        if it is not the sys.argv list then the first element will be interpretted
        as an argument (and will NOT be skipped).
    '''

    _START_TIME = time.time()
    usage = """usage: %prog [options] <settings_file1> <settings_file2> ..."""
    parser = optparse.OptionParser(usage=usage,
                                    description=PROGRAM_LONG_DESCRIPTION,
                                    formatter=IndentedHelpFormatterWithNL(),
                                    version="%s v%s" % (PROGRAM_NAME, PROGRAM_VERSION))

    user_config = get_configuration()
    command_line_group = user_config.get('commandline')
    command_line_group.add_to_optparser(parser)
    sate_group = user_config.get('sate')
    sate_group.add_to_optparser(parser)
    if argv == sys.argv:
        (options, args) = parser.parse_args(argv[1:])
    else:
        (options, args) = parser.parse_args(argv)
    if options.multilocus:
        sys.exit("SATe: Multilocus mode is disabled in this release.")
    config_filenames = list(args)
    for fn in config_filenames:
        if fn[0] == '"' and fn[-1] == '"':
            fn = fn[1:-1]
        if not os.path.exists(fn):
            raise Exception('The configuration (settings) file "%s" does not exist' % fn)
        try:
            user_config.read_config_filepath(fn)
        except:
            raise Exception('The file "%s" does not appear to be a valid configuration file format. It lacks section headers.' % fn)
    user_config.set_values_from_dict(options.__dict__)
    command_line_group.job = coerce_string_to_nice_outfilename(command_line_group.job, 'Job', 'satejob')

    exportconfig = command_line_group.exportconfig
    if exportconfig:
        command_line_group.exportconfig = None
        user_config.save_to_filepath(exportconfig)

        ### TODO: wrap up in messaging system
        sys.stdout.write('Configuration written to "%s". Exiting successfully.' % exportconfig )

        return True, None

    if user_config.commandline.input is None:
        sys.exit("ERROR: Input file(s) not specified.")

    # note: need to read sequence files first to allow SateProducts to
    # correctly self-configure
    seq_filename_list = user_config.read_seq_filepaths(src=user_config.commandline.input,
            multilocus=user_config.commandline.multilocus)
    sate_products = filemgr.SateProducts(user_config)
    #sate_products.setup()

    MESSENGER.run_log_streams.append(sate_products.run_log_stream)
    MESSENGER.err_log_streams.append(sate_products.err_log_stream)
    temp_dir = run_sate_from_config(seq_filename_list, user_config, sate_products)
    _TIME_SPENT = time.time() - _START_TIME
    MESSENGER.send_info("Total time spent: %ss" % _TIME_SPENT)
    return True, temp_dir

