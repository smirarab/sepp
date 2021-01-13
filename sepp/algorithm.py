"""
Created on Oct 2, 2012

@author: smirarab
"""
from sepp.config import options
from abc import abstractmethod
from sepp.scheduler import JobPool
from sepp.filemgr import directory_has_files_with_prefix, get_temp_file
from sepp.filemgr import remove_temp, get_root_temp_dir
from sepp.alignment import MutableAlignment
from sepp.tree import PhylogeneticTree
import dendropy
from sepp import get_logger
import sys
import os
from sepp.problem import RootProblem
import time
from sepp.checkpointing import CheckPointManager

_LOG = get_logger(__name__)


class AbstractAlgorithm(object):
    """
    This class provides the interface for an abstract algorithm. All different
    ways of implementing SEPP are considered to be different 'algorithms'.

    New algorithms should be implemented by extending AbastractAlgorithm and
    implementing its abstract methods.
    """

    def __init__(self):
        """
        Constructor
        """
        self.root_problem = None
        self.results = None
        self.options = options()
        self.outchecked = False  # for ease of access

    def check_options(self, supply=[]):
        """
        This method should check the input values stored in config.option to
        make sure every necessary argument is provided, and that the provided
        values are all fine.

        In the event of recognizing invalid or missing input, an Exception
        (maybe an ArgumentError) should be raised.

        By default expects tree_file, raxml_file, and fragment_file. Overwrite
        if required.
        """
        if options().tree_file is None:
            supply = supply + ["tree file"]
        if options().alignment_file is None:
            supply = supply + ["alignment file"]
        if options().fragment_file is None:
            supply = supply + ["fragment file"]
        if len(supply) != 0:
            raise ValueError(
                ("Failed to supply: %s\nRun with -h option to see a list of"
                 " options.")
                % " , ".join(supply))
        self.check_outputprefix()

    @abstractmethod
    def build_subproblems(self):
        """
        This method should read the config.options() and build a hierarchy of
        subproblems (see sepp.problem). This hierarchy will be used by
        build_job_dag to create a DAG of jobs (see sepp.jobs and
        sepp.scheduler)
        """
        raise NotImplementedError()

    @abstractmethod
    def connect_jobs(self):
        """
        Join Jobs are joined together to form a DAG using Joins (see
        sepp.scheduler)
        and call_back functions (see sepp.scheduler.Job).

        Once the first level of jobs (those with no dependency) are enqueued
        (using enqueue_firstlevel_job) everything else should be automatically
        enqueued when all their dependencies are satisfied.
        For example, if the conceptual DAG is:
        A----> E---|--> F
                   |
        B--|-> D---|
        C--|
        then enqueue_firstlevel_job should enqueue A B and C. When A finishes,
        it should enqueue E (either using a callback or a Join with only one
        dependency). B and C should be joined together using a Join, and that
        Join needs to enqueue D. Also E and D need to be joined, and their join
        needs to enqueue F.
        build_jobs is called before this, and all jobs are saved as part of the
        subproblem hierarchy. This function only connects those jobs.
        This is called after recovery from a checkpoint, because joins and
        callbacks are not checkpointed.
        """
        raise NotImplementedError()

    @abstractmethod
    def build_jobs(self):
        """
        Builds separate jobs for different tasks that need to be done.
        This just creates jobs, without connecting them. connect_jobs is used
        to create jobs.
        build_job is not called after checkpoints are recovered, because the
        jobs are checkpointed (their connections are not).
        """
        raise NotImplementedError()

    @abstractmethod
    def enqueue_firstlevel_job(self):
        """
        This is called after the DAG is created (see buld_job_dag) to enqueue
        the first level of jobs (those with no dependency). These jobs should
        automatically enqueue the rest of the DAG upon completion.
        """
        raise NotImplementedError()

    @abstractmethod
    def merge_results(self):
        """
        This method should read results and prepare all the final results
        for output. Resutls should be saved in the problem hierarchy, and
        by setting ?? attribute
        """
        raise NotImplementedError()

    @abstractmethod
    def output_results(self):
        """
        This method should output the final results. Final results can be found
        in the problem hierarchy
        """
        raise NotImplementedError()

    def run(self):
        checkpoint_manager = options().checkpoint
        assert isinstance(checkpoint_manager, CheckPointManager)

        t = time.time()

        if checkpoint_manager.is_recovering:
            checkpoint_manager.restore_checkpoint()
            self.root_problem = \
                checkpoint_manager.checkpoint_state.root_problem
            self.check_outputprefix()
        else:
            '''check input arguments'''
            self.check_options()

            '''build the problem structure'''
            self.root_problem = self.build_subproblems()

            '''build jobs'''
            self.build_jobs()

        '''connect jobs into a DAG'''
        self.connect_jobs()

        '''Queue up first level jobs (i.e. those with no dependency).
        Once these run, they should automatically enqueue the rest of the
        DAG through joins and callbacks '''
        self.enqueue_firstlevel_job()

        '''start the checkpointing (has any effects only in
           checkpointing mode)'''
        checkpoint_manager.start_checkpointing(self.root_problem)

        '''Wait for all jobs to finish'''
        if not JobPool().wait_for_all_jobs():
            _LOG.exception(
                "There have been errors in executed jobs. Terminating.")
            sys.exit(1)

        ''' terminate The job pool and release memory'''
        JobPool().terminate()

        ''' Pause Checkpointing'''
        checkpoint_manager.pause_checkpointing()
        # checkpoint_manager.force_checkpoint()

        '''Merge results into final outputs'''
        self.merge_results()

        '''Output final results'''
        self.output_results()

        ''' Pause Checkpointing'''
        checkpoint_manager.stop_checkpointing()

        _LOG.info("Current execution Finished in %d seconds"
                  % (time.time() - t))
        _LOG.info(
            "All checkpointed executions Finished in %d cumulative time" %
            (checkpoint_manager.get_total_time()))

        if self.options.remtemp:
            _LOG.info("Now removing temp files from %s" % str(
                get_root_temp_dir()))
            remove_temp(get_root_temp_dir())
        else:
            _LOG.info("Temp files are left behind at %s" % str(
                get_root_temp_dir()))

    ''' The following are a bunch of helper methods that will be needed in
    most implementations of sepp'''

    def check_and_set_sizes(self, total):
        # If sizes are not set, then use 10% rule
        my_options = self.options
        if (my_options.maxDiam is not None) and \
           (my_options.decomp_strategy not in ["midpoint", "centroid"]):
            raise Exception(("The max diameter option can be used only with "
                             "the midpoint or the centroid decomposition "
                             "specified using -S"))
        if my_options.alignment_size is None:
            if my_options.placement_size is None:
                my_options.alignment_size = int(total * .10)
            else:
                my_options.alignment_size = my_options.placement_size
        if my_options.placement_size is None:
            my_options.placement_size = max(int(total * .10),
                                            my_options.alignment_size)
        if my_options.placement_size is not None and \
           my_options.placement_size < my_options.alignment_size:
            _LOG.warning(
                ("alignment_size (%d) cannot be larger than placement_size "
                 "(%d).   Setting alignment_size to be placement_size") %
                (my_options.alignment_size, my_options.placement_size))
            my_options.alignment_size = my_options.placement_size
        if my_options.placement_size > total:
            my_options.placement_size = total
            _LOG.warning(
                ("Input placement size larger than total number of sequences!"
                 "  Setting placement size to %d") % total)
        if my_options.alignment_size > total:
            my_options.alignment_size = total
            _LOG.warning(
                ("Input alignment size larger than total number of sequences!"
                 "  Setting alignment size to %d") % total)

        _LOG.info(
            "Decomposition Sizes are set to alignment: %d placement: %d" %
            (my_options.alignment_size, my_options.placement_size))

    def check_outputprefix(self):
        if self.outchecked:
            return
        self.outchecked = True
        if directory_has_files_with_prefix(
           self.options.outdir, self.options.output):
            raise ValueError(
                ("Output directory [%s] already contains files with prefix "
                 "[%s]...\nTerminating to avoid loss of existing files.") % (
                    self.options.outdir, self.options.output))

    def get_output_filename(self, name):
        return os.path.join(
            self.options.outdir, "%s_%s" % (self.options.output, name))

    def read_alignment_and_tree(self):
        _LOG.info("Reading input alignment: %s" % (
            self.options.alignment_file))
        alignment = MutableAlignment()
        alignment.read_file_object(self.options.alignment_file)

        # fragments = MutableAlignment()
        # fragments.read_file_object(self.options.fragment_file);
        _LOG.info("Reading input tree: %s" % self.options.tree_file)
        tree = PhylogeneticTree(
            dendropy.Tree.get_from_stream(self.options.tree_file,
                                          schema="newick",
                                          preserve_underscores=True))

        return alignment, tree

    def read_and_divide_fragments(self, chunks, extra_frags={}):
        max_chunk_size = self.options.max_chunk_size
        _LOG.debug(
            ("start reading fragment files and breaking to at least %s chunks"
             " but at most %s sequences ") % (
                str(chunks), str(max_chunk_size)))
        self.root_problem.fragments = MutableAlignment()
        self.root_problem.fragments.read_file_object(
            self.options.fragment_file)

        # test if input fragment names might collide with reference names.
        # code contribution by Stefan Janssen (June 13th, 2018)
        ids_reference = set(self.root_problem.subalignment.keys())
        ids_inputfragments = set(self.root_problem.fragments.keys())
        ids_overlap = ids_reference & ids_inputfragments
        if len(ids_overlap) > 0 and not self.options.ignore_overlap:
            raise ValueError((
                "Your input fragment file contains %i sequences, whose names "
                "overlap with names in your reference. Please rename your inp"
                "ut fragments and re-start. Duplicate names are:\n  '%s'") %
                (len(ids_overlap), "'\n  '".join(ids_overlap)))
        elif len(ids_overlap) > 0:
            _LOG.debug("Ignoring following %i query sequences present "
                       "in the backbone: \n '%s'"
                       % (len(ids_overlap), "' , '".join(ids_overlap)))

            self.root_problem.fragments = self.root_problem.fragments.\
                get_soft_sub_alignment(ids_inputfragments - ids_reference)

        # test if input fragment names contain whitespaces / tabs which would
        # cause hmmsearch to fail.
        # code contribution by Stefan Janssen (June 22nd, 2018)
        ids_inputfragments_spaces = [
            id_
            for id_ in ids_inputfragments
            if (' ' in id_) or ('\t' in id_)]
        if len(ids_inputfragments_spaces) > 0:
            raise ValueError((
                "Your input fragment file contains %i sequences, whose names "
                "contain either whitespaces: ' ' or tabulator '\\t' symbols. "
                "Please rename your input fragments and re-start. Affected "
                "names are:\n  '%s'") %
                (len(ids_inputfragments_spaces),
                 "'\n  '".join(ids_inputfragments_spaces)))

        for (k, v) in extra_frags.items():
            self.root_problem.fragments[k] = v.replace("-", "")
        alg_chunks = self.root_problem.fragments.divide_to_equal_chunks(
            chunks, max_chunk_size)
        ret = []
        for i in range(0, len(alg_chunks)):
            temp_file = None
            if alg_chunks[i]:
                temp_file = get_temp_file(
                    "fragment_chunk_%d" % i, "fragment_chunks", ".fasta")
                alg_chunks[i].write_to_path(temp_file)
            ret.append(temp_file)
        _LOG.debug("fragment files read and divided.")
        return ret

    def _create_root_problem(self, tree, alignment):
        """ Create the root problem"""
        self.root_problem = RootProblem(tree.leaf_node_names())
        self.root_problem.label = "root"
        self.root_problem.subalignment = alignment
        self.root_problem.subtree = tree
