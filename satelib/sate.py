#!/usr/bin/env python

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


import math
import os
import random
import time
import sys
import copy
from threading import Lock
from satelib import get_logger
_LOG = get_logger(__name__)

from satelib.treeholder import TreeHolder
from satelib.satealignerjob import SateAlignerJob
from satelib.alignment import concatenate_alignments
from satelib import get_logger
from satelib.utility import record_timestamp
from satelib.scheduler import jobq
from satelib.filemgr import  TempFS

class SateTeam (object):
    '''A blob for holding the appropriate merger, alignment, and tree_estimator tools
    as well as the TempFS object that keeps track of the directories that have
    been created by this process.

    '''
    def __init__(self, config):
        """Uses a configuration object `cfg` to get reference for the tools the user
        has chosen.
        """
        try:
            self._temp_fs = TempFS()
            self.aligner = config.create_aligner(temp_fs=self._temp_fs)
            self.merger = config.create_merger(temp_fs=self._temp_fs)
            self.tree_estimator = config.create_tree_estimator(temp_fs=self._temp_fs)
        except AttributeError:
            raise
            raise ValueError("config cannot be None unless all of the tools are passed in.")
    def get_temp_fs(self):
        return self._temp_fs
    temp_fs = property(get_temp_fs)

class AcceptMode:
    BLIND_MODE, NONBLIND_MODE = range(2)

class SateJob (TreeHolder):
    """The top-level SATe algorithm class.  The run_sate method iteratively
    decomposes the tree and does tree searching until the termination criterion
    is reached.
    """
    BEHAVIOUR_DEFAULTS = {  'time_limit' : 24*60*60.0,
                            'iter_limit' : -1,
                            'time_without_imp_limit' : -1,
                            'iter_without_imp_limit' : -1,
                            'blind_after_total_time' : -1,
                            'blind_after_total_iter' : -1,
                            'blind_after_time_without_imp' : -1,
                            'blind_after_iter_without_imp' : -1,
                            'after_blind_time_term_limit' : 24*60*60.0,
                            'after_blind_iter_term_limit' : -1,
                            'after_blind_time_without_imp_limit': -1,
                            'after_blind_iter_without_imp_limit': 1,
                            'move_to_blind_on_worse_score' : False,
                            'blind_mode_is_final' : True,
                            'max_subproblem_size' : 200,
                            'max_subproblem_frac' : 0.2,
                            'start_tree_search_from_current' : False,
                            'keep_realignment_temporaries' : False,
                            'keep_iteration_temporaries' : False,
                            'return_final_tree_and_alignment' : False
                        }
    def configuration(self):
        d = {}
        for k in SateAlignerJob.BEHAVIOUR_DEFAULTS.keys():
            d[k] = getattr(self, k)
        for k in SateJob.BEHAVIOUR_DEFAULTS.keys():
            d[k] = getattr(self, k)
        return d

    def __init__(self, alignments, sate_team, multilocus_dataset, tree=None, name=None, **kwargs):
        TreeHolder.__init__(self, multilocus_dataset[0].dataset)
        self._job_lock = Lock()
        self.alignments = alignments
        self.sate_team = sate_team
        self.datasets = [sd.dataset for sd in multilocus_dataset]
        #self.model = kwargs.get("model")
        self.tree = tree
        self.score = None
        self.blind_mode_is_final = True
        self.is_stuck_in_blind = False

        self._tree_build_job = None
        self._sate_decomp_job = None
        self._reset_jobs()

        self._status_message_func = kwargs.get('status_messages')

        behavior = copy.copy(SateJob.BEHAVIOUR_DEFAULTS)
        behavior.update(SateAlignerJob.BEHAVIOUR_DEFAULTS)
        behavior.update(kwargs)
        self.__dict__.update(behavior)

        # right now max_subproblem_size can be an integer or proportion.
        # These two run modes should be separate parameters!

        # if self.max_subproblem_size <= 1.0:
        #   self.max_subproblem_size *= alignment.get_num_taxa()
        # self.max_subproblem_size = max(1, int(round(self.max_subproblem_size)))

        if isinstance(self.break_strategy, str):
            self.break_strategy = [self.break_strategy]

        self.start_time = None
        self.current_iteration = 0
        self.last_improvement_iter = 0
        self.last_improvement_time = None
        self.best_alignment = None
        self.best_tree_str = self.get_tree_str()
        self.best_score = self.score
        self.is_stuck_in_blind = False
        self.switch_to_blind_iter = None
        self.switch_to_blind_timestamp = None
        self._termination_trigger = None
        self._blindmode_trigger = None
        self._reset_current_run_settings()
        self.killed = False
        self.partitions = None

    def _reset_current_run_settings(self):
        self.start_time = None
        self.current_iteration = 0
        self.best_alignment = None
        self.best_tree_str = self.get_tree_str()
        self.best_score = self.score
        self.num_iter_since_imp = 0
        self.is_stuck_in_blind = False
        self.switch_to_blind_iter = None
        self.switch_to_blind_timestamp = None
        self._termination_trigger = None
        self._blindmode_trigger = None
        self._sate_alignment_job = None
        self._tree_build_job = None

    def _reset_jobs(self):
        self.tree_build_job = None
        self.sate_decomp_job = None

    # We lock the job objects because kill might be called from another thread
    def get_sate_alignment_job(self):
        return self._sate_alignment_job

    def set_sate_alignment_job(self, val):
        self._job_lock.acquire()
        self._sate_alignment_job = val
        self._job_lock.release()
    sate_aligner_job = property(get_sate_alignment_job, set_sate_alignment_job)

    def get_tree_build_job(self):
        return self._tree_build_job

    def set_tree_build_job(self, val):
        self._job_lock.acquire()
        self._tree_build_job = val
        self._job_lock.release()
    tree_build_job = property(get_tree_build_job, set_tree_build_job)

    def _curr_running_times(self):
        "Duration of the current run (in seconds)"
        if self.start_time is None:
            return 0.0
        curr_time = time.time()
        return curr_time - self.start_time

    def _curr_time_since_imp(self):
        "Duration of the current run (in seconds)"
        if self.start_time is None:
            return 0.0
        curr_time = time.time()
        return curr_time - self.last_improvement_time

    def _keep_iterating(self):

        # SATe has run for at least 'iter_limit' iterations
        if self.iter_limit >= 0 and self.current_iteration >= self.iter_limit:
            self._termination_trigger = 'iter_limit'
            return False

        # SATe has run for 'iter_without_imp_limit' iterations and the likelihood score does not improve
        if self.iter_without_imp_limit >= 0 and self.num_iter_since_imp >= self.iter_without_imp_limit:
            self._termination_trigger = 'iter_without_imp_limit'
            return False

        # SATe is in BlIND mode and it has run 'after_blind_iter_term_limit' iterations since then
        if (self.switch_to_blind_iter is not None) and self.after_blind_iter_term_limit >= 0:
            if (self.current_iteration - self.switch_to_blind_iter) >= self.after_blind_iter_term_limit:
                self._termination_trigger = 'after_blind_iter_term_limit'
                return False

        # SATe has run 'time_limit' seconds
        if self.time_limit >= 0.0:
            running_time = self._curr_running_times()
            if running_time > self.time_limit:
                self._termination_trigger = 'time_limit'
                return False

        # SATe has run 'time_without_imp_limit' seconds and the likelihood score does not improve
        if self.time_without_imp_limit >= 0.0:
            time_since_last_imp = self._curr_time_since_imp()
            if time_since_last_imp > self.time_without_imp_limit:
                self._termination_trigger = 'time_without_imp_limit'
                return False

        # SATe is in BlIND mode and it has run 'after_blind_time_term_limit' seconds since then
        if (self.switch_to_blind_timestamp is not None) and self.after_blind_time_term_limit >= 0.0:
            time_since_switching_to_blind = time.time() - self.switch_to_blind_timestamp
            if time_since_switching_to_blind > self.after_blind_time_term_limit:
                self._termination_trigger = 'after_blind_time_term_limit'
                return False

        # SATe is in BlIND mode and it has run 'self.after_blind_time_without_imp_limit' seconds since then
        # without improvements in likelihood score
        if (self.switch_to_blind_timestamp is not None) and self.after_blind_time_without_imp_limit >= 0.0:
            time_since_last_imp = self._curr_time_since_imp()
            time_since_switching_to_blind = time.time() - self.switch_to_blind_timestamp

            if time_since_switching_to_blind > self.after_blind_time_without_imp_limit:
                self._termination_trigger = 'after_blind_time_without_imp_limit'
                return False

        # SATe is in BlIND mode and it has run 'self.after_blind_iter_without_imp_limit' iterations since then
        # without improvements in likelihood score
        if (self.switch_to_blind_timestamp is not None) and self.after_blind_iter_without_imp_limit >= 0:
            if self.num_iter_since_imp >= self.after_blind_iter_without_imp_limit:
                self._termination_trigger = 'after_blind_iter_without_imp_limit'
                return False

        return True

    def _get_break_strategy(self, break_strategy_index):
        try:
            return self.break_strategy[break_strategy_index]
        except:
            return None

    def _get_accept_mode(self, new_score, break_strategy_index=0):
        # If we have several break_strategies,  we do not go to "BLIND_MODE" until we have tried them all out
        if break_strategy_index < (len(self.break_strategy) - 1):
            return AcceptMode.NONBLIND_MODE
        if self.is_stuck_in_blind:
            return AcceptMode.BLIND_MODE
        if self.move_to_blind_on_worse_score and (new_score <= self.score):
            self._blindmode_trigger = 'move_to_blind_on_worse_score'
            return AcceptMode.BLIND_MODE
        if (self.blind_after_total_iter >= 0) and (self.current_iteration >= self.blind_after_total_iter):
            self._blindmode_trigger = 'blind_after_total_iter'
            return AcceptMode.BLIND_MODE
        if (self.blind_after_iter_without_imp >= 0) and (self.num_iter_since_imp >= self.blind_after_iter_without_imp):
            self._blindmode_trigger = 'blind_after_iter_without_imp'
            return AcceptMode.BLIND_MODE
        if self.blind_after_total_time >= 0.0:
            running_time = self._curr_running_times()
            if running_time > self.blind_after_total_time:
                self._blindmode_trigger = 'blind_after_total_time'
                return AcceptMode.BLIND_MODE
        if self.blind_after_time_without_imp >= 0.0:
            time_since_last_imp = self._curr_time_since_imp()
            if time_since_last_imp > self.blind_after_time_without_imp:
                self._blindmode_trigger = 'blind_after_time_without_imp'
                _LOG.debug("blind_after_time_without_imp reached")

        return AcceptMode.NONBLIND_MODE

    def run(self, tmp_dir_par):
        assert(os.path.exists(tmp_dir_par))

        self._reset_current_run_settings()
        self._reset_jobs()

        self.start_time = time.time()
        self.last_improvement_time = self.start_time

        num_non_update_iter = 0

        configuration = self.configuration()
        # Here we check if the max_subproblem_frac is more stringent than max_subproblem_size
        frac_max = int(math.ceil(self.max_subproblem_frac*self.tree.n_leaves))
        if frac_max > self.max_subproblem_size:
            configuration['max_subproblem_size'] = frac_max

        delete_iteration_temps = not self.keep_iteration_temporaries
        delete_realignment_temps = delete_iteration_temps or (not self.keep_realignment_temporaries)
        configuration['delete_temps'] = delete_realignment_temps

        while self._keep_iterating():
            record_timestamp(os.path.join(tmp_dir_par, 'start_sateiter_timestamp.txt'))

            # create a subdirectory for this iteration
            curr_iter_tmp_dir_par = os.path.join(tmp_dir_par, 'step' + str(self.current_iteration))
            curr_iter_tmp_dir_par = self.sate_team.temp_fs.create_subdir(curr_iter_tmp_dir_par)
            _LOG.debug('directory %s created' % curr_iter_tmp_dir_par)
            break_strategy_index = 0
            this_iter_score_improved = False

            while True:
                break_strategy =  self._get_break_strategy(break_strategy_index)
                if not bool(break_strategy):
                    break
                context_str = "iter%d-%s" % (self.current_iteration, break_strategy)
                # create a subdirectory for this iteration/break_strategy
                curr_tmp_dir_par = os.path.join(curr_iter_tmp_dir_par, break_strategy)
                curr_tmp_dir_par = self.sate_team.temp_fs.create_subdir(curr_tmp_dir_par)
                record_timestamp(os.path.join(curr_tmp_dir_par, 'start_align_timestamp.txt'))
                # Align (with decomposition...)
                self.status('Step %d. Realigning with decomposition strategy set to %s' % (self.current_iteration, break_strategy))
                if self.killed:
                    raise RuntimeError("SATe Job killed")
                tree_for_aligner = self.get_tree_copy()
                new_alignments = []
                for alignment_component in self.alignments:
                    tree_for_aligner = self.get_tree_copy()
                    i = self.alignments.index(alignment_component)
                    aligner = SateAlignerJob(alignment=alignment_component,
                                             sate_team=self.sate_team,
                                             dataset=self.datasets[i],
                                             tree=tree_for_aligner,
                                             tmp_dir_par=curr_tmp_dir_par,
                                             **configuration)
                    self.sate_aligner_job = aligner
                    aligner.launch_alignment(break_strategy=break_strategy,
                                             context_str=context_str)

                    x = aligner.get_results()
                    new_alignments.append(x)
                    self.sate_aligner_job = None
                    del aligner
                (new_alignment, new_partitions) = concatenate_alignments(new_alignments)

                new_alignment.datatype = self.alignments[0].datatype
                record_timestamp(os.path.join(curr_tmp_dir_par, 'start_treeinference_timestamp.txt'))
                # Tree inference
                if self.start_tree_search_from_current:
                    start_from = self.tree
                else:
                    start_from = None
                self.status('Step %d. Alignment obtained. Tree inference beginning...' % (self.current_iteration))
                if self.killed:
                    raise RuntimeError("SATe Job killed")
                tbj = self.sate_team.tree_estimator.create_job(new_alignment,
                                                               starting_tree=start_from,
                                                               num_cpus=self.num_cpus,
                                                               partitions=new_partitions,
                                                               context_str=context_str + " tree",
                                                               tmp_dir_par=curr_tmp_dir_par,
                                                               delete_temps=delete_iteration_temps)
                self.tree_build_job = tbj
                jobq.put(tbj)
                new_score, new_tree_str = tbj.get_results()
                self.tree_build_job = None
                del tbj
                if self.killed:
                    raise RuntimeError("SATe Job killed")

                record_timestamp(os.path.join(curr_tmp_dir_par, 'end_treeinference_timestamp.txt'))
                curr_timestamp = time.time()
                accept_iteration = False

                if self.best_score is None or new_score > self.best_score:
                    self.best_alignment = copy.deepcopy(new_alignment)
                    self.alignments = copy.deepcopy(new_alignments)
                    self.best_tree_str = new_tree_str
                    self.best_score = new_score
                    self.last_improvement_time = curr_timestamp
                    self.num_iter_since_imp = 0
                    this_iter_score_improved = True
                    accept_iteration = True

                if self._get_accept_mode(new_score=new_score, break_strategy_index=break_strategy_index) == AcceptMode.BLIND_MODE:
                    if self.blind_mode_is_final:
                        self.is_stuck_in_blind = True
                        if self.switch_to_blind_timestamp is None:
                            if self._blindmode_trigger:
                                _LOG.debug("Blind runmode trigger = %s" % self._blindmode_trigger)
                            self.switch_to_blind_iter = self.current_iteration
                            self.switch_to_blind_timestamp = curr_timestamp
                            if self.best_alignment is None:
                                self.best_alignment = copy.deepcopy(new_alignment)
                    accept_iteration = True

                if accept_iteration:
                    self.score = new_score
                    self.alignment = copy.deepcopy(new_alignment)
                    self.alignments = copy.deepcopy(new_alignments)
                    self.partitions = new_partitions
                    self.tree_str = new_tree_str
                    self.status('realignment accepted.')
                    # we do not want to continue to try different breaking strategies for this iteration so we break
                    self.status('current score: %s, best score: %s' % (self.score, self.best_score) )
                    break
                else:
                    self.status('realignment NOT accepted.')
                break_strategy_index += 1
                # self.status('current score: %s, best score: %s' % (self.score, self.best_score) )

            if not this_iter_score_improved:
                self.num_iter_since_imp += 1

            self.current_iteration += 1

        if self._termination_trigger:
            _LOG.debug("Termination trigger = %s" % self._termination_trigger)
        record_timestamp(os.path.join(tmp_dir_par, 'end_sateiter_timestamp.txt'))

        ### TODO: if configuration is 'return_final_iter_T&Apair', then skip the following three lines
        if not self.return_final_tree_and_alignment:
            self.alignment = copy.deepcopy(self.best_alignment)
            self.tree_str = self.best_tree_str
            self.score = self.best_score
        else:
            assert self.alignment is not None
            assert self.tree_str is not None
            assert self.score is not None
        # smirarab: I had also added the following commented code. Test to see if needed. 
        # if self.best_alignment is None:
        #    self.best_alignment = copy.deepcopy(new_alignment)


    def status(self, message):
        if self._status_message_func:
            self._status_message_func(message)

    def kill(self):
        self.killed = True
        j = self.sate_aligner_job
        if j:
            self.sate_aligner_job = None
            j.kill()
        j = self.tree_build_job
        if j:
            self.tree_build_job = None
            j.kill()
