"""
Created on Aug 22, 2012

This is meant to emulate a Hierarchical SEPP. We have four types of jobs that
emulate SEPP jobs. These are called in parallel using a a DAG that would be
structurally identical to the hierarchical SEPP and can map to that problem
directly.

@author: smirarab
"""
from sepp.scheduler import Job, JobPool, Join
from random import random
import sys
import os
from multiprocessing import Lock, cpu_count
from sepp.problem import Problem
import unittest


# Depth of problem hierarchy (for decomposition). Tips are equivalent
# of alignment subsets
DEPTH = 2
# The level subproblems are aggregated (like placement subsets)
SUMMERIZE_LEVEL = 1


def trimstr(i):
    s = str(i)
    return s if len(s) < 101 else "%s ..." % (s[0:100])


class TestProblem(Problem):
    def __init__(self, parent):
        Problem.__init__(self, parent)
        self.jobs = dict(buildmodel=None, searchfragment=None, applymodel=None,
                         summarize=None)
        self.fragments = []

    def get_node_label(self):
        # str(self.level),str(self.result)
        return "%s" % (str(self.get_job_result_by_name("buildmodel")))


class GenericJob(Job):
    """
    Base class of all other job types
    """
    def __init__(self, problem):
        Job.__init__(self)
        problem.add_job(self.type, self)
        self.problem_name = "level " + str(problem.level)
        # self.add_call_Back(lambda res: problem.add_result_to_probl
        # em_object(res, self.type))

    def print_result(self, result):
        print("Process [%s]: %s finished with results: %s" % (
              os.getpid(), self.type, trimstr(result)), file=sys.stderr)


class BuildModelJob(GenericJob):
    """
    Build a Model (equivalent of HMMERBUILD)
    """
    def __init__(self, problem):
        self.type = "buildmodel"
        GenericJob.__init__(self, problem)
        self.add_call_Back(self.print_result)
        # self.someclass

    def run(self):
        """
        The model is going to be a random number calculated in an unnecessarily
        lengthy fashion.
        """
        print("Process [%s]: buildmodel running %s" % (
            os.getpid(), self.problem_name), file=sys.stderr)
        h = 0
        step = random() / 10000
        for i in range(0, 10000):
            h += step * i
            # time.sleep(step/100)

        return int(h)


class SearchJob(GenericJob):
    def __init__(self, problem):
        """
        Search fragments against a Model (equivalent of HMMERSEARCH).

        This jobs depend on BuildModel job of the same subproblem, and also
        Search job __parent problems (if any)
        """
        self.type = "searchfragment"
        GenericJob.__init__(self, problem)
        '''
        Both the model and fragments list are inputs to this jobs. These are
        not calculated yet at the time the job object is *created*. But they
        should be set before the job in *enqueued*. By the time this job is
        enqueued, the dependencies are run, and hence the model is available
        (from BuildModel) and fragments are known (from Search of parents)
        '''
        self.model = None
        self.fragments = None

    def run(self):
        """
        Simply find the difference between each fragment value and the model
        associated with this job. Do this in an inefficient way.
        """
        print("Process [%s]: %s running %s with model %d" % (
            os.getpid(), self.type, self.problem_name, self.model),
            file=sys.stderr)
        # time.sleep(random()/10)
        # self.state = step
        for i in range(1, 2000):

            ret = [abs(self.model-fragment) for fragment in self.fragments]
        return ret


class ApplyModelJob(GenericJob):
    def __init__(self, problem):
        """
        Apply the model to fragments associated with each subproblem
        (equivalent of HMMERALIGN)

        This job is run only on tip subproblems, and depends on the Search
        jobs. Input again is a model and a set of fragments. These are both
        available only after Search job is run. These values should be set
        before an ApplyModel job is enqueued.
        """
        self.type = "applymodel"
        GenericJob.__init__(self, problem)
        self.fragments = None
        self.model = None

    def run(self):
        """
        Do some mathematical calculation involving each fragment and the model.
        """
        print("Process [%s]: %s running %s with model %d" % (
            os.getpid(), self.type, self.problem_name, self.model),
            file=sys.stderr)
        # time.sleep(random()/20)
        # self.state = step
        for i in range(1, 2000):
            ret = [self.model * fragment / 2000.0
                   for fragment in self.fragments]

        return ret


class SummarizeJob(GenericJob):
    def __init__(self, problem):
        """
        Use model-applied fragment values (equivalent of extended alignments)
        to compute some new results (equivalent of pplacer)

        This jobs
        """
        self.type = "summarize"
        GenericJob.__init__(self, problem)
        self.resultsPerTipSubproblem = None

    def run(self):
        """
        Simply normalize the model-applied fragments across
        """
        print("Process [%s]: %s running %s with results from tips %s" % (
                os.getpid(), self.type, self.problem_name,
                trimstr(self.resultsPerTipSubproblem)), file=sys.stderr)
        # time.sleep(random()/20)
        # self.state = step
        all_fragments = []
        for fragments in self.resultsPerTipSubproblem:
            all_fragments.extend(fragments)
        for i in range(1, 20000):
            fsum = sum(all_fragments)
            ret = [(fragment + .0) / fsum for fragment in all_fragments]
        return ret


class Join_BuildModel_SearchFragment(Join):
    """
    A join object used to join buildmodel jobs with their __parent
    searchfragments jobs. The input to the constructor is the
    grandparent problem; i.e., the problem two level above buildmodel jobs,
    and one level above searchfragment jobs to be joined together
    """
    def __init__(self, problem):
        Join.__init__(self)
        self.grandparent_problem = problem
        '''automatically add appropriate jobs to this join'''
        for l1 in self.grandparent_problem.children:
            self.add_job(l1.jobs["searchfragment"])
            for l2 in l1.children:
                self.add_job(l2.jobs["buildmodel"])

    def perform(self):
        print("Process [%s]: Join_BuildModel_SearchFragment joining %s" % (
            os.getpid(), self.grandparent_problem), file=sys.stderr)
        '''
        1 - start from grandparent fragments.
        2 - Based on results from joined search operations, figure out for each
            grandparent fragments whether it is closer to parent1 or parent2
            (grandparent's children), and divide the fragments accordingly
        3 - Set fragments attribute of parent1 and parent2 based on preceding
            calculation
        4 - For each of grandchildrens of the grandparent problem,
            4-1 Set the fragments attribute of their searchfragment job to the
                fragments of its parents (set in step 3)
            4-2 Set the model attribute of their searchfragment job to the
                model computed in its buildmodel job
            4-3 Enqueue its searchfragment job

        The above procedure is equivalent of figuring out the child HMM model
        that a fragment is closer to, and classifying it accordingly.
        '''
        fragments = self.grandparent_problem.fragments
        frags_range = list(range(0, len(fragments)))
        c1_res = self.grandparent_problem.children[0].get_job_result_by_name(
            "searchfragment")
        c2_res = self.grandparent_problem.children[1].get_job_result_by_name(
            "searchfragment")
        model_search = [c1_res[i] - c2_res[i] for i in frags_range]
        self.grandparent_problem.children[0].fragments = [
            fragments[i] for i in frags_range if model_search[i] < 0]
        self.grandparent_problem.children[1].fragments = [
            fragments[i] for i in frags_range if model_search[i] >= 0]
        for l1 in self.grandparent_problem.children:
            fragments = l1.fragments
            for l2 in l1.children:
                l2j = l2.jobs["searchfragment"]
                l2j.model = l2.get_job_result_by_name("buildmodel")
                l2j.fragments = fragments
                JobPool().enqueue_job(l2j)


class Join_ApplyModel_Summarize(Join):
    """
    This joins applymodel jobs with summirze jobs
    """
    def __init__(self, problem):
        Join.__init__(self)
        self.summarylevel_problem = problem
        '''Add jobs of applymodel of the grand __parent's children (parents)'''
        for tip in self.summarylevel_problem.iter_leaves():
            self.add_job(tip.jobs["applymodel"])

    def perform(self):
        """
        Aggregate fragments from tips to the SUMMERIZE_LEVEL level problem,
        and enqueue a summarize job
        """
        print("Process [%s]: Join_ApplyModel_Summarize joining %s" % (
            os.getpid(), self.summarylevel_problem), file=sys.stderr)
        resultsPerTipSubproblem = []
        for tip in self.summarylevel_problem.iter_leaves():
            resultsPerTipSubproblem.append(tip.get_job_result_by_name(
                "applymodel"))
        self.summarylevel_problem.jobs["summarize"].resultsPerTipSubproblem = \
            resultsPerTipSubproblem
        JobPool().enqueue_job(self.summarylevel_problem.jobs["summarize"])


class Join_tip_searchfragment(Join):
    """
    After all search jobs have finished on tips, we need to start applying
    the model of each tip subproblem to its fragments. This join takes care of
    that.
    """
    def __init__(self, root_problem):
        Join.__init__(self)
        self.root_problem = root_problem
        for p in root_problem.iter_leaves():
            self.add_job(p.jobs["searchfragment"])

    def perform(self):
        """
        First print out some summary of everything up to here.
        Then update applymodel jobs with correct fragment and model, and
        then enqueue them.
        """
        print("Process [%s]: Join_tip_searchfragment joining %s" % (
            os.getpid(), self.root_problem), file=sys.stderr)

        def print_fragments(problem):
            print(
                "\t" * (problem.level),
                trimstr(problem.get_job_result_by_name("buildmodel")),
                trimstr(problem.fragments))
            for c in problem.children:
                print_fragments(c)
        print_fragments(self.root_problem)

        for p in root_problem.iter_leaves():
            j = p.jobs["applymodel"]
            j.fragments = p.jobs["searchfragment"].fragments
            j.model = p.jobs["searchfragment"].model
            JobPool().enqueue_job(j)

    def __str__(self):
        return "join Join_tip_searchfragment for", self.root_problem


def build_subproblems(problem=None):
    """ Makes a fully balanced problem (binary) tree upto level DEPTH"""
    if problem is None:
        problem = TestProblem(None)
        ''' Root problem needs to have fragments assigned to it'''
        problem.fragments = [int(random() * 1000000) for i in range(1, 2000)]
    if problem.level < DEPTH:
        subproblem = TestProblem(problem)
        build_subproblems(subproblem)
        problem.add_child(subproblem)
        subproblem = TestProblem(problem)
        build_subproblems(subproblem)
        problem.add_child(subproblem)
        # problem.ad
    return problem


def build_job_dag(root_problem):
    def recursive_add_buildmodel_and_searchfragment_job(problem):
        """build a search job and a buildmodel job for each problem in
        the problem tree."""
        if problem.parent is not None:
            BuildModelJob(problem)
            SearchJob(problem)
        for subp in problem.children:
            recursive_add_buildmodel_and_searchfragment_job(subp)

    def connect_buildmodel_and_searchfragment_jobs(problem):
        """
        Connect all searchfragment jobs for nodes under a certain node
        ("grandparent") with all the buildmodels jobs under all those nodes
        (grandchildren).
        """
        # and sum((len(c.children) for c in problem.children)) > 0:
        if len(problem.children) > 0:
            Join_BuildModel_SearchFragment(problem)
            for c in problem.children:
                connect_buildmodel_and_searchfragment_jobs(c)

    recursive_add_buildmodel_and_searchfragment_job(root_problem)
    connect_buildmodel_and_searchfragment_jobs(root_problem)

    ''' All the tip search problems should be joined'''
    Join_tip_searchfragment(root_problem)

    ''' only for leaf problems create a apply model job'''
    for problem in root_problem.iter_leaves():
        ApplyModelJob(problem)
    ''' only for problems at SUMMERIZE_LEVEL create a summarize job'''
    for problem in root_problem.iter_nodes_at_level(SUMMERIZE_LEVEL):
        SummarizeJob(problem)
        '''join these summarizejobs with all the leave applymodel
           jobs under them'''
        Join_ApplyModel_Summarize(problem)

    '''
    highest level needs special handling. Make sure you add a call_back
    function to call searchfragment after buildmodel for the top level
    buildmodel jobs. Note that this could have been also achived using a join
    that has only one job.'''
    for c in root_problem.children:
        def enq_job_searchfragment(result, next_job):
            next_job.model = result
            next_job.fragments = root_problem.fragments
            JobPool().enqueue_job(next_job)
        c.jobs["buildmodel"].add_call_Back(
            lambda result, next_job=c.jobs["searchfragment"]:
                enq_job_searchfragment(result, next_job))


s = 0
lock = Lock()
root_problem = None


def run():
    global root_problem
    JobPool().terminate()
    JobPool().__init__(2)
    JobPool(2)

    '''build the problem structure'''
    root_problem = build_subproblems()

    '''build the dat of jobs'''
    build_job_dag(root_problem)

    '''All buildmodel jobs are ready to be started (i.e. no dependency).
    Queue them up. Once they run, they will automatically enqueue the rest of
    the DAG through joins and callbacks '''
    def enqueue_buildmodel_job(problem):
        if problem.parent is not None:
            JobPool().enqueue_job(problem.jobs["buildmodel"])
        for child in problem.children:
            enqueue_buildmodel_job(child)
    enqueue_buildmodel_job(root_problem)

    '''Wait for all jobs to finish'''
    JobPool().wait_for_all_jobs()

    ''' print out the results of summarize jobs. We could merge the
    results from all summarize jobs here if we wanted '''
    for problem in root_problem.iter_nodes_at_level(SUMMERIZE_LEVEL):
        print(trimstr(problem.get_job_result_by_name("buildmodel")),
              trimstr(problem.get_job_result_by_name("summarize")))
    # print [str(x) for x in root_problem.iter_leaves()]


class Test(unittest.TestCase):
    def tearDown(self):
        # clean up JobPool for other unit tests
        JobPool().terminate()
        JobPool().__init__(cpu_count())

    def test_me(self):
        run()


if __name__ == '__main__':
    unittest.main()
