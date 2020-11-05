"""
Created on Aug 22, 2012

@author: smirarab
"""

from sepp.scheduler import Job, JobPool
from multiprocessing import cpu_count, set_start_method
from random import random
from multiprocessing import Lock
import time
import unittest

set_start_method('fork')


class TestJob(Job):
    def __init__(self, jobname):
        Job.__init__(self)
        global s
        lock.acquire()
        s += 1
        lock.release()
        self.jobname = jobname
        self.state = None

        def add_a_child(parent):
            # print("Adding a child job for %s" % (parent), file=sys.stderr)
            JobPool().enqueue_job(TestJob("%s.child" % parent))
            # print "Added a child for: ",parent

        def a_very_bad_callback():
            raise Exception("A bad callback raises an exception.")

        ''' Make a random number of these jobs, child-bearing'''
        if random() < 0.75:
            self.add_call_Back(
                lambda result, parent=self.jobname: add_a_child(parent))

        if random() < 0.1:
            self.add_call_Back(lambda result: a_very_bad_callback())

        self.add_call_Back(self.print_res)

    def print_res(self, result):
        global s
        # print("%s returned: %d" % (self.jobname, self.result))
        lock.acquire()
        s -= 1
        lock.release()

    def run(self):
        # print("Process [%s]: running %s" % (os.getpid(), self.jobname),
        #       file=sys.stderr)
        # Adding the following line results in a failure:
        # AssertionError: daemonic processes are not allowed to have children.
        # This is expected because the new process is going to try to start a
        # new JobPool, which is not allowed.
        # print "cpus ",JobPool().cpus
        h = 0
        step = random()
        if step < 0.1:
            raise Exception("Some (truly) random error occurred in job %s." %
                            self.jobname)
        for i in range(0, 100):
            h += step * i
            time.sleep(step / 100)
        # self.state = step
        return h


s = 0
lock = Lock()


def run():
    global pool
    pool1 = JobPool(2)
    pool2 = JobPool()
    if pool1 != pool2:
        raise Exception("hmmm, I thought JobPool is 'Singleton'")
    try:
        JobPool(4)
    except Exception as e:
        print(("As expected, making a new JobPool with a"
               " different cpu count failed: %s") % e)

    pool = JobPool()
    jobs = []
    for j in range(1, 20):
        job = TestJob(str(j))
        jobs.append(job)
        pool.enqueue_job(job)

    sample_job = pool.get_asynch_result_object(jobs[3])

    # pool.terminate()

    pool.wait_for_all_jobs(ignore_error=True)

    # Test one of the jobs, to see if it is successful
    if sample_job.ready() and sample_job.successful():
        assert(jobs[3].result_set is True)
    else:
        assert(jobs[3].result_set is False)

    errors = pool.get_all_job_errors()
    # print("Following job errors were raised:", errors)

    try:
        pool.wait_for_all_jobs(ignore_error=False)
    except Exception as e:
        print("Seems we have some jobs that failed (expected): ", e)

    errs = [pool.get_job_error(job) for job in pool.get_failed_jobs()]

    # print(errs)

    assert len(errs) == len(errors), \
        "Number of errors from failed jobs: %d. Number of errors: %d" % \
        (len(errs), len(errors))
    assert False not in [x in errors for x in errs]

    # print [job.state for job in jobs]
    # print("Number of started jobs - number of printed results:", s)
    # print("Number of failed jobs:", len(errors))
    assert s == len(errors), \
        "Parallelization Error, what happened to the rest?"

    # print("Everything seems fine!")


class Test(unittest.TestCase):
    def tearDown(self):
        # clean up JobPool for other unit tests
        JobPool().terminate()
        JobPool().__init__(cpu_count())

    def test_me(self):
        run()


if __name__ == '__main__':
    unittest.main()
