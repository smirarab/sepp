'''
Created on Oct 2, 2012

@author: smirarab
'''
import unittest
from sepp.jobs import ExternalSeppJob
from sepp.scheduler import JobPool, JobError
from multiprocessing import cpu_count


class TestExternalJob(ExternalSeppJob):
    def tearDown(self):
        # clean up JobPool for other unit tests
        JobPool().terminate()
        JobPool().__init__(cpu_count())

    def __init__(self, pipe=0, **kwargs):
        self.pipe = pipe
        if pipe == 0:
            '''By default errors and output is piped to this process.
            output/err will be available at as self.stdoutdata and
            self.stderrdata and also will be logged'''
            pass
        elif pipe == 1:
            '''The following line forces output to be saved to a file. Note
            that
            1- A file object/descriptor cannot be used here because this is
               executed in parent process, and in child process those file
               handlers will not be available.
            2- self.stdoutdata and self.stderrdata will not be automatically
               filled (to avoid unnecessary I/O). Read stderr and stdout in
               read_results if necessary.
            '''
            kwargs['stdout'] = "./.stdout"
            kwargs['stderr'] = "./.stderr"

        '''always call parent constructor'''
        ExternalSeppJob.__init__(self, "test_find_job", path='dummy', **kwargs)

        '''inputs that need to be set before running this job'''
        self.pattern = None
        self.options = None

    def get_invocation(self):
        return ['find', self.pattern] + (
            [] if self.options is None else self.options.split())

    def characterize_input(self):
        return "pattern: %s options: %s (%s)" % (
            self.pattern, self.options, "pipe: %d" % self.pipe)

    def read_results(self):
        if self.pipe == 2:
            return open("./.stdout").read()
        else:
            if self.stdoutdata is not None:
                return self.stdoutdata
            else:
                return "Nothing captured"


class Test(unittest.TestCase):
    def tearDown(self):
        # clean up JobPool for other unit tests
        JobPool().terminate()
        JobPool().__init__(cpu_count())

    def testSuccess(self):
        find_job = TestExternalJob()
        find_job.pattern = "."
        find_job.options = "-name test*.py"
        JobPool().enqueue_job(find_job)

        JobPool().wait_for_all_jobs()

        res0 = find_job.result.split('\n')[0]

        assert res0 != ""

        find_job = TestExternalJob()
        find_job.pattern = ".."
        find_job.options = "-name %s" % res0.split('/')[-1]
        JobPool().enqueue_job(find_job)

        JobPool().wait_for_all_jobs()

        res1 = find_job.result.split('\n')[0]

        assert res1 != ""

        assert res1.endswith(res0[2:])

    def testError(self):
        find_job = TestExternalJob()
        find_job.pattern = "somerandomdirectorywewillneverhavehere_ordowe"
        ''' Let's ignore this error in subsequent test cases. '''
        find_job.ignore_error = True
        try:
            JobPool().enqueue_job(find_job)
            assert JobPool().wait_for_all_jobs() is False
        except JobError as e:
            assert str(e).find("No such file or directory") != -1, \
                "The error we expected is no such file or directory"

        assert JobPool().get_asynch_result_object(find_job).successful() is \
            False, "We expected the job to fail"

    def testNoPipe(self):
        find_job = TestExternalJob(pipe=1)
        find_job.pattern = "."
        find_job.options = "-name *.py"
        JobPool().enqueue_job(find_job)
        JobPool().wait_for_all_jobs()
        res0 = find_job.result.split('\n')[0]
        assert res0 != ""

        find_job = TestExternalJob(pipe=1)
        find_job.pattern = "somerandomdirectorywewillneverhavehere_ordowe"
        ''' Let's ignore this error in subsequent test cases. '''
        find_job.ignore_error = True
        try:
            JobPool().enqueue_job(find_job)
            assert JobPool().wait_for_all_jobs() is False
            JobPool().get_asynch_result_object(find_job).get()
        except JobError as e:
            assert str(e).find("No such file or directory") != -1, \
                "The error we expected is no such file or directory"

        assert JobPool().get_asynch_result_object(find_job).successful() is \
            False, "We expected the job to fail"


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testBaseExternalJob']
    unittest.main()
