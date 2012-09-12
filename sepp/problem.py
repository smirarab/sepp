'''
Created on Sep 12, 2012

@author: smirarab
'''
from scheduler import Job

class Problem(object):
    '''
    This class gives a hierarchy of Problems. Each problem includes a list
    of subproblems, and a pointer to the parent problem. Root problem has no
    parent problem, and tips have empty lists as their children. 
    
    This class provides the data structure required to decompose a problem
    recursively into subproblems. 
    
    A problem can have many Jobs associated with it. These job objects are kept
    in a dictionary. 
    '''
    def __init__(self, parent):        
        self.children = []
        self.parent = parent
        self.level = 0 if parent is None else (parent.level + 1)
        '''
        This dictionary refers to jobs associated with this subproblem,
        using a name.
        '''
        self.jobs = dict()
               
    def add_child(self, subproblem):
        '''
        Add a child job.
        '''
        self.children.append(subproblem)
    
    def iter_leaves(self):
        '''
        Returns an iterator over tips of the problem hierarchy below self.
        '''
        if len(self.children) == 0:
            yield self
        else:
            for c in self.children:
                for tip in c.iter_leaves():
                    yield tip

    def iter_nodes_at_level(self, level):
        '''
        a generator function, returning all Problems under the current problem,
        with a label tag equal to the given level.
        
        Note that this is NOT returning nodes that are 'level' levels below 
        self. 
        '''
        if self.level == level:
            yield self
        else:
            for c in self.children:
                for tip in c.iter_nodes_at_level(level):
                    yield tip

    def get_job_result_by_name(self, jobname, ignore_missing_results = False):
        '''
        returns results of a job, given the job name
        '''
        job = self.jobs[jobname]
        if job is None:
            return None
        assert isinstance(job, Job), "job '%s' is not a Job object" %str(job)
        if job.result_set:
            return job.result
        else:
            raise Exception ("job '%s' has no results set yet: %s" %(jobname,str(job)))

    def add_job(self, jobname, job):
        '''
        A a new job to this problem, to be saved with a name
        '''
        self.jobs[jobname] = job
    
    def get_node_label(self):
        '''
         returns a label for this Problem. Used in __str__ function. 
        '''
        return None
    
    def __str__(self, *args, **kwargs):
        me  = self.get_node_label()
        if len(self.children) == 0:
            return me
        return "(%s)%s" %(",".join((str(child) for child in self.children)),me)
