'''
Created on Aug 22, 2012

@author: smirarab
'''
from scheduler import Job, JobPool, Join
from random import random
import sys
import os
from multiprocessing import Lock
import time

class Subproblem(object):
    
    def __init__(self, parent):        
        self.result = 0
        self.fragments = []
        self.children = []
        self.parent = parent
        self.level = 0 if parent is None else (parent.level + 1)        
        self.jobs = dict(type1=None,type2=None,type3=None,type4=None)

    def get_job_type_result(self, jobtype):
        job = self.jobs[jobtype]
        if job is None:
            return None
        return job.result
                
    def add_child(self, subproblem):
        self.children.append(subproblem)
    
    def get_tip_problems(self):
        if len(self.children) == 0:
            yield self
        else:
            for c in self.children:
                for tip in c.get_tip_problems():
                    yield tip
    
#    def add_result_to_problem_object(self, result, jobtype):
#        self.job_results[jobtype] = result
            
    def __str__(self, *args, **kwargs):
        me  = "%s:%s:%s" %(str(self.level),str(self.result),str(self.get_job_type_result("type1")))
        if len(self.children) == 0:
            return me
        return "(%s)%s" %(",".join((str(child) for child in self.children)),me)

class GenericJob(Job):
    def __init__(self, problem):        
        Job.__init__(self)        
        problem.jobs[self.type] = self
        self.problem_name = str(problem)
        #self.add_call_Back(lambda res: problem.add_result_to_problem_object(res, self.type))
        
                   
class Type1Job(GenericJob):    
    def __init__(self, problem):
        self.type = "type1"
        GenericJob.__init__(self,problem)        
        #self.subproblem = sp  
        #self.add_call_Back(self.print_res)
        #self.state = None
    
    def run(self):
        print >>sys.stderr, "Process [%s]: typ1 running %s" %(os.getpid(),self.problem_name)                 
        time.sleep(random()/10)
        return int(random()*100)

class Type2Job(GenericJob):    
    def __init__(self, problem):
        self.type = "type2"        
        GenericJob.__init__(self,problem)
        self.model = None
        self.fragments = None
        
        #self.subproblem = sp  
        #self.add_call_Back(self.print_res)
        #self.state = None
    
    def run(self):
        print >>sys.stderr, "Process [%s]: typ2 running %s" %(os.getpid(),self.problem_name)
        #self.model = self.subproblem.job_results["type1"]        
        time.sleep(random()/2)    
        #self.state = step
        return [abs(self.model-fragment) for fragment in self.fragments]
            
def build_subproblems(problem = None):
    if problem is None:
        problem = Subproblem(None)
        problem.fragments = [int(random()*100) for i in range(1,20)] 
    if problem.level < 5:
        subproblem = Subproblem(problem)
        build_subproblems(subproblem)
        problem.add_child(subproblem)
        subproblem = Subproblem(problem)
        build_subproblems(subproblem)
        problem.add_child(subproblem)
        #problem.ad
    return problem
    
def recursive_add_type1_and_type2_job(problem):
    if problem.parent is not None:
        Type1Job(problem)
        Type2Job(problem)
    for subp in problem.children:
        recursive_add_type1_and_type2_job(subp)

class Join_type1_type2(Join):
    def __init__(self, problem):
        Join.__init__(self)      
        self.grandparent_problem = problem
        '''Add jobs of type2 of the grand parent's children (parents)''' 
        for l1  in self.grandparent_problem.children:
            self.add_job(l1.jobs["type2"])
            for l2 in l1.children:
                self.add_job(l2.jobs["type1"])               

    def perform(self):
        #print "A join on", self.grandparent_problem
        fragments = self.grandparent_problem.fragments
        frags_range = range(0,len(fragments))
        c1_res = self.grandparent_problem.children[0].get_job_type_result("type2")
        c2_res = self.grandparent_problem.children[1].get_job_type_result("type2")
        model_search = [c1_res[i] - c2_res[i] for i in frags_range]
        self.grandparent_problem.children[0].fragments = [fragments[i] for i in frags_range if model_search[i] < 0]
        self.grandparent_problem.children[1].fragments = [fragments[i] for i in frags_range if model_search[i] >= 0]
        for l1  in self.grandparent_problem.children:
            fragments = l1.fragments                            
            for l2 in l1.children:
                l2j = l2.jobs["type2"]
                l2j.model = l2.get_job_type_result("type1")
                l2j.fragments = fragments 
                JobPool().enqueue_job(l2j)
        
def connect_type1_and_typ2_jobs(problem):
    if len(problem.children) > 0:# and sum((len(c.children) for c in problem.children)) > 0:
        Join_type1_type2(problem)
        for c in problem.children:
            connect_type1_and_typ2_jobs(c)
            

class Join_tip_type2(Join):
    def __init__(self,root_problem):
        Join.__init__(self)
        self.root_problem = root_problem            
        for p in root_problem.get_tip_problems():
            self.add_job(p.jobs["type2"])         
    
    def perform(self):                
        def print_fragments(problem):
            print str(problem.get_job_type_result("type1")), problem.level, problem.fragments
            for c in problem.children:
                print_fragments(c) 
                
        print_fragments(self.root_problem)
                
    def __str__(self):
        return "join tip_type for", self.root_problem            
                        
#def add_a_type2_child(parent,join):
#    #print >>sys.stderr, "Adding a child job for %s" %(parent)
#    new_job = GenericJob()
#    if join is not None:
#        join.replace_job(parent, new_job)
#    JobPool().enqueue_job(new_job)
#    #print "Added a child for: ",parent       
        
s = 0    
lock = Lock()
if __name__ == '__main__':
    pool = JobPool(2)
    
    root_problem = build_subproblems()
    
    recursive_add_type1_and_type2_job (root_problem)
    
    connect_type1_and_typ2_jobs(root_problem)
    
    Join_tip_type2(root_problem)
    
    for c in root_problem.children:
        def enq_job_type2(result, next_job):
            next_job.model = result
            next_job.fragments = root_problem.fragments
            JobPool().enqueue_job(next_job)
        c.jobs["type1"].add_call_Back(lambda result, next_job = c.jobs["type2"]: enq_job_type2(result, next_job))

    def enqueue_type1_job(problem):
        if problem.parent is not None:
            JobPool().enqueue_job(problem.jobs["type1"])
        for child in problem.children:
            enqueue_type1_job(child)
            
    enqueue_type1_job(root_problem)
    
    JobPool().wait_for_all_jobs()            
    #print [str(x) for x in root_problem.get_tip_problems()]
    
#    jobs = {}
#    joins = []
#    for level in range(1,3):
#        jobs[level] = []
#        for i in range(0,2 << (level-1)):
#            job = TestJob(level,1,i,0)            
#            jobs[level].append(job)
#        if level> 1: 
#            this_and_previous = jobs[level] + jobs[level-1]
#            join = Join_1_2()
#            join.add_jobs(this_and_previous)
#            joins.append[join]
#            
#    for job in jobs[1]:
#        job.add_call_Back(lambda result, parent = job, join = joins[0]: add_a_type2_child(parent))
#        
#    pool.enqueue_job(job)
#    
#    
#    sample_job = pool.get_asynch_result_object(jobs[3])
#    
#    #pool.terminate()
#    
#    pool.wait_for_all_jobs(ignore_error=True)
#    
#    # Test one of the jobs, to see if it is successful
#    if sample_job.ready() and sample_job.successful():
#        assert jobs[3].callbacks_finished == True
#    else:
#        assert jobs[3].callbacks_finished == False
#    
#    errors = pool.get_all_job_errors()
#    print "Following job errors were raised:", errors 
#    
#    try:
#        pool.wait_for_all_jobs(ignore_error=False)
#    except Exception as e:
#        print "Seems we have some jobs that failed (expected): ", e    
#    
#    errs = [pool.get_job_error(job) for job in pool.get_failed_jobs()]
#    
#    assert len(errs) == len(errors) and False not in [x in errors for x in errs]
#    
#    #print [job.state for job in jobs]
#    print "Number of started jobs - number of printed results:", s
#    print "Number of failed jobs:", len(errors)
#    assert s == len (errors), "Parallelization Error, what happened to the rest?"
