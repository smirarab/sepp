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

class Problem(object):
    
    def __init__(self, parent):        
        self.result = 0
        self.fragments = []
        self.children = []
        self.parent = parent
        self.level = 0 if parent is None else (parent.level + 1)        
        self.jobs = dict(buildmodel=None, searchfragment=None, applymodel=None, summarize=None)

    def get_job_type_result(self, jobtype):
        job = self.jobs[jobtype]
        if job is None:
            return None
        return job.result
                
    def add_child(self, subproblem):
        self.children.append(subproblem)
    
    def tip_iterator(self):
        if len(self.children) == 0:
            yield self
        else:
            for c in self.children:
                for tip in c.tip_iterator():
                    yield tip

    def get_problems_at_level(self, level):
        if self.level == level:
            yield self
        else:
            for c in self.children:
                for tip in c.get_problems_at_level(level):
                    yield tip
                                
    def __str__(self, *args, **kwargs):
        #str(self.level),str(self.result)
        me  = "%s" %(str(self.get_job_type_result("buildmodel")))
        if len(self.children) == 0:
            return me
        return "(%s)%s" %(",".join((str(child) for child in self.children)),me)

class GenericJob(Job):
    def __init__(self, problem):        
        Job.__init__(self)        
        problem.jobs[self.type] = self
        self.problem_name = "level " + str(problem.level)
        #self.add_call_Back(lambda res: problem.add_result_to_problem_object(res, self.type))
        
    def print_result(self,result):
        print >>sys.stderr, "Process [%s]: %s finished with results: %s" %(os.getpid(),self.type,result)
                
                   
class BuildModelJob(GenericJob):    
    def __init__(self, problem):
        self.type = "buildmodel"
        GenericJob.__init__(self,problem)        
        self.add_call_Back(self.print_result)    
    def run(self):
        print >>sys.stderr, "Process [%s]: buildmodel running %s" %(os.getpid(),self.problem_name)                 
        h=0
        step = random()/1000
        for i in xrange(0,100000):
            h+=step*i
            #time.sleep(step/100)
        return int(h)

        
class SearchJob(GenericJob):    
    def __init__(self, problem):
        self.type = "searchfragment"        
        GenericJob.__init__(self,problem)
        self.model = None
        self.fragments = None    
    def run(self):
        print >>sys.stderr, "Process [%s]: searchfragment running %s with model %d" %(os.getpid(),self.problem_name, self.model)  
        #time.sleep(random()/10)    
        #self.state = step
        return [abs(self.model-fragment) for fragment in self.fragments]

class ApplyModelJob(GenericJob):    
    def __init__(self, problem):
        self.type = "applymodel"        
        GenericJob.__init__(self,problem)
        self.fragments = None  
        self.model = None  
    def run(self):
        print >>sys.stderr, "Process [%s]: %s running %s with model %d" %(os.getpid(),self.type,self.problem_name, self.model)  
        #time.sleep(random()/20)    
        #self.state = step
        return [self.model*fragment for fragment in self.fragments]
    
    
class SummarizeJob(GenericJob):
    def __init__(self, problem):
        self.type = "summarize"        
        GenericJob.__init__(self,problem)
        self.resultsPerTipSubproblem = None          
    def run(self):
        print >>sys.stderr, "Process [%s]: %s running %s with results from tips %s" %(os.getpid(),self.type,self.problem_name, str(self.resultsPerTipSubproblem))  
        #time.sleep(random()/20)    
        #self.state = step
        all_fragments = []
        for fragments in self.resultsPerTipSubproblem: all_fragments.extend(fragments)
        fsum=sum(all_fragments)
        return [(fragment+.0)/fsum for fragment in all_fragments]    


class Join_BuildModel_SearchFragment(Join):
    def __init__(self, problem):
        Join.__init__(self)      
        self.grandparent_problem = problem
        '''Add jobs of searchfragment of the grand parent's children (parents)''' 
        for l1  in self.grandparent_problem.children:
            self.add_job(l1.jobs["searchfragment"])
            for l2 in l1.children:
                self.add_job(l2.jobs["buildmodel"])               

    def perform(self):
        print >>sys.stderr, "Process [%s]: Join_BuildModel_SearchFragment joining %s" %(os.getpid(),self.grandparent_problem)
        fragments = self.grandparent_problem.fragments
        frags_range = range(0,len(fragments))
        c1_res = self.grandparent_problem.children[0].get_job_type_result("searchfragment")
        c2_res = self.grandparent_problem.children[1].get_job_type_result("searchfragment")
        model_search = [c1_res[i] - c2_res[i] for i in frags_range]
        self.grandparent_problem.children[0].fragments = [fragments[i] for i in frags_range if model_search[i] < 0]
        self.grandparent_problem.children[1].fragments = [fragments[i] for i in frags_range if model_search[i] >= 0]
        for l1  in self.grandparent_problem.children:
            fragments = l1.fragments                            
            for l2 in l1.children:
                l2j = l2.jobs["searchfragment"]
                l2j.model = l2.get_job_type_result("buildmodel")
                l2j.fragments = fragments 
                JobPool().enqueue_job(l2j)

class Join_ApplyModel_Summarize(Join):
    def __init__(self, problem):
        Join.__init__(self)      
        self.summarylevel_problem = problem
        '''Add jobs of applymodel of the grand parent's children (parents)''' 
        for tip  in self.summarylevel_problem.tip_iterator():
            self.add_job(tip.jobs["applymodel"])
    
    def perform(self):
        print >>sys.stderr, "Process [%s]: Join_ApplyModel_Summarize joining %s" %(os.getpid(),self.summarylevel_problem)
        resultsPerTipSubproblem = []
        for tip  in self.summarylevel_problem.tip_iterator():
            resultsPerTipSubproblem.append(tip.get_job_type_result("applymodel"))
        self.summarylevel_problem.jobs["summarize"].resultsPerTipSubproblem = resultsPerTipSubproblem
        JobPool().enqueue_job(self.summarylevel_problem.jobs["summarize"])
          

class Join_tip_searchfragment(Join):
    def __init__(self,root_problem):
        Join.__init__(self)
        self.root_problem = root_problem            
        for p in root_problem.tip_iterator():
            self.add_job(p.jobs["searchfragment"])         
    
    def perform(self):
        print >>sys.stderr, "Process [%s]: Join_tip_searchfragment joining %s" %(os.getpid(),self.root_problem)                
        def print_fragments(problem):
            print "level " + str(problem.level), str(problem.get_job_type_result("buildmodel")),  problem.fragments
            for c in problem.children:
                print_fragments(c)                 
        print_fragments(self.root_problem)
        
        for p in root_problem.tip_iterator():
            j = p.jobs["applymodel"]
            j.fragments = p.jobs["searchfragment"].fragments
            j.model = p.jobs["searchfragment"].model
            JobPool().enqueue_job(j)
                
    def __str__(self):
        return "join Join_tip_searchfragment for", self.root_problem            

                        
def connect_buildmodel_and_searchfragment_jobs(problem):
    if len(problem.children) > 0: # and sum((len(c.children) for c in problem.children)) > 0:
        Join_BuildModel_SearchFragment(problem)
        for c in problem.children:
            connect_buildmodel_and_searchfragment_jobs(c)
                
def build_subproblems(problem = None):
    if problem is None:
        problem = Problem(None)
        problem.fragments = [int(random()*10000000) for i in range(1,2000)] 
    if problem.level < 8:
        subproblem = Problem(problem)
        build_subproblems(subproblem)
        problem.add_child(subproblem)
        subproblem = Problem(problem)
        build_subproblems(subproblem)
        problem.add_child(subproblem)
        #problem.ad
    return problem
    
def recursive_add_buildmodel_and_searchfragment_job(problem):
    if problem.parent is not None:
        BuildModelJob(problem)
        SearchJob(problem)
    for subp in problem.children:
        recursive_add_buildmodel_and_searchfragment_job(subp)                          
s = 0    
lock = Lock()
if __name__ == '__main__':
    pool = JobPool(2)
    
    root_problem = build_subproblems()    
    
    recursive_add_buildmodel_and_searchfragment_job (root_problem)      
    connect_buildmodel_and_searchfragment_jobs(root_problem)    
    Join_tip_searchfragment(root_problem)
    for problem in root_problem.tip_iterator():
        ApplyModelJob(problem)
    for problem in root_problem.get_problems_at_level(3):
        SummarizeJob(problem)
        Join_ApplyModel_Summarize(problem)
    
    for c in root_problem.children:
        def enq_job_searchfragment(result, next_job):
            next_job.model = result
            next_job.fragments = root_problem.fragments
            JobPool().enqueue_job(next_job)
        c.jobs["buildmodel"].add_call_Back(lambda result, next_job = c.jobs["searchfragment"]: enq_job_searchfragment(result, next_job))

    def enqueue_buildmodel_job(problem):
        if problem.parent is not None:
            JobPool().enqueue_job(problem.jobs["buildmodel"])
        for child in problem.children:
            enqueue_buildmodel_job(child)
            
    enqueue_buildmodel_job(root_problem)
    
    JobPool().wait_for_all_jobs()     
    
    
    for problem in root_problem.get_problems_at_level(3):
        print  problem.get_job_type_result("buildmodel"), problem.get_job_type_result("summarize")      
    #print [str(x) for x in root_problem.tip_iterator()]
    
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
#        job.add_call_Back(lambda result, parent = job, join = joins[0]: add_a_searchfragment_child(parent))
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
