'''
Created on Sep 19, 2012

@author: smirarab
'''
from sepp.scheduler import Job, JobError
from sepp import get_logger
from abc import abstractmethod, ABCMeta
from subprocess import Popen

import os
import subprocess
import sepp
import stat
import re
from sepp.tree import PhylogeneticTree
import traceback

_LOG = get_logger(__name__)
    
class ExternalSeppJob(Job):
    '''
    All Sepp jobs that run external programs 
    should extend this abstract class.
    This class handles executing external jobs, error handling, and more.     
    '''
    __metaclass__ = ABCMeta
    
    def __init__(self, jobtype, **kwargs):        
        Job.__init__(self)
        self.job_type = jobtype
        self._id = None #_process id for this job
        self._kwargs = dict(kwargs)
        self._process = None
        ''' The following will contain stdout and stderr values *only if* those
        are piped. If not, user should read the files to find the output if needed.
        '''
        self.stdoutdata = self.stderrdata = None        
        self.ignore_error = False # setting this variable tell JobPoll that errors in this job can be ignored when waiting for reults of all jobs to finish
        self.fake_run = False
        self.attributes=dict()
        self.path = sepp.config.options().__getattribute__(self.job_type).path         
    
    def get_id(self):
        return self._id
    id = property(get_id)               

    def get_process(self):
        return self._process
    process = property(get_process)      
    
    def run(self):        
        ''' Runs the external job, and handles errors, piping, etc. 
        get_invocation() needs to be implemented in child classes.
        '''
        if self.fake_run:
            return self.read_results()        
        try:
            _LOG.info("Starting %s Job with input: %s" %(self.job_type, self.characterize_input()))        
            
            assert self.result_set is False, "Job is already run."
            
            #By default discard standard otuput and error from external programs        
            if 'stdout' not in self._kwargs:      
                self._kwargs['stdout'] = subprocess.PIPE
            elif isinstance(self._kwargs['stdout'],str):
                self._kwargs['stdout'] = file(self._kwargs['stdout'], 'w')
    
            if 'stderr' not in self._kwargs:            
                self._kwargs['stderr'] = subprocess.PIPE
                #k['stderr'] = file(os.devnull, 'w')    
            elif isinstance(self._kwargs['stderr'],str):
                self._kwargs['stderr'] = file(self._kwargs['stderr'], 'w')
                    
            self._process = Popen(self.get_invocation(), **self._kwargs)        
            self._id = self._process.pid
            
            (self.stdoutdata, self.stderrdata) = self._process.communicate()         
            
    #        if self.process.stdout:
    #            self.process.stdout.close()
    #        if self.process.stderr:
    #            self.process.stderr.close()
    
    #        except:
    #            raise Exception("\n".join(["Failed to run the following execution:",
    #                                       '" "'.join(self._invocation)])
                      
               
            if self.process.returncode == 0:
                _LOG.info("Finished %s Job with input: %s with:\n"
                      " return code: %s\n output: %s" 
                      %(self.job_type, self.characterize_input(),
                        self.process.returncode, self.stdoutdata))
    
            else:         
                _LOG.info("Finished %s Job with input: %s with:\n"
                      " return code: %s\n output: %s\n error:%s" 
                      %(self.job_type, self.characterize_input(),
                        self.process.returncode, self.stdoutdata, self.read_stderr()))              
                raise JobError("\n".join(["The following execution failed:",
                                           ' '.join(self.get_invocation()),
                                           self.read_stderr() if self.read_stderr() else 'No error messages available']))
                          
    
            return self.read_results()
        except:
            traceback.print_exc()
            raise
        
    def read_stderr(self):    
        '''
        Used for reading standard error when an error is detected. 
        '''    
        if self.stderrdata is not None:
            return self.stderrdata
        elif self._kwargs.has_key("stderr") and isinstance(self._kwargs["stderr"],file):
            return open(self._kwargs["stderr"].name,'r').read()
        else:
            return None
    
    @abstractmethod
    def get_invocation(self):
        '''
        The method needs to return a list with the first argument giving
        the executable, and the rest giving the arguments.
        '''
        raise NotImplementedError("get_invocation should be implemented by subclasses")
            
    @abstractmethod        
    def characterize_input(self):
        '''
        Need to implement this method to help with automatic logging.
        Output a string characterizing the input to this job
        '''
        return ""
    
    @abstractmethod 
    def read_results(self):
        '''
        This method should read the results of an external execution, and turn
        the results into a python object (could be simply the path to an output
        file) and return that python object. This is the result that will be
        ultimately pickled and sent to the main process, and will be accessible 
        to the other jobs, joins, etc. Better not to pass around large files. 
        This results should be picklable. 
        '''
        raise NotImplementedError("read_results should be implemented by subclasses")        
    
    def get_attribute(self,key):
        ''' each job maintains a dictionary of free form attributes. 
        '''
        return self.attributes[key]

    def set_attribute(self,key,val):
        ''' each job maintains a dictionary of free form attributes. 
        '''
        self.attributes[key] = val

class HMMBuildJob(ExternalSeppJob):
    '''
    The Job class that executes a HMM build
    '''

    def __init__(self, **kwargs):
        self.job_type = "hmmbuild"
        ExternalSeppJob.__init__(self, self.job_type, **kwargs)
        self.infile = None #input reference alignment
        self.informat = None #format of input reference alignment
        self.outfile = None #location of output file
        
    def setup(self, infile, outfile, informat = "fasta",**kwargs):
        '''
        Use this to setup the job if you already have input file written to a file.
        Use setup_for_subproblem when possible. 
        '''
        self.infile = infile
        self.informat = informat
        self.outfile = outfile
        self._kwargs = kwargs

    def setup_for_subproblem(self, subproblem ,**kwargs):
        '''
        Automatically sets up a job given a subproblem object. It outputs the
        right alignment subset to a temporary file.
        '''
        assert isinstance(subproblem, sepp.problem.SeppProblem)
        assert isinstance(subproblem.subalignment, sepp.problem.ReadonlySubalignment)
        
        self.infile = sepp.filemgr.tempfile_for_subproblem("hmmbuild.input.", 
                                                       subproblem,
                                                       ".fasta")
        
        subproblem.write_subalignment_without_allgap_columns(self.infile)
        
        self.informat = "fasta"
        self.outfile = sepp.filemgr.tempfile_for_subproblem("hmmbuild.model.", 
                                                       subproblem)
        self._kwargs = kwargs
        
    def get_invocation(self):
        invoc = [self.path, "--symfrac" ,"0.0" ,"--dna"]
        if self._kwargs.has_key("user_options"):
            invoc.extend(self._kwargs["user_options"].split())
        if self.informat == "fasta":
            invoc.extend(['--informat', 'afa'])
        invoc.extend([self.outfile, self.infile])
        return invoc

    def characterize_input(self):
        return self.infile

    def read_results(self):
        '''
        Simply make sure the file exists and is not empty. Don't need to load
        the file into memory or anything else. Just return the location of the
        file. 
        '''
        assert os.path.exists(self.outfile)
        assert os.stat(self.outfile)[stat.ST_SIZE] != 0
        return self.outfile        

class HMMAlignJob(ExternalSeppJob):

    def __init__(self, **kwargs):       
        self.job_type = 'hmmalign'
        ExternalSeppJob.__init__(self, self.job_type, **kwargs)
        self.hmmmodel = None
        self.fragments = None 
        self.outfile = None
        self.base_alignment = None
        self.trim = None
                
    def setup(self,hmmmodel, fragments, output_file, base_alignment=None, trim=True, **kwargs):
        '''
        Setup job parameters when those are externally decided.
        Use setup_for_subproblem when possible.
        '''
        self.hmmmodel = hmmmodel
        self.fragments = fragments 
        self.outfile = output_file
        self.base_alignment = base_alignment
        self.trim = trim
        self._kwargs = kwargs                        
        
    def partial_setup_for_subproblem(self, subproblem, 
                              trim=False, **kwargs):
        '''Automatically sets up a job given a subproblem object. Note that 
        hmmmodel is not set and fragments is just a filename that needs to be
        created later. base_alignment is not set either.  
        '''
        assert isinstance(subproblem, sepp.problem.SeppProblem)
        
        self.outfile = sepp.filemgr.tempfile_for_subproblem("hmmalign.results.",
                                                             subproblem)
        self.fragments = sepp.filemgr.tempfile_for_subproblem("hmmalign.frag.", 
                                                   subproblem,".fasta")
              
        self.trim = trim
        self._kwargs = kwargs  
                    
    def get_invocation(self):
        invoc = [self.path, "--allcol", "--dna",
                 "-o", self.outfile]
        if self.trim:
            invoc.extend(["--trim"])
        #if self.base_alignment:
        #    invoc.extend["--mapali" , self.base_alignment]
        if self._kwargs.has_key("user_options"):
            invoc.extend(self._kwargs["user_options"].split())        
        invoc.extend([self.hmmmodel, self.fragments])
        return invoc

    def characterize_input(self):
        return "model:%s, fragments:%s, trim:%s, base_alignment:%s" %(self.hmmmodel,self.fragments,self.trim, self.base_alignment)

    def read_results(self):
        '''
        Since the output file can be huge, we don't want to read it here, because
        it will need to get pickled and unpickled. Instead, we just send
        back the file name, and will let the caller figure out what to do with it. 
        '''
        if os.path.exists(self.outfile):
            assert os.stat(self.outfile)[stat.ST_SIZE] != 0        
            return self.outfile
        else:
            return None
    
    
    
class HMMSearchJob(ExternalSeppJob):

    def __init__(self, **kwargs):       
        self.job_type = 'hmmsearch'
        ExternalSeppJob.__init__(self, self.job_type, **kwargs)
        self.hmmmodel = None
        self.fragments = None 
        self.outfile = None
        self.elim = None
        self.filters = None
        
    def setup(self, hmmmodel, fragments, output_file, elim=None, filters=True, **kwargs):
        self.hmmmodel = hmmmodel
        self.fragments = fragments 
        self.outfile = output_file
        self.elim = elim
        self.filters = filters
        self._kwargs = kwargs                        
    
    def partial_setup_for_subproblem(self, fragments_file, subproblem, 
                              elim=None, filters=True, **kwargs):
        '''
        Automatically sets up a job given a subproblem object. 
        Note that hmmmodel is not setup and needs to be set separately. 
        '''
        assert isinstance(subproblem, sepp.problem.SeppProblem)
        
        self.outfile = sepp.filemgr.tempfile_for_subproblem("hmmsearch.results.", 
                                                       subproblem)              
        self.fragments = fragments_file
        self.elim = elim
        self.filters = filters
        self._kwargs = kwargs  
    
            
    def get_invocation(self):
        invoc = [self.path, "-o", self.outfile, "--noali", "--cpu", "1"]
        if self.elim is not None:
            invoc.extend(["-E", str(self.elim)])
        if not self.filters:
            invoc.extend(["--max"])
        if self._kwargs.has_key("user_options"):
            invoc.extend(self._kwargs["user_options"].split())        
        invoc.extend([self.hmmmodel, self.fragments])      
        return invoc

    def characterize_input(self):
        return "model:%s, fragments:%s, elim:%s, filter:%s, output:%s" %(self.hmmmodel,self.fragments,self.elim, self.filters, self.outfile)

    def read_results(self):
        '''
           Reads the search output file and returns a dictionary that contains 
           the e-values of the searched fragments
        '''
        assert os.path.exists(self.outfile)
        assert os.stat(self.outfile)[stat.ST_SIZE] != 0        
        outfile = open(self.outfile, 'r');
        results = {}
    
        #Group 1 (e-value) 2 (bitscore) and 9 (taxon name) contain the relevant information, other ones can be ignored unless we plan to do something later
        pattern = re.compile(r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)")
        start_reading = False
        for line in outfile:
            line = line.strip()
            if (not start_reading and line.startswith("E-value") == True):
                start_reading = True
            elif (start_reading and line == ""):
                start_reading = False
                break
            elif (start_reading):             
                matches = pattern.search(line)
                if (matches is not None and matches.group(0).find("--") == -1):
                    results[matches.group(9).strip()] = (float(matches.group(1).strip()), float(matches.group(2).strip()))
        return results
    


class PplacerJob(ExternalSeppJob):

    def __init__(self, **kwargs):       
        self.job_type = 'pplacer'
        ExternalSeppJob.__init__(self, self.job_type, **kwargs)
        '''The following is just a string indicating the type of input provided 
        to pplacer job. Different setup methods can setup this job with different
        types of input (i.e. reference package, separate refernce alignment, etc.)
        '''
        self.setup_setting = None 
        self.tree_file = None         
        self.info_file = None          
        self.extended_alignment_file = None 
        self.out_file = None 
        
    def setup(self, tree_file, info_file, extended_alignment_file, output_file, **kwargs):
        self.tree_file = tree_file        
        self.info_file = info_file         
        self.extended_alignment_file = extended_alignment_file
        self.out_file = output_file
        self._kwargs = kwargs      
        self.setup_setting = "File:TrInEx"                  

    def partial_setup_for_subproblem(self, subproblem, info_file, **kwargs):
        ''' Automatically sets up a job given a subproblem object. 
        Note that extended alignment is just a file name, referring to an empty
        file at this point. This file needs to be created before the job is queued. 
        '''
        assert isinstance(subproblem, sepp.problem.SeppProblem)        
        self.tree_file = sepp.filemgr.tempfile_for_subproblem("pplacer.tree.", 
                                                       subproblem, ".tre")
        self.extended_alignment_file = \
                         sepp.filemgr.tempfile_for_subproblem("pplacer.extended.", 
                                                       subproblem, ".fasta")
        self.out_file = os.path.join(sepp.filemgr.tempdir_for_subproblem(subproblem),
                             self.extended_alignment_file.replace("fasta","json"))                       
        assert isinstance(subproblem.subtree, PhylogeneticTree)
        subproblem.subtree.write_newick_to_path(self.tree_file)
             
        self.info_file = info_file         
        self._kwargs = kwargs      
        self.setup_setting = "File:TrInEx"  

        
    def get_invocation(self):
        invoc = [self.path, 
                 "--out-dir", os.path.dirname(self.out_file)]   
        if self._kwargs.has_key("user_options"):
            invoc.extend(self._kwargs["user_options"].split())
        
        if self.setup_setting == "File:TrInEx":
            invoc.extend(["-s", self.info_file, 
                         "-t", self.tree_file,
                         self.extended_alignment_file])
        return invoc

    def characterize_input(self):
        if self.setup_setting == "File:TrInEx": 
            return "tree_file:%s, info_file:%s, extended alignment:%s, output:%s"\
                %(self.tree_file, self.info_file,self.extended_alignment_file,self.out_file)
        else:
            return "Not setup properly"

    def read_results(self):
        '''
        Since the output file can be huge, we don't want to read it here, because
        it will need to get pickled and unpickled. Instead, we just send
        back the file name, and will let the caller figure out what to do with it. 
        '''
        assert os.path.exists(self.out_file)
        assert os.stat(self.out_file)[stat.ST_SIZE] != 0        
        return self.out_file
        