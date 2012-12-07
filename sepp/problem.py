'''
Created on Sep 12, 2012

@author: smirarab
'''
from scheduler import Job
from sepp.alignment import ReadonlySubalignment, ExtendedAlignment
from sepp import get_logger

_LOG = get_logger(__name__)

class Problem(object):
    '''
    This class gives a hierarchy of Problems. Each problem includes a list
    of subproblems, and a pointer to the __parent problem. Root problem has no
    __parent problem, and tips have empty lists as their __children. 
    
    This class provides the data structure required to decompose a problem
    recursively into subproblems. 
    
    A problem can have many Jobs associated with it. These job objects are kept
    in a dictionary. 
    '''
    def __init__(self, parent):        
        self.__children = []
        self.__parent = parent
        if parent:
            if not isinstance(parent, Problem):
                raise TypeError("Expecting Problem but founding %s instead" %type(parent))
            parent.add_child(self)
        '''
        This dictionary refers to jobs associated with this subproblem,
        using a name.
        '''
        self.jobs = dict()
        self.annotations = {}
    
    def get_parent(self):
        return self.__parent    
    parent = property(get_parent)
               
    
    def get_level(self):
        return 0 if self.parent is None else (self.parent.level + 1)
    level = property(get_level)
               
    def get_children(self):
        return self.__children
    children = property(get_children)           
    
    def add_child(self, subproblem):
        '''
        Add a child job.
        '''
        self.__children.append(subproblem)        
        
    def iter_leaves(self):
        '''
        Returns an iterator over tips of the problem hierarchy below self.
        '''
        if len(self.__children) == 0:
            yield self
        else:
            for c in self.__children:
                for tip in c.iter_leaves():
                    yield tip

    def iter_nodes_at_level(self, level):
        '''
        a generator function, returning all Problems under the current problem,
        with a label tag equal to the given __level.
        
        Note that this is NOT returning nodes that are '__level' levels below 
        self. 
        '''
        if self.__level == level:
            yield self
        else:
            for c in self.__children:
                for tip in c.iter_nodes_at_level(level):
                    yield tip

    def add_job(self, jobname, job):
        '''
        A a new job to this problem, to be saved with a name
        '''
        self.jobs[jobname] = job

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
    
    def get_node_label(self):
        '''
         returns a label for this Problem. Used in __str__ function. 
        '''
        return None
    
    def get_path_to_root(self):
        ret = [self]
        n = self
        while n.parent is not None:
            ret.append(n.parent)
            n = n.parent
        return ret 
        
    def __str__(self, *args, **kwargs):
        me  = self.get_node_label()
        if len(self.__children) == 0:
            return me
        return "(%s)%s" %(",".join((str(child) for child in self.__children)),me)
    
    
class SeppProblem(Problem):
    ''' A typical Sepp subproblem, defined by a set of taxa, and a parent
    
        Adds few basic attributes of a Sepp Problem to the Problem class:
        1- a list of taxa associated with this sub problem
        2- a Readonly subalignment of the parent problem's alginemnt, induced by this probelms's taxa (with all gaps columns included)
        3- the subtree of parent problem's tree, induced by this problem's taxa. This is a deep copy. 
        4- a set of fragments associated with this subproblem
        
        List of taxa is a required field, provided upon initiailzation. 
        Subalignment is automatically computed, and cannot be changed. 
        Subtree is automatically computed, but can be overwriteen by the user
        Fragments is not automatically computed, and needs to be set by the user
    '''     
    def __init__(self, taxa, parent=None):
        Problem.__init__(self, parent)
        self.__taxa = taxa
        self.__subalignment = None
        self.__subtree = None
        self.fragments = None   
        self.label = None             
    
    def get_taxa(self):
        return self.__taxa
    taxa = property(get_taxa)
    
    def get_subalignment(self):
        '''
        If subalignment is not set before, automatically builds a readonly 
        subalignment for this subproblem, based on the taxa assigned to this subproblem.
        Otherwise, returns the saves subalignment.  
        '''
        if self.__subalignment == None:
            if isinstance(self.parent, SeppProblem) or hasattr(self.parent, "subalignment"):
                self.__subalignment = ReadonlySubalignment(self.taxa, self.parent.subalignment)
                #print sorted(self.__subalignment.keys())
                #print sorted(self.parent.subalignment.keys())
                #print "..... %s %s" %(self.label, self.parent.label) 
            else:
                raise TypeError("parent object %s has no subalignment attribute" %self.parent)        
        return self.__subalignment    
    def set_subalignment(self, alignment):
        self.__subalignment = alignment
    subalignment = property(get_subalignment, set_subalignment)
    
    def get_subtree(self):
        '''
        If subtree is not assigned before, automatically buils a subtree based
        on the taxa assigned to this subproblem. Otherwise, returns the saves
        subtree.
        '''
        if self.__subtree == None:
            if isinstance(self.parent, SeppProblem) or hasattr(self.parent, "subtree"):
                self.__subtree = self.parent.subtree.get_subtree(self.taxa)
            else:
                raise TypeError("parent object %s has no subtree attribute" %self.parent)
        return self.__subtree    
    def set_subtree(self, subtree):
        self.__subtree = subtree
    subtree = property(get_subtree,set_subtree)
    
    def get_node_label(self):
        return self.label
    
    def write_subalignment_without_allgap_columns(self, path):
        '''
        Writes out the subalignment associated with this subproblem to a file,
        removing all gap columns, but also keeping track of what was removed and
        what was kept. This method keeps track of column names that were actually
        written to file, so that later on column names could be set to the original
        value. This is crucial for a correct merge. 
        (see read_extendend_alignment_and_relabel_columns) 
        '''
        mut_subalg = self.subalignment.get_mutable_alignment()
        remaining_cols = mut_subalg.delete_all_gap()
        mut_subalg.write_to_path(path)
        self.annotations["ref.alignment.columns"] = remaining_cols
        _LOG.debug("subalignment without allgap columns written to %s; %d columns remaining." %(path, len(remaining_cols)))
        return remaining_cols         
        
    def read_extendend_alignment_and_relabel_columns(self, orig_path, extension_path, convert_to_string = True):
        '''
        This method goes with write_subalignment_without_allgap_columns method.
        It enables reading back an alignment that was previously written to disk, 
        and relabeling its columns with the original labels. 
        extension_path is a path to an .sto file (or a list of paths).
        Alignments from these .sto files are also read, and merged with the 
        original (base) alignment.  
        '''
        remaining_cols = self.annotations["ref.alignment.columns"]
        assert remaining_cols is not None and len(remaining_cols) != 0, ("Subproblem"
        " needs to have a proper list of alignment columns associated with it") 
        
        _LOG.debug("Reading %s %s and relabeling it based on %d orig column labels." %(orig_path, extension_path, len(remaining_cols)))
        
        ap_alg = ExtendedAlignment(self.fragments.keys())
        ap_alg.build_extended_alignment(orig_path, extension_path, convert_to_string)
        ap_alg.relabel_original_columns(remaining_cols)
        return ap_alg
                            