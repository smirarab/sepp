'''
Created on Nov 26, 2012

@author: smirarab
'''
import os
import pickle
import sys
import threading
from sepp import get_logger

_LOG = get_logger(__name__)

class CheckPointState(object):
    '''
    The current state as relevant to the checkpointign feature
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        self.root_problem = None
        
def save_checkpoint(checkpoint_manager):
    #assert os.path.exists(self.checkpoint_path)
    if checkpoint_manager.is_checkpointing:
        _LOG.info("Checkpoint is being updated: %s" %str(checkpoint_manager.checkpoint_state))
        currenlimit = sys.getrecursionlimit()
        sys.setrecursionlimit(100000)
        pickle.dump(checkpoint_manager.checkpoint_state, open(checkpoint_manager.checkpoint_path,"w"))
        sys.setrecursionlimit(currenlimit)
        _LOG.info("Checkpoint Saved to: %s" %str(checkpoint_manager.checkpoint_path))
        threading.Timer(60, save_checkpoint, args={checkpoint_manager}).start() 
            
class CheckPointManager:
    
    def __init__(self, checkpoint_file):            
        self.checkpoint_path = checkpoint_file
        self.checkpoint_state = None
        self.checkpoint_path = None
        self.is_recovering = False
        self.is_checkpointing = False
        self._init_state_and_file()

    def _init_state_and_file(self):
        if self.checkpoint_path is None:
            return                    
        if not os.path.exists(self.checkpoint_path):
            open("w").close(self.checkpoint_file)       
            self.checkpoint_state = CheckPointState()
            self.is_checkpointing = True
        else:                            
            self.restore_checkpoint()
            self.is_recovering = True
            self.is_checkpointing = True
                                
    def restore_checkpoint(self):
        assert os.path.exists(self.checkpoint_path)
        self.checkpoint_state = pickle.load(open(self.checkpoint_path))        
    
    def remove_checkpoint_file(self):
        os.remove(self.checkpoint_path)
        
    def start_checkpointing(self, root_problem):
        if self.is_checkpointing:
            self.checkpoint_state.root_problem = root_problem 
            save_checkpoint(self)
#    def backup_temp_directory(self, path):
#        assert(os.path.exists(path))
#        idx = 0
#        while os.path.exists("%s_back-%d" %(path,idx)): idx += 1                                        
#        new_name = "%s_back-%d" %(path,idx)
#        os.rename(path, new_name)
#        return new_name  