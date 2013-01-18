'''
Created on Nov 26, 2012

@author: smirarab
'''
import os
import cPickle
import sys
import threading
from sepp import get_logger
import stat
from sepp.filemgr import get_temp_file, get_root_temp_dir, set_root_temp_dir
import gzip
import time
from sepp.config import options
import datetime
#from threading import Lock

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
        self.temp_root = None
        self.cumulative_time = 0
        
def save_checkpoint(checkpoint_manager):    
    '''
    This is the callback function that is called periodically to save the
    current state of the system.  
    '''
    # Note: this module is not bullet proof in terms of race conditions. 
    # Most importantly, it is possible (though extremely unlikely) that 
    # while the new temp path is being written (f.write...) 
    if checkpoint_manager.is_checkpointing:
        #checkpoint_manager.lock.acquire()
        checkpoint_manager.saving = True
        newTmpDest = get_temp_file("dump", "checkpoints")
        _LOG.info("Checkpoint is being updated: %s" %newTmpDest)
        oldTmpFile = open(checkpoint_manager.checkpoint_path).readlines()
        oldTmpFile = None if len(oldTmpFile) == 0 else oldTmpFile[-1].split(",")[0]

        checkpoint_manager.update_time()

        currenlimit = sys.getrecursionlimit()
        sys.setrecursionlimit(100000)
        picklefile = gzip.GzipFile(newTmpDest, 'wb')
        cPickle.dump(checkpoint_manager.checkpoint_state, picklefile, 2)    
        picklefile.close()  
        sys.setrecursionlimit(currenlimit)

        f = open(checkpoint_manager.checkpoint_path,"a")
        f.write("%s, %s\n"%(newTmpDest,datetime.datetime.now()))
        f.close()
        if oldTmpFile is not None:
            os.remove(oldTmpFile)        
        _LOG.info("Checkpoint Saved to: %s and linked in %s." %(newTmpDest,checkpoint_manager.checkpoint_path))
        checkpoint_manager.saving = False      
        #checkpoint_manager.lock.release()
        checkpoint_manager.timer = threading.Timer(options().checkpoint_interval, save_checkpoint, args=[checkpoint_manager])
        checkpoint_manager.timer.setDaemon(True)
        checkpoint_manager.timer.start() 
            
class CheckPointManager:
    
    def __init__(self, checkpoint_file):            
        self.checkpoint_path = checkpoint_file
        self.checkpoint_state = None
        self.is_recovering = False
        self.is_checkpointing = False
        self.timer = None
        self._init_state_and_file()
        self.last_checkpoint_time = time.time()
        self.saving = False
        #self.lock = Lock()
          
    def _init_state_and_file(self):
        if self.checkpoint_path is None:
            return                    
        if (not os.path.exists(self.checkpoint_path) ) or (os.stat(self.checkpoint_path)[stat.ST_SIZE] == 0):
            open(self.checkpoint_path,"w").close()       
            self.checkpoint_state = CheckPointState()
            self.is_checkpointing = True
        else:                            
            self.is_recovering = True
            self.is_checkpointing = True
                                
    def restore_checkpoint(self):
        _LOG.info("Checkpoint is being restored: %s" %str(self.checkpoint_path))
        assert os.path.exists(self.checkpoint_path)
        lastPath = open(self.checkpoint_path).readlines()[-1].split(",")[0]
        picklefile = gzip.GzipFile(lastPath, 'rb')
        self.checkpoint_state = cPickle.load(picklefile)
        picklefile.close()
        set_root_temp_dir(self.checkpoint_state.temp_root)
        _LOG.info("Checkpoint restore finished: %s" %lastPath)
    
    def remove_checkpoint_file(self):
        if self.checkpoint_path is not None and os.path.exists(self.checkpoint_path):
            os.remove(self.checkpoint_path)
        
    def start_checkpointing(self, root_problem):
        if self.is_checkpointing:
            _LOG.info("Checkpoint every %d seconds" %options().checkpoint_interval)
            self.checkpoint_state.root_problem = root_problem 
            self.checkpoint_state.temp_root = get_root_temp_dir()
            if self.checkpoint_state.cumulative_time is None:
                self.checkpoint_state.cumulative_time = 0
            save_checkpoint(self)
            
    def stop_checkpointing(self):
        if self.is_checkpointing:
            self.is_checkpointing = False
            if self.timer is not None:
                self.timer.cancel()
            self.remove_checkpoint_file()
            self.update_time()
            _LOG.info("Stop Checkpointing. Cumulative time: %d" %self.checkpoint_state.cumulative_time)
     
    def pause_checkpointing(self):
        if self.is_checkpointing:           
            if self.timer is not None:
                if self.saving:
                    self.timer.join()
                    self.timer.cancel()
                else:
                    self.timer.cancel()
            _LOG.info("Paused Checkpointing. Cumulative time: %d" %self.checkpoint_state.cumulative_time)
     
    def update_time(self):
        if self.is_checkpointing:
            self.checkpoint_state.cumulative_time += (time.time() - self.last_checkpoint_time)
        self.last_checkpoint_time = time.time()     
        
        
    def get_total_time(self):
        if self.checkpoint_state is not None:
            return self.checkpoint_state.cumulative_time + time.time() - self.last_checkpoint_time
        else:
            return time.time() - self.last_checkpoint_time
        
#    def force_checkpoint(self):
#        if self.is_checkpointing:
#            self.timer = threading.Timer(0, save_checkpoint, args=[self,False])
#            self.timer.setDaemon(True)
#            self.timer.start()  
#            self.timer.join()
