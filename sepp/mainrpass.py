'''
Created on Jun 8, 2011

@author: smirarab
'''
from sepp.filemgr import temp_fs
import sys
import os
from sepp.scheduler import start_worker
from sepp.settings import settings

if __name__ == '__main__':
    
    tmp_dir = sys.argv[1]
    job_name = sys.argv[2]
    
    if tmp_dir is None:
        tmp_dir = os.path.join(os.path.expanduser('~'), '.rpas')
    subdir = job_name
    tmp_dir = os.path.abspath(os.path.join(tmp_dir, subdir))
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
            
    temporaries_dir = temp_fs.create_top_level_temp(parent=tmp_dir, prefix='temp')
    
    start_worker(settings.get_setting("num.cpus"))