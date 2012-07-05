'''
Created on Jun 8, 2011

@author: smirarab
'''

class Settings(object):
    '''
    classdocs
    '''    

    def __init__(self):
       self.settings = {
                        "hmmbuild":{"path":"/projects/sate7/tools/bin/hmmbuild-32"},
                        "hmmalign":{"path":"/projects/sate7/tools/bin/hmmalign-32"},
                        "num.cpus":1,
                        "delete_temps":False}
        
    def get_setting(self, attribute):
        return self.settings.get(attribute)
    
settings=Settings()