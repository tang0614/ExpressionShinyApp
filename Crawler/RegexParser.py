'''
Created on Nov 22, 2019

@author: Shan.Jiang
'''
import re
from datetime import datetime

class RegexParser:
    def __init__(self,flags = re.S|re.I, log = open("ParserLog.txt",'a'), task_id = "NA"): # @UndefinedVariable
        self.flags      = flags
        self.log        = log
        self.task_id    = task_id
        
        print("\n=================="+str(datetime.now().strftime("%Y-%M-%D %H:%M:%S"))+"==========================\n", file = self.log)
        
    def parse(self, regex_str, source):
        matches = re.compile(regex_str,self.flags).findall(source)
        #if(len(matches) == 0):
            #print("Task ["+str(self.task_id)+"]"+regex_str+" failed @ "+source,file = self.log)
        return matches

        
        
    def set_flags(self,flags):
        self.flags = flags