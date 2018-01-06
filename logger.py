#!/usr/bin/env python

# comprehensive logging all stdout to a file

import sys

class Logger(object):
    def __init__(self,filename):
        self.filename = filename
        out = filename.rsplit(".",1)[0]
        self.terminal = sys.stdout
        self.log = open(str(out)+"-log.txt", "a")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass 
