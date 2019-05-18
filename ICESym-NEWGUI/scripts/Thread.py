# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 18:01:53 2016

@author: santiago
"""

from PyQt5.QtCore import QThread
import subprocess

class Thread(QThread):
    def __init__(self, cmd):
        QThread.__init__(self)
        self.cmd = cmd
        self.pid = -1
        return
        
    def __del__(self):
        self.wait()
        return
    
    def get_pid(self):
        return self.pid

    def run(self):
        p=subprocess.Popen([self.cmd],shell=True)
        self.pid = p.pid
        p.wait()
        return
        