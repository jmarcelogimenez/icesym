#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 30 15:48:21 2019

@author: etekken
"""

import inspect
from utils import show_message

CURRENT_EXCEPTION = False

def init_exception_handling():
    global CURRENT_EXCEPTION, LOWER_EXCEPTION
    CURRENT_EXCEPTION = False
    LOWER_EXCEPTION = False
    return

def handle_exception(message):
    global CURRENT_EXCEPTION
    global LOWER_EXCEPTION
    parent_f = inspect.getouterframes( inspect.currentframe() )[2][3]
#    import time
#    print 'exception at %s %s'%(time.time(), parent_f)
    
    if not LOWER_EXCEPTION:
        show_message(message)
        LOWER_EXCEPTION = True
    if parent_f == 'main':
        show_message(message)
        CURRENT_EXCEPTION = False
        LOWER_EXCEPTION = False
    else:
        CURRENT_EXCEPTION = True
    return