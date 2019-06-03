#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 16:47:03 2019

@author: etekken
"""

import os
from sys import platform
from PyQt5 import QtWidgets, QtGui
from usageDialog_ui import Ui_UsageDialog
from utils import DOCS_PATH,INSTALL_PATH

class UsageDialog(QtWidgets.QDialog):
    def __init__(self):
        QtWidgets.QDialog.__init__(self)
        self.ui_ud = Ui_UsageDialog()
        self.ui_ud.setupUi(self)
        self.current_dict = None
        self.setWindowTitle('ICESym Usage')
        self.setBaseSize(760,800)        
        archive = os.path.join(DOCS_PATH,'About.html')
        self.spacer = '\\' if 'win' in platform else '/'
        self.set_textEdit(archive)
        return
    
    def updateCaseSetup(self, QTreeWidgetItem, index):
        if not QTreeWidgetItem:
            return
        item = QTreeWidgetItem.text(0)
        archive = os.path.join(DOCS_PATH,'%s.html'%item)
        self.set_textEdit(archive)
        return
    
    def set_textEdit(self,archive):
        text = ''
        if os.path.isfile(archive):
            with open(archive) as f:
                for line in f:
                    if 'ICESYM_IMG_DIR' in line:
                        line = line.replace('ICESYM_IMG_DIR',os.path.join(INSTALL_PATH,'images','about')+self.spacer)
                    text += line
                    
                        

        if text!='':
            self.ui_ud.usage_textEdit.clear()
            self.ui_ud.usage_textEdit.append(text)
        self.ui_ud.usage_textEdit.moveCursor(QtGui.QTextCursor.Start)
        return
        