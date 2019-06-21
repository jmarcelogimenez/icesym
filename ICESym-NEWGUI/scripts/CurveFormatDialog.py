#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 18:07:31 2019

@author: etekken
"""

from PyQt5 import QtWidgets
from curveFormatDialog_ui import Ui_CurveFormatDialog

class CurveFormatDialog(QtWidgets.QDialog):
    """
    """
    def __init__(self):
        QtWidgets.QDialog.__init__(self)
        self.ui_cfd = Ui_CurveFormatDialog()
        self.ui_cfd.setupUi(self)
        return