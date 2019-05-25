#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 16:43:33 2018

@author: etekken
"""

import sys, os
from PyQt5 import QtWidgets, QtGui

def main():
    app = QtWidgets.QApplication(sys.argv)
    try:
        app.setStyle(QtWidgets.QStyleFactory.create('Fusion'))
    except:
        app.setStyle(QtWidgets.QStyleFactory.keys()[0])

    if "ICESYM_INST_DIR" not in os.environ.keys():
        # no uso show_message porque esta en utils y ahi se necesita la 
        # variable de entorno que estoy verificando si existe
        font = QtGui.QFont()
        font.setFamily("Ubuntu")
        font.setPointSize(10)
        content = "The enviroment variable of the instalation directory is not set. Please, execute the bashrc file"
        w = QtWidgets.QMessageBox(3,'Error',content,QtWidgets.QMessageBox.Ok)
        w.setFont(font)
        w.exec_()
        sys.exit(1)

    from NewCaseDialog import NewCaseDialog
    open_dialog = NewCaseDialog()
    ret = open_dialog.exec_()
    if ret:
        from ICESymMainWindow import ICESymMainWindow
        ismw = ICESymMainWindow(open_dialog.case_name, open_dialog.case_dir, open_dialog.case_type)
        ismw.show()
        sys.exit(app.exec_())
    return

if __name__ == '__main__':
    main()