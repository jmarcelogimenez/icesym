#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 18:38:26 2018

@author: etekken
"""

import subprocess, os
from PyQt5 import QtCore, QtGui, QtWidgets
from logTabWidget_ui import Ui_LogTabWidget
from utils import remove_values_from_list, show_message, show_errors
from Thread import Thread

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)
    
class LogTabWidget(QtWidgets.QWidget):
    def __init__(self, simulator_dir, case_name, case_dir, current_dir, enable_ppw):
        QtWidgets.QWidget.__init__(self)
        self.ui_tab_widget = Ui_LogTabWidget()
        self.ui_tab_widget.setupUi(self)
        self.set_palette()
        self.simulator_dir = simulator_dir
        self.case_name = case_name
        self.case_dir = case_dir
        self.current_dir = current_dir
        self.enable_ppw = enable_ppw
        self.current_log = ''
        self.lastPos = 0
        self.lastlastPos = -1
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(50)
        self.timer.timeout.connect(self.file_changed)
        self.thread_runSimulation = None
        self.process_killed = False
        return
    
    def set_palette(self):
        """
        Setear el formato del QTextEdit.
        """
        palette = QtGui.QPalette()
        bgc = QtGui.QColor(0, 0, 0)
        palette.setColor(QtGui.QPalette.Base, bgc)
        textc = QtGui.QColor(255, 255, 255)
        palette.setColor(QtGui.QPalette.Text, textc)
        self.setPalette(palette)
        self.ui_tab_widget.logPlainTextEdit.setPalette(palette)
        return
    
    def load_log(self):
        """
        Cargar un log y mostrarlo en el QTextEdit.
        """
        self.ui_tab_widget.logPlainTextEdit.clear()
        dialog = QtWidgets.QFileDialog(self)
        dialog.setWindowTitle('Select Log')
        if dialog.exec_():
            filename = dialog.selectedFiles()[0]
            with  open(filename) as open_log:
                all_text = open_log.readlines()
                all_text = remove_values_from_list(all_text, '\n')
                for line in all_text:
                    self.ui_tab_widget.logPlainTextEdit.appendPlainText(line)
        return

    def file_changed(self):
        if self.current_log == '':
            return

        if self.lastPos == self.lastlastPos:
            return

        with open(self.current_log, 'r') as openFile:
            openFile.seek(self.lastPos)
            newTexto = openFile.read()
            if(len(newTexto)>1):
                self.ui_tab_widget.logPlainTextEdit.appendPlainText(newTexto)
                QtWidgets.QApplication.processEvents()
                self.lastlastPos = self.lastPos
                self.lastPos = self.lastPos + len(newTexto)
        return            
    
    def init_log(self,path):
        self.current_log = path
        self.timer.start()
        return
    
    def reset_log(self):
        self.file_changed()
        self.current_log = ''
        self.lastPos = 0
        self.lastlastPos = -1
        self.timer.stop()
        return
    
    def kill_process(self):
        """
        Verifica si existe una simulacion corriendo y termina el proceso
        """
        if self.thread_runSimulation:
            self.current_pid = self.thread_runSimulation.get_pid() + 2
            if self.current_pid != -1:
                self.process_killed = True
                command = "kill -KILL %s 2> %s/error.log"%(self.current_pid,self.simulator_dir)
                p = subprocess.Popen([command],shell=True)
                p.wait()
                if not show_errors(self.simulator_dir):
                    show_message('Process Killed', 1)
                    QtWidgets.QApplication.processEvents()
                    self.reset_pid()
            else:
                show_message('Cannot find the pid simulation. Please report this bug to the developer.')
        else:
            show_message('There is no simulation running now')
        return
    
    def reset_pid(self):
        self.current_pid = -1
        self.thread_runSimulation = None
        return
    
    def success_simulation(self):
        if not show_errors(self.simulator_dir) and not self.process_killed:
            show_message('Simulation Finished', 1)
        self.process_killed = False
        return
    
    def run_simulation(self):
        """
        Iniciar el simulador de ICESym
        """
        
        # TODO: Verificar si hay componentes desconectados
        simulator = self.simulator_dir+"/exec"
        if not os.path.isfile(simulator):
            show_message('Cannot find exec file. Please, create this file to run the simulator')
            return

        logfile = '%s/run.log'%(self.simulator_dir)
        command = "%s/exec %s %s 1> %s/run.log 2> %s/error.log"%(self.simulator_dir,self.case_dir,self.case_name,self.simulator_dir,self.simulator_dir)
        self.thread_runSimulation = Thread(command)
        # limpio el log y agrego el archivo al watcher
        self.ui_tab_widget.logPlainTextEdit.clear()
        self.thread_runSimulation.started.connect(lambda: self.init_log(logfile))
        self.thread_runSimulation.finished.connect(self.reset_log)
        self.thread_runSimulation.finished.connect(self.reset_pid)
        self.thread_runSimulation.finished.connect(self.enable_ppw)
        self.thread_runSimulation.finished.connect(self.success_simulation)
        # inicio el thread
        self.thread_runSimulation.start()
        return