#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 18:38:26 2018

@author: etekken
"""

import subprocess, os
from sys import platform
from PyQt5 import QtCore, QtGui, QtWidgets
from logTabWidget_ui import Ui_LogTabWidget
from utils import remove_values_from_list, show_message, show_errors, SIMULATOR_PATH, RUNS_PATH
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
    def __init__(self, case_name, case_dir, folder_name, rpms, save_data):
        QtWidgets.QWidget.__init__(self)
        self.ui_tab_widget = Ui_LogTabWidget()
        self.ui_tab_widget.setupUi(self)
        self.set_palette()
        self.case_name          = case_name   # Nombre del caso (permanente)
        self.case_dir           = case_dir    # Directorio donde se ubica el .py del caso (variable)
        self.current_run_dir    = os.path.join(RUNS_PATH,folder_name) # Directorio donde se guarda la corrida (variable)
        self.rpms               = rpms # RPMS a correr (variable)
        self.save_data_f = save_data
        self.current_log = ''
        self.lastPos = 0
        self.lastlastPos = -1
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(50)
        self.timer.timeout.connect(self.file_changed)
        self.thread_runSimulation = None
        self.process_killed = False
        return

    def change_attributes(self, folder_name, rpms):
        """
        On saving the case, change the attributes used to simulate (important
        if any of them has changed)
        """
        self.rpms = rpms
        self.current_run_dir = os.path.join(RUNS_PATH,folder_name)
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
                errorlog_path = os.path.join(SIMULATOR_PATH,"error.log")
                if 'linux' in platform:
                    command = "kill -KILL %s 2> %s"%(self.current_pid,errorlog_path)
                elif 'win' in platform:
                    command_kill = "taskkill /F /IM exec.exe"
                    execwfile = os.path.join(SIMULATOR_PATH,"execw_k.bat")
                    command = "echo %s > %s"%(command_kill,execwfile)
                    os.system(command)
                    command = os.path.join(SIMULATOR_PATH,"execw_k.bat")
                p = subprocess.Popen([command],shell=True)
                p.wait()
                if not show_errors(SIMULATOR_PATH):
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
        if not show_errors(SIMULATOR_PATH) and not self.process_killed:
            show_message('Simulation Finished', 1)
        self.process_killed = False
        return
    
    def show_info(self):
        """
        Shows the current rpm folder located and the ones to calculate
        """
        # Folders located
        msg = 'RPM(s) already simulated: \n'
        if os.path.isdir(self.current_run_dir):
            rpm_folders = [f.replace('RPM_','') for f in os.listdir(self.current_run_dir) if 'RPM_' in f]
            for irpm in rpm_folders:
                msg += '- %s \n'%irpm
        else:
            msg += '- None \n'
        msg += '\n RPM(s) to simulate: \n'
        for irpm in self.rpms:
            msg += '- %s \n'%irpm
        show_message(msg,1)
        return

    def run_simulation(self):
        """
        Iniciar el simulador de ICESym
        """

        # TODO: Verificar si hay componentes desconectados
        
        simulator = None
        if 'linux' in platform:
            simulator = os.path.join(SIMULATOR_PATH,"exec")
        elif 'win' in platform:
            simulator = os.path.join(SIMULATOR_PATH,"exec.exe")
        if not simulator or not os.path.isfile(simulator):
            show_message('Cannot find exec file. Please, create this file to run the simulator')
            return

        # Salvar el caso
        self.save_data_f(None,True)

        logfile = os.path.join(SIMULATOR_PATH,'run.log')
        errorlogfile = os.path.join(SIMULATOR_PATH,'error.log')
        if 'linux' in platform:
            command = "%s %s %s 1> %s 2> %s"%(simulator,self.case_dir,self.case_name,logfile,errorlogfile)
        elif 'win' in platform:            
            command = "echo.> %s"%logfile
            os.system(command)
            command_ex = "%s %s %s 1^> %s 2^> %s"%(simulator,self.case_dir,self.case_name,logfile,errorlogfile)
            execwfile = os.path.join(SIMULATOR_PATH,'execw.bat')
            command = "echo %s > %s"%(command_ex,execwfile)
            os.system(command)
            command = os.path.join(SIMULATOR_PATH,'execw.bat')

        self.thread_runSimulation = Thread(command)
        # limpio el log y agrego el archivo al watcher
        self.ui_tab_widget.logPlainTextEdit.clear()
        self.thread_runSimulation.started.connect(lambda: self.init_log(logfile))
        self.thread_runSimulation.started.connect(self.show_info)
        self.thread_runSimulation.finished.connect(self.reset_log)
        self.thread_runSimulation.finished.connect(self.reset_pid)
        self.thread_runSimulation.finished.connect(self.success_simulation)
        # inicio el thread
        self.thread_runSimulation.start()
        return