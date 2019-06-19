#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 18:26:36 2019

@author: etekken
"""

import os
from PyQt5 import QtCore, QtGui, QtWidgets
from plotTypeTwoWidget_ui import Ui_PlotTypeTwoWidget
from utils import show_message
from numpy import loadtxt, take
from units import UNITS, CONVERSIONS

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
    
LISTNDOFA = ['Density', 'Velocity', 'Pressure']
    
class PlotTypeTwoWidget(QtWidgets.QWidget):
    def __init__(self, plot_function, current_run_dir, run_attributes, current_objects, plot_type, get_oa, set_oa):
        QtWidgets.QWidget.__init__(self)
        self.ui = Ui_PlotTypeTwoWidget()
        self.ui.setupUi(self)
        self.current_run_dir        = current_run_dir
        self.plot_function          = plot_function
        self.plot_type              = plot_type
        self.current_index_element  = -1
        self.change_attributes(run_attributes,current_objects)
        # Funciones para obtener y setear los archivos abiertos, existe una 
        # unica instancia compartida por todos los PlotTypeOneWidget y la tiene
        # el padre postProcessWidget
        self.get_open_archives = get_oa
        self.set_open_archives = set_oa        
        return
    
    def change_attributes(self, run_attributes, current_objects):
        self.run_attributes = run_attributes
        self.current_objects = current_objects
        self.set_rpms()
        self.set_time_restrictions(0)
        self.set_elements()
        self.set_variables()
        return
    
    def set_time_restrictions(self, index_rpm):
        if self.run_attributes['final_times']:
            time = self.run_attributes['final_times'][index_rpm]
            self.ui.time.setText(str(time))
            self.ui.time.setValidator(QtGui.QDoubleValidator(0.0, float(time), 3))
        return

    def set_rpms(self):
        self.ui.rpms.clear()
        for irpm in self.run_attributes['rpms']:
            self.ui.rpms.addItem(str(irpm))
        return
    
    def set_elements(self):
        for index,itube in enumerate(self.current_objects['Tubes']):
            if itube.object['typeSave']==0:
                self.ui.element.addItem(itube.object['label'])
        return
    
    def set_variables(self):
        self.ui.variable.addItems(LISTNDOFA)
        return
    
    def read_normal_txt(self, archive, time, variable_index, unit):
        open_archives = self.get_open_archives()
        # Nodo, ciclo, angulo, tiempo .. variables
        if archive not in open_archives.keys():
            A = loadtxt(archive)
            self.set_open_archives(archive,A)
        else:
            A = open_archives[archive]

        scale = CONVERSIONS[unit]
        times = A[:,3]
        new_time = min(times, key=lambda x:abs(x-time)) # el mas cercano en la lista

        A_time_filtered = A[A[:,3] == new_time]
        data = take(A_time_filtered, [0,4+variable_index], axis=1)
        for idata in data:
            idata[1] = idata[1]*scale
        return (data,new_time)
    
    def get_plot_attributes(self):
        plot_attributes = {}
        plot_attributes['rpm'] = str(self.ui.rpms.currentText())
        plot_attributes['time'] = float(self.ui.time.text())
        element = self.ui.element.currentText()
        plot_attributes['index_element'] = str(element[-1])
        plot_attributes['legend'] = str(self.ui.legend.text())
        plot_attributes['variable_index'] = int(self.ui.variable.currentIndex())
        plot_attributes['variable'] = str(self.ui.variable.currentText())
        plot_attributes['title'] = str(self.ui.title.text())
        plot_attributes['figure_number'] = self.ui.figure_number.currentIndex()-1
        plot_attributes['units'] = str(self.ui.units.currentText())
        return plot_attributes
    
    def change_variable(self, variable):
        if variable == '':
            return
        self.ui.units.clear()
        self.ui.units.addItems(UNITS[variable])
        return

    def prepare_plot(self):
        datas = []
        legends = []
        self.ui.plot_pushButton.setEnabled(False)
        QtWidgets.QApplication.processEvents()
        try:
            plot_attributes = self.get_plot_attributes()
            rpm_folder = os.path.join(self.current_run_dir,"RPM_%s"%plot_attributes['rpm'])
            archive = os.path.join(rpm_folder,"tube_"+plot_attributes['index_element']+".txt")
            (data,new_time) = self.read_normal_txt(archive,plot_attributes['time'],plot_attributes['variable_index'],plot_attributes['units'])
            datas.append(data)
            legends.append(plot_attributes['legend'] + '_time_' + str(new_time) + '_tube_' + plot_attributes['index_element'])
            n_plots = self.plot_function(datas, plot_attributes['title'], legends, 'Nodes', plot_attributes['variable'], '', \
                                         plot_attributes['units'], plot_attributes['figure_number'], 4)
            if plot_attributes['figure_number']==-1:
                self.ui.figure_number.addItem('Figure '+str(n_plots-1))
        except:
            show_message('Error trying to plot the current selection')
        self.ui.plot_pushButton.setEnabled(True)
        return