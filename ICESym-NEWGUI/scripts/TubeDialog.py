#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 00:42:08 2018

@author: etekken
"""

import os
import json
import copy
import numpy as np
from time import localtime, strftime
from PyQt5 import QtGui, QtWidgets, QtCore
from tubeDialog_ui import Ui_TubeDialog
from utils import check_if_float, show_message, convert_string, INSTALL_PATH,\
                  LOADS_PATH, save_current_configuration_aux, load_configuration_aux

JSON_TUBE_KEYS = ['diameter', 'nnod', 'label', 'twall', 'ndof', 'state_ini', \
                   'histo', 'xnod', 'longitud', "tleft", "nleft", "tright", \
                   "nright", 'typeSave', 'numNorm', 'posNorm', 'histo']

def configure_default_tube(object_tube, nstroke = -1):
    """
    input: a dictionary corresponding to a tube
    Function to create the arrays of values depending on the number of nodes
    (nnod) in the tube on default input
    """
    nnod = object_tube['nnod']
    for inn in range(0,nnod-1):
        object_tube['diameter'].append(object_tube['diameter'][0])
        object_tube['twall'].append(object_tube['twall'][0])
        object_tube['state_ini'].append(object_tube['state_ini'][0])
    vals = np.linspace(0, object_tube['longitud'], nnod)
    for iv in vals:
        object_tube['xnod'].append(iv)
    return object_tube

class TubeDialog(QtWidgets.QDialog):
    """
    class to manage the valve atributes. If current_valve is None, we are creating
    a new one. On the other hand, we are modifying an old one.
    """
    def __init__(self, current_tube = None, item_index = 0):
        QtWidgets.QDialog.__init__(self)
        self.ui_td = Ui_TubeDialog()
        self.ui_td.setupUi(self)
        self.setBaseSize(560, 725)        
        self.current_dict = None # default valve dictionary
        self.setWindowTitle( self.windowTitle() + " " + str(item_index) )
        if not current_tube:
            self.load_default()
        else:
            self.current_dict = current_tube
            self.set_parameters()
        self.set_restrictions()
        self.change_typeSave(current_tube['typeSave'])
        return

    def set_restrictions(self):
        """
        set the initial restrictions on lineEdits
        """
        self.ui_td.length_lineEdit.setValidator(QtGui.QDoubleValidator(0, 1000, 3))
        self.ui_td.numNorm.setMaximum(self.current_dict['nnod'])
        return
    
    def load_default(self):
        """
        load the default atmosphere template
        """
        filename = os.path.join(INSTALL_PATH,"templates","tube_default.json")
        if not os.path.isfile(filename):
            show_message("Cannot find default tube configuration file")
            return
        with open(filename) as openedfile:
            try:
                self.current_dict = json.load(openedfile)
            except ValueError as error:
                show_message('JSON object issue: %s'%error)
                return
        if self.current_dict is None:
            show_message("An error ocurred with the json default tube archive")
            return
        if not self.check_json_keys():
            return
        self.set_parameters()
        return

    def save_current_configuration(self):
        dict_to_save = copy.deepcopy(self.current_dict)
        dict_to_save['nleft']   = -1
        dict_to_save['nright']  = -1
        dict_to_save['tleft']   = "<none>"
        dict_to_save['tright']  = "<none>"
        time_label = strftime("%Y%m%d", localtime())
        dict_to_save['label']   = "%s_save_%s"%(dict_to_save['label'],time_label)
        save_current_configuration_aux(self,dict_to_save)
        return

    def load_configuration(self):
        (success,new_configuration) = load_configuration_aux(self,'tube')
        if not success or not self.check_json_keys(new_configuration):
            return
        self.current_dict = new_configuration
        self.set_parameters()
        return

    def check_json_keys(self, configuration_to_check = None):
        """
        check if the loaded json include all the obligatory keys
        """
        if not configuration_to_check:
            configuration_to_check = self.current_dict
        for ikey in JSON_TUBE_KEYS:
            if ikey not in configuration_to_check.keys():
                show_message("Wrong number of keys in json default tube archive")
                return False
        return True

    def set_parameters(self):
        """
        puts the loaded tube parameters in the items
        """
        if not self.current_dict:
            return
        self.ui_td.name_lineEdit.setText(self.current_dict['label'])
        self.ui_td.nnodes_spinBox.setValue(self.current_dict['nnod'])
        self.ui_td.freedeg_spinBox.setValue(self.current_dict['ndof'])
        self.ui_td.numNorm.setValue(self.current_dict['numNorm'])
        self.ui_td.typeSave.setCurrentIndex(self.current_dict['typeSave'])
        self.ui_td.length_lineEdit.setText(str(self.current_dict['longitud']*1e3))
        if 'curvature' in self.current_dict.keys():
            if self.current_dict['curvature'] != []:
                self.ui_td.has_curvature_checkBox.setChecked(True)
                self.set_table('curvature')
        if 'xnod' in self.current_dict.keys():
            self.set_table('xnod')
        if 'diameter' in self.current_dict.keys():
            self.set_table('diameter')
        if 'twall' in self.current_dict.keys():
            self.set_table('twall')
        if 'state_ini' in self.current_dict.keys():
            self.set_table('state_ini')
        if 'histo' in self.current_dict.keys():
            self.set_table('histo')
        if 'posNorm' in self.current_dict.keys():
            self.set_table('posNorm')
        self.check_curvature()
        return

    def change_nnodes(self):
        tables = [self.ui_td.nodalcoord_tableWidget, self.ui_td.diameter_tableWidget, \
                  self.ui_td.temperature_tableWidget, self.ui_td.curvature_tableWidget, \
                  self.ui_td.state_tableWidget]
        nrows = self.ui_td.nnodes_spinBox.value()
        for itable in tables:
            itable.clearContents()
            itable.setRowCount(0)
            for irow in range(0, nrows):
                it = QtWidgets.QTableWidgetItem()
                itable.setItem(irow, 0, it)
                itable.setRowCount(irow + 1)
        self.set_equispaced()
        return

    def set_equispaced(self):
        if not self.ui_td.equispaced_checkBox.isChecked() or self.ui_td.length_lineEdit.text() == '':
               return
        try:
            longitud = float ( self.ui_td.length_lineEdit.text() )
            nnodes = self.ui_td.nnodes_spinBox.value()
        except:
            show_message("Bad Value in Length or Number of Nodes")
            return
        self.ui_td.nodalcoord_tableWidget.clearContents()
        nrows = self.ui_td.nodalcoord_tableWidget.rowCount()
        vals = np.linspace(0, longitud, nnodes)
        for irow in range(0, nrows):
            it = QtWidgets.QTableWidgetItem(str(vals[irow]))
            it.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable) 
            self.ui_td.nodalcoord_tableWidget.setItem(irow, 0, it)
        return
    
    def input_dialog(self, dialog_label):
        from inputDialog_ui import Ui_InputDialog
        cid = Ui_InputDialog()
        dialog = QtWidgets.QDialog()
        dialog.setFixedSize(327, 80)
        cid.setupUi(dialog)
        cid.constant_lineEdit.setValidator( QtGui.QDoubleValidator(0, 1000, 3) )
        cid.dialog_label.setText(dialog_label)
        ret = dialog.exec_()
        if not ret:
            return -1
        return cid.constant_lineEdit.text()
    
    
    def set_linear_table(self):
        """
        set a values defined with linear function (min/max user input) 
        in all the rows of the table
        """
        if self.sender()==self.ui_td.linear_diameter_pushButton:
            table = self.ui_td.diameter_tableWidget
            dialog_label1 = 'Insert First Diameter [mm]:'
            dialog_label2 = 'Insert Second Diameter [mm]:'
        elif self.sender()==self.ui_td.linear_temperature_pushButton:
            table = self.ui_td.temperature_tableWidget
            dialog_label1 = 'Insert First Temperature [K]:'
            dialog_label2 = 'Insert Second Temperature [K]:'
        else:
            return
        
        vals = []
        for i in range(0,2):
            textval = self.input_dialog(dialog_label1 if i==0 else dialog_label2)
            if textval == -1:
                return
            try:
                vals.append( float( textval ) )
            except:
                show_message("Wrong input value")
                return
        nrows = self.ui_td.nodalcoord_tableWidget.rowCount()
        vals = np.linspace( vals[0], vals[1], nrows )
        table.clearContents()
        for irow in range(0, nrows):
            it = QtWidgets.QTableWidgetItem(str(vals[irow]))
            it.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable) 
            table.setItem(irow, 0, it)
        return
    
    def load_table(self):
        """
        loads a defined archive with values for the
        angle-lift or the lift-dc tables
        """
        key = ''
        if self.sender().objectName()=="load_nodalcoord_pushButton":
            key = 'xnod'
            shape = 1
        elif self.sender().objectName()=="load_temperature_pushButton":
            key = 'twall'
            shape = 1
        elif self.sender().objectName()=="load_curvature_pushButton":
            key = 'curvature'
            shape = 1
        elif self.sender().objectName()=="load_diameter_pushButton":
            key = 'diameter'
            shape = 1
        elif self.sender().objectName()=="load_state_pushButton":
            key = 'state_ini'
            shape = 2
        else:
            return
            
        assert(key!='')

        dialog = QtWidgets.QFileDialog(self)
        dialog.setNameFilter(" (*.txt)")
        dialog.setWindowTitle('Select a File to Open')
        dialog.setDirectory(LOADS_PATH)
        if dialog.exec_():
            filename = dialog.selectedFiles()[0]
        else:
            show_message("Nothing selected", 2)
            return

        if key in self.current_dict.keys():
            msg = "You will lose the current values. Do you want to continue?"
            reply = show_message(msg,4,QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No)
            if reply == QtWidgets.QMessageBox.Yes:
                self.current_dict[key] = np.loadtxt(filename)
            else:
                show_message("Operation cancelled", 1)
                return
        else:
            self.current_dict[key] = np.loadtxt(filename)

        nrows = self.ui_td.nnodes_spinBox.value()
        results = (len(self.current_dict[key].shape) == shape, self.current_dict[key].shape[0] == nrows)
        if all(results):
            pass
        else:
            if not results[0]:
                error_msg = "Wrong type of data"
            elif not results[1]:
                archi_length = self.current_dict[key].shape[0] 
                error_msg = "Number of nodes don't match with archive length (%s)"%archi_length
            show_message(error_msg)
            self.current_dict.pop(key, None)
            return
        self.set_table(key)
        return
    
    def set_constant_table(self):
        """
        set a constant value (user input) in all the rows of the table
        """
        dialog_labels = []
        if self.sender()==self.ui_td.const_diameter_pushButton:
            table = self.ui_td.diameter_tableWidget
            dialog_labels.append('Insert Diameter [mm]:')
        elif self.sender()==self.ui_td.constant_temperature_pushButton:
            table = self.ui_td.temperature_tableWidget
            dialog_labels.append('Insert Temperature [K]:')
        elif self.sender()==self.ui_td.constant_state_pushButton:
            table = self.ui_td.state_tableWidget
            dialog_labels.append('Insert Density [kg/m^3]:')
            dialog_labels.append('Insert Velocity [m/s]:')
            dialog_labels.append('Insert Pressure [Pa]:')            
        else:
            return
        
        constvals = []
        for ilabel in dialog_labels:
            textval = self.input_dialog(ilabel)
            if textval == -1:
                return
            try:
                constvals.append( float( textval ) )
            except:
                show_message("Wrong input value")
                return
        nrows = self.ui_td.nodalcoord_tableWidget.rowCount()
        ncols = len(dialog_labels)
        vals = []
        for icol in range(0, ncols):
            vals.append( np.linspace( constvals[icol], constvals[icol], nrows ) )
        table.clearContents()
        for icol in range(0, ncols):
            for irow in range(0, nrows):
                it = QtWidgets.QTableWidgetItem(str(vals[icol][irow]))
                it.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable) 
                table.setItem(irow, icol, it)
        return
    
    def check_curvature(self):
        """
        if the checkBox is enabled, this function enables the table
        """
        enabled_curv = self.ui_td.has_curvature_checkBox.isChecked()
        self.ui_td.curvature_tableWidget.setEnabled(enabled_curv)
        return
        

    def set_table(self, key):
        """
        set the values of the 5 available tables. <table> is the current table widget
        to set, <ncols> the number of columns to set and <factor> a unit conversion 
        parameter
        """
        table = None
        if key=='xnod':
            table = self.ui_td.nodalcoord_tableWidget
            ncols = 1
            factor = 1e3
            flags = QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable
        elif key=='diameter':
            table = self.ui_td.diameter_tableWidget
            ncols = 1
            factor = 1e3
            flags = QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable
        elif key=='twall':
            table = self.ui_td.temperature_tableWidget
            ncols = 1
            factor = 1
            flags = QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable
        elif key=='state_ini':
            table = self.ui_td.state_tableWidget
            ncols = 3
            factor = 1
            flags = QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable
        elif key=='curvature':
            table = self.ui_td.curvature_tableWidget
            ncols = 1
            factor = 1
            flags = QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable
        elif key=='posNorm':
            table = self.ui_td.posNorm
            ncols = 1
            factor = 1
            flags = QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable
        elif key=='histo':
            table = self.ui_td.histo
            ncols = 1
            factor = 1
            flags = QtCore.Qt.ItemIsEnabled

        assert(table!=None)

        table.clearContents()
        table.setRowCount(0)
        for ituple in self.current_dict[key]:
            current_row = table.rowCount()
            table.setRowCount(current_row + 1)
            for icol in range(0, ncols):
                val = -1
                # tengo que verificar si es una lista, array de numpy, etc, 
                # porque sino no funciona el indexado
                try:
                    val = ituple[icol]
                except:
                    val = ituple
                it = QtWidgets.QTableWidgetItem(str(val*factor))
                it.setFlags(flags)
                table.setItem(current_row, icol, it)
        return
    
    def change_numNorm_max(self, nnod):
        self.ui_td.numNorm.setMaximum(nnod)
        self.ui_td.posNorm.clear()
        self.ui_td.histo.clear()
        return
    
    def change_numNorm(self, numNorm):
        pre_row_count = self.ui_td.posNorm.rowCount()        
        self.ui_td.posNorm.setRowCount(numNorm)
        self.ui_td.histo.setRowCount(numNorm)
        for irow in range(pre_row_count,numNorm):
            it = QtWidgets.QTableWidgetItem('0.0')
            it.setFlags(QtCore.Qt.ItemIsEnabled)
            self.ui_td.histo.setItem(irow,0,it)
            it2 = QtWidgets.QTableWidgetItem('0.0')
            it2.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable)
            self.ui_td.posNorm.setItem(irow,0,it2)
        return
    
    def transform_normPoints(self, item):
        if not check_if_float(item.text(),'Normalized Points'):
            show_message('Please insert a valid value')
            item.setText('0.0')
            return
        number = float(item.text())
        if number>1.0 or number<0.0:
            show_message('Please insert a value in the range [0,1]')
            item.setText('0.0')
            return
        nrows = self.ui_td.posNorm.rowCount()
        nnodes = self.ui_td.nnodes_spinBox.value()-1
        for irow in range(0, nrows):
            it_posNorm  = self.ui_td.posNorm.item(irow,0)
            it_histo    = self.ui_td.histo.item(irow,0)                
            if it_posNorm and it_histo:
                it_histo.setText(str(int(float(it_posNorm.text())*nnodes)))
        return
    
    def get_table(self, ikey, table, current_dict, factor):
        if not table.isEnabled():
            return
        ncols = table.columnCount()
        nrows = table.rowCount()
        current_dict[ikey] = []
        for irow in range(0, nrows):
            dict_row = []
            for icol in range(0, ncols):
                if not table.item(irow,icol):
                    current_dict[ikey] = []
                    show_message('Please provide all the fields for the %s table'%ikey)
                    # dejo que salte el error en el try de accept
                text = table.item(irow,icol).text()
                text_parsed = convert_string(text)
                text_parsed = text_parsed*factor
                if ncols==1:
                    current_dict[ikey].append(text_parsed)
                else:
                    dict_row.append(text_parsed)
            if dict_row!=[]:
                current_dict[ikey].append(dict_row)
        return
    
    def change_typeSave(self, index):
        state = True if index==1 else False
        self.ui_td.numNorm.setEnabled(state)
        self.ui_td.posNorm.setEnabled(state)
        self.ui_td.histo.setEnabled(state)
        if index==0:
            self.current_dict['histo'] = []
            for inode in range(self.current_dict['nnod']):
                self.current_dict['histo'].append(inode)
        return

    def accept(self):
        # verificaciones        
        if not check_if_float(self.ui_td.length_lineEdit.text(), 'length'):
            return

        success = False
        # argumentos salvados
        self.current_dict['label']      = self.ui_td.name_lineEdit.text()
        self.current_dict['nnod']       = self.ui_td.nnodes_spinBox.value()
        self.current_dict['ndof']       = self.ui_td.freedeg_spinBox.value()
        self.current_dict['longitud']   = float(self.ui_td.length_lineEdit.text())*1e-3
        self.current_dict['typeSave']   = int(self.ui_td.typeSave.currentIndex())
        try:
            self.get_table('diameter', self.ui_td.diameter_tableWidget, self.current_dict, 1e-3)
            self.get_table('xnod', self.ui_td.nodalcoord_tableWidget, self.current_dict, 1e-3)
            self.get_table('twall', self.ui_td.temperature_tableWidget, self.current_dict, 1.0)
            self.get_table('state_ini', self.ui_td.state_tableWidget, self.current_dict, 1.0)
            if self.ui_td.has_curvature_checkBox.isChecked():
                self.get_table('curvature', self.ui_td.curvature_tableWidget, self.current_dict, 1.0)
            else:
                self.current_dict['curvature'] = []
            if self.ui_td.typeSave.currentIndex()==1:
                self.get_table('histo', self.ui_td.histo, self.current_dict, 1)
                self.get_table('posNorm', self.ui_td.posNorm, self.current_dict, 1.0)
                self.current_dict['numNorm'] = self.ui_td.numNorm.value()
            success = True
        except:
            success = False
        if success:
            self.close()
        return

    def cancel(self):
        self.close()
        return