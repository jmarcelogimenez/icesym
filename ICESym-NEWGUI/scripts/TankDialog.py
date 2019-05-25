#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 16:35:51 2018

@author: etekken
"""

import os, json
from PyQt5 import QtWidgets, QtCore, QtGui
from tankDialog_ui import Ui_TankDialog
from utils import M2MM, SQM2SQCM, CM2CCM, CCM2CM, MM2M, SQCM2SQM, check_if_float,\
                    show_message, convert_string, INSTALL_PATH, LOADS_PATH
import numpy as np

JSON_TANK_KEYS = ['label','nnod','ndof','Volume','mass','h_film','Area_wall','T_wall','Cd_ports','state_ini',\
                  'exh2tube','int2tube','histo', 'extras','state_ini_0']

PARSED = {}
PARSED['Volume']    = [['Volume',CM2CCM],[]]
PARSED['mass']      = [['mass',M2MM],[]]
PARSED['Area_wall'] = [['Area_wall',SQM2SQCM],[]]

EXTRAS = {}

class TankDialog(QtWidgets.QDialog):
    """
    class to manage the tank atributes. If current_tank is None, we are 
    creating a new one. On the other hand, we are modifying an old one.
    """
    def __init__(self, current_tank = None, item_index = 0):
        QtWidgets.QDialog.__init__(self)
        self.ui_td = Ui_TankDialog()
        self.ui_td.setupUi(self)
        self.setFixedSize(452, 578)
        self.set_restrictions()
        self.current_dict = None # default tank dictionary
        self.setWindowTitle( self.windowTitle() + " " + str(item_index) )
        if not current_tank:
            self.load_default()
        else:
            self.current_dict = current_tank
        self.set_tubes_connections()
        self.set_parameters()
        return
    
    def set_restrictions(self):
        """
        set the initial restrictions on lineEdits
        """
        self.ui_td.Volume.setValidator(QtGui.QDoubleValidator(0, 100000000, 3))
        self.ui_td.mass.setValidator(QtGui.QDoubleValidator(0, 100000000, 3))
        self.ui_td.h_film.setValidator(QtGui.QDoubleValidator(0, 100000000, 3))
        self.ui_td.Area_wall.setValidator(QtGui.QDoubleValidator(0, 100000000, 3))
        self.ui_td.T_wall.setValidator(QtGui.QDoubleValidator(0, 100000000, 3))
        return
    
    def load_default(self):
        """
        load the default tank template
        """
        filename = os.path.join(INSTALL_PATH,"templates","tank_default.json")
        if not os.path.isfile(filename):
            show_message("Cannot find default tank configuration file")
            return
        with open(filename) as openedfile:
            try:
                self.current_dict = json.load(openedfile)
            except ValueError as error:
                show_message('JSON object issue: %s'%error)
                return
        if self.current_dict is None:
            show_message("An error ocurred with the json default tank archive")
            return
        if not self.check_json_keys():
            self.close()
            return
        self.set_parameters()
        return

    def check_json_keys(self):
        """
        check if the loaded json include all the obligatory keys
        """
        for ikey in JSON_TANK_KEYS:
            if ikey not in self.current_dict.keys():
                show_message("Wrong number of keys in json default tank archive")
                return False
        return True
    
    def set_tubes_connections(self):
        self.ui_td.histo.setSpacing(3)
        for itype_tube in ('int2tube','exh2tube'):
            for itube in self.current_dict[itype_tube]:
                name = "Connection with tube %s"%itube
                it = QtWidgets.QListWidgetItem(name)
                it.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsUserCheckable)
                it.setCheckState(0)
                self.ui_td.histo.addItem(it)
        return
    
    def load_table(self):
        """
        loads a defined archive with values for the
        Cd_ports or the state_ini tables
        """
        key = ''
        ncols = -1
        if self.sender().objectName()=="load_Cd_ports_pushButton":
            key = 'Cd_ports'
            ncols = 1
            off = 1
            table = self.ui_td.Cd_ports
        elif self.sender().objectName()=="load_state_ini_pushButton":
            key = 'state_ini'
            ncols = 3
            off = 0
            table = self.ui_td.state_ini
            
        assert(key!='' and ncols>0)

        dialog = QtWidgets.QFileDialog(self)
        dialog.setNameFilter(" (*.txt)")
        dialog.setWindowTitle('Select a File to Open')
        dialog.setDirectory(LOADS_PATH)
        if dialog.exec_():
            filename = dialog.selectedFiles()[0]
        else:
            show_message("Nothing selected", 2)
            return

        if key in self.current_dict.keys() and self.current_dict[key]!=[]:
            msg = "You will lose the current values. Do you want to continue?"
            reply = show_message(msg,4,QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No)
            if reply == QtWidgets.QMessageBox.Yes:
                self.current_dict[key] = np.loadtxt(filename)
            else:
                show_message("Operation cancelled", 1)
                return
        else:
            self.current_dict[key] = np.loadtxt(filename)
        
        if self.current_dict[key].shape[0] != self.current_dict['nnod']-off:
            show_message("The number of Cd Ports must match the Number of Nodes")
            self.current_dict[key] = []
            return
        
        self.current_dict[key] = list(self.current_dict[key]) # para que no sea numpy array
        self.set_table(key,table,self.current_dict)
        return
    
    
    def match_list(self, new_number, cur_list, new_item):
        while new_number!=len(cur_list):
            if new_number<len(cur_list):
                cur_list.pop()
            else:
                cur_list.append(new_item)
        return
    
    def modify_nnodes(self, ncd):
        self.match_list(ncd-1,self.current_dict['Cd_ports'],0.8)
        self.match_list(ncd-1,self.current_dict['state_ini'],[1.1769, 0.1, 101330.0])
        
        self.set_table('Cd_ports',self.ui_td.Cd_ports,self.current_dict)
        self.set_table('state_ini',self.ui_td.state_ini,self.current_dict)
        return
    
    def parse_data(self, dict_t, ikey):
        """
        parses the dict data to gui values
        """
        if ikey not in PARSED.keys():
            return dict_t[ikey]
        newvalue = 1.0
        # aquellos que se multiplican
        for ii in PARSED[ikey][0]:
            if (type(ii)==str):
                newvalue = newvalue*dict_t[ii]
            else:
                newvalue = newvalue*ii
        # aquellos que se dividen
        for ii in PARSED[ikey][1]:
            if (type(ii)==str):
                newvalue = newvalue/dict_t[ii]
            else:
                newvalue = newvalue/ii
        if ikey in EXTRAS.keys():
            newvalue = newvalue + EXTRAS[ikey]
        return newvalue
    
    def set_table(self, ikey, table, current_dict):
        table.clearContents()
        table.setRowCount(0)
        ncols = table.columnCount()
        for ituple in current_dict[ikey]:
            current_row = table.rowCount()
            table.insertRow(current_row)
            table.setRowCount(current_row + 1)
            for icol in range(0, ncols):
                val = -1
                # tengo que verificar si es una lista, array de numpy, etc, 
                # porque sino no funciona el indexado
                try:
                    val = ituple[icol]
                except:
                    val = ituple
                it = QtWidgets.QTableWidgetItem(str(val))
                it.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable)
                table.setItem(current_row, icol, it)
        return
    
    def set_list(self, ikey, item_tw, current_dict):
        for index_item in current_dict[ikey]:
            it = item_tw.item(index_item)
            if it is not None:
                it.setCheckState(2)
        return
    
    def set_parameters(self):
        """
        puts the loaded tank parameters in the items
        """
        if not self.current_dict:
            return
        
        # El state_ini tiene el inicial para el nodo 0 del tanque, luego el resto
        # son de tubos. Deben separarse en dos tablas pero el diseño original contemplaba
        # todo en una sola lista. Por ello las divido aca, creando este array "nuevo"
        self.current_dict['state_ini_0'] = []
        self.current_dict['state_ini_0'].append(self.current_dict['state_ini'][0])
        self.current_dict['state_ini'].pop(0)

        for itab in range(3): # por cada tab del widget

            current_dict = self.current_dict

            for ikey in current_dict.keys():
                if ikey not in JSON_TANK_KEYS:
                    continue
                item_le  = self.ui_td.tabWidget.widget(itab).findChild(QtWidgets.QLineEdit,ikey)
                item_cob = self.ui_td.tabWidget.widget(itab).findChild(QtWidgets.QComboBox,ikey)
                item_chb = self.ui_td.tabWidget.widget(itab).findChild(QtWidgets.QCheckBox,ikey)
                item_sb  = self.ui_td.tabWidget.widget(itab).findChild(QtWidgets.QSpinBox,ikey)
                item_dsb = self.ui_td.tabWidget.widget(itab).findChild(QtWidgets.QDoubleSpinBox,ikey)
                item_tw  = self.ui_td.tabWidget.widget(itab).findChild(QtWidgets.QTableWidget,ikey)
                item_lw  = self.ui_td.tabWidget.widget(itab).findChild(QtWidgets.QListWidget,ikey)
                if item_le is not None:
                    parsed_ikey = self.parse_data(current_dict,ikey)
                    # intentar convertir la cadena a string (si es un numero), de otro
                    # modo es directamente un string (por ejemplo, label)
                    try:
                        text = str(parsed_ikey)
                    except:
                        text = parsed_ikey
                    item_le.setText(text)
                elif item_cob is not None:
                    parsed_ikey = self.parse_data(current_dict,ikey)
                    assert(type(parsed_ikey) == int or type(parsed_ikey) == float)
                    item_cob.setCurrentIndex(parsed_ikey)
                elif item_sb is not None:
                    assert(type(current_dict[ikey]) == int)
                    item_sb.setValue(current_dict[ikey])
                elif item_dsb is not None:
                    parsed_ikey = self.parse_data(current_dict,ikey)
                    assert(type(parsed_ikey) == float)
                    item_dsb.setValue(parsed_ikey)
                elif item_tw is not None:
                    assert(type(current_dict[ikey]) == list)
                    self.set_table(ikey, item_tw, current_dict)
                elif item_chb is not None:
                    typ = type(current_dict[ikey])
                    assert(typ==float or typ==int)
                    item_chb.setChecked(current_dict[ikey])
                elif item_lw is not None:
                    assert(type(current_dict[ikey]) == list)
                    self.set_list(ikey, item_lw, current_dict)
                else:
                    continue
        return
    
    def get_table(self, ikey, table, current_dict, factor):
        if not table.isEnabled():
            return
        ncols = table.columnCount()
        nrows = table.rowCount()
        for irow in range(0, nrows):
            dict_row = []
            for icol in range(0, ncols):
                if not table.item(irow,icol):
                    current_dict[ikey] = []
                    show_message('Please provide all the fields for the %s table'%ikey)
                    # dejo que salte el error en el try de accept
                text = table.item(irow,icol).text()
                text_parsed = convert_string(text)
                text_parsed = text_parsed*factor[icol]
                if ncols==1:
                    current_dict[ikey].append(text_parsed)
                else:
                    dict_row.append(text_parsed)
            if dict_row!=[]:
                current_dict[ikey].append(dict_row)
        return
    
    def get_list(self, ikey, list_widget, current_dict):
        current_dict[ikey] = []
        for index_item in range(list_widget.count()):
            if list_widget.item(index_item).checkState()==2:
                current_dict[ikey].append(index_item)
        return

    def accept(self):
        # verificar campos
        to_verify = []
        to_verify.append( check_if_float( self.ui_td.Volume.text(),'Volume') )
        to_verify.append( check_if_float( self.ui_td.mass.text(),'mass') )
        to_verify.append( check_if_float( self.ui_td.h_film.text(),'h_film') )
        to_verify.append( check_if_float( self.ui_td.Area_wall.text(),'Area_wall') )
        to_verify.append( check_if_float( self.ui_td.T_wall.text(),'T_wall') )
        if False in to_verify:
            return
        # lineEdits
        self.current_dict['label']          = str(self.ui_td.label.text())
        self.current_dict['nnod']           = self.ui_td.nnod.value()
        self.current_dict['ndof']           = self.ui_td.ndof.value()
        self.current_dict['Volume']         = float(self.ui_td.Volume.text())*CCM2CM
        self.current_dict['mass']           = float(self.ui_td.mass.text())*MM2M
        self.current_dict['h_film']         = float(self.ui_td.h_film.text())
        self.current_dict['Area_wall']      = float(self.ui_td.Area_wall.text())*SQCM2SQM
        self.current_dict['T_wall']         = float(self.ui_td.T_wall.text())
        self.current_dict['extras']         = int(self.ui_td.extras.isChecked())

        success = False
        # tablas
        try:
            self.current_dict['Cd_ports'] = []
            self.get_table('Cd_ports', self.ui_td.Cd_ports, self.current_dict, [1.0])
            self.current_dict['state_ini'] = []
            self.get_table('state_ini', self.ui_td.state_ini_0, self.current_dict, [1.0,1.0,1.0])
            self.get_table('state_ini', self.ui_td.state_ini, self.current_dict, [1.0,1.0,1.0])
            self.get_list('histo', self.ui_td.histo, self.current_dict)
            success = True
        except:
            success = False
            
        if success:
            self.close()
        return
    
    def cancel(self):
        self.close()
        return
    
    def closeEvent(self, event):
        if not self.sender() or self.sender().objectName()=='cancel_pushButton':
            self.current_dict['state_ini'].insert(0,self.current_dict['state_ini_0'][0])
        self.close()
        return