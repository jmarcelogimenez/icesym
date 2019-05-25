#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 17:13:17 2018

@author: etekken
"""

import os, json
from PyQt5 import QtWidgets, QtCore
from junctionDialog_ui import Ui_JunctionDialog
from utils import show_message, INSTALL_PATH

JSON_JUNCTION_KEYS = ['label','ndof','nnod','type_end','node2tube','histo', 'extras']

PARSED = {}
PARSED['modelo_junc'] = [['modelo_junc'],[]]

EXTRAS = {}
EXTRAS['modelo_junc'] = -1

class JunctionDialog(QtWidgets.QDialog):
    """
    class to manage the junction atributes. If current_tank is None, we are 
    creating a new one. On the other hand, we are modifying an old one.
    """
    def __init__(self, current_junction = None, item_index = 0):
        QtWidgets.QDialog.__init__(self)
        self.ui_jd = Ui_JunctionDialog()
        self.ui_jd.setupUi(self)
        self.setFixedSize(358, 300)
        self.current_dict = None # default junction dictionary
        self.setWindowTitle( self.windowTitle() + " " + str(item_index) )
        if not current_junction:
            self.load_default()
        else:
            self.current_dict = current_junction
        self.set_tubes_connections()
        self.set_parameters()
        return
    
    def load_default(self):
        """
        load the default junction template
        """
        filename = os.path.join(INSTALL_PATH,"templates","junction_default.json")
        if not os.path.isfile(filename):
            show_message("Cannot find default junction configuration file")
            return
        with open(filename) as openedfile:
            try:
                self.current_dict = json.load(openedfile)
            except ValueError, error:
                print 'JSON object issue: %s'%error
        if self.current_dict is None:
            show_message("An error ocurred with the json default junction archive")
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
        for ikey in JSON_JUNCTION_KEYS:
            if ikey not in self.current_dict.keys():
                show_message("Wrong number of keys in json default junction archive")
                return False
        return True
    
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
                it.setFlags(QtCore.Qt.ItemIsEnabled)
                table.setItem(current_row, icol, it)
        return

    def set_tubes_connections(self):
        self.ui_jd.histo.setSpacing(3)
        for itube in self.current_dict['node2tube']:
            name = "Connection with tube %s"%itube
            it = QtWidgets.QListWidgetItem(name)
            it.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsUserCheckable)
            it.setCheckState(0)
            self.ui_jd.histo.addItem(it)
        return

    def set_list(self, ikey, item_tw, current_dict):
        for index_item in current_dict[ikey]:
            it = item_tw.item(index_item)
            if it is not None:
                it.setCheckState(2)
        return
    
    def set_parameters(self):
        """
        puts the loaded junction parameters in the items
        """
        if not self.current_dict:
            return

        for itab in range(2): # por cada tab del widget

            current_dict = self.current_dict

            for ikey in JSON_JUNCTION_KEYS and current_dict.keys():
                item_le  = self.ui_jd.tabWidget.widget(itab).findChild(QtWidgets.QLineEdit,ikey)
                item_cob = self.ui_jd.tabWidget.widget(itab).findChild(QtWidgets.QComboBox,ikey)
                item_chb = self.ui_jd.tabWidget.widget(itab).findChild(QtWidgets.QCheckBox,ikey)
                item_sb  = self.ui_jd.tabWidget.widget(itab).findChild(QtWidgets.QSpinBox,ikey)
                item_dsb = self.ui_jd.tabWidget.widget(itab).findChild(QtWidgets.QDoubleSpinBox,ikey)
                item_tw  = self.ui_jd.tabWidget.widget(itab).findChild(QtWidgets.QTableWidget,ikey)
                item_lw  = self.ui_jd.tabWidget.widget(itab).findChild(QtWidgets.QListWidget,ikey)
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

    def get_list(self, ikey, list_widget, current_dict):
        current_dict[ikey] = []
        for index_item in range(list_widget.count()):
            if list_widget.item(index_item).checkState()==2:
                current_dict[ikey].append(index_item)
        return

    def accept(self):
        self.current_dict['label']          = str(self.ui_jd.label.text())
        self.current_dict['nnod']           = self.ui_jd.nnod.value()
        self.current_dict['ndof']           = self.ui_jd.ndof.value()
        self.current_dict['modelo_junc']    = self.ui_jd.modelo_junc.currentIndex()+1
        self.current_dict['extras']         = int(self.ui_jd.extras.isChecked())

        success = False
        # tablas
        try:
            self.get_list('histo', self.ui_jd.histo, self.current_dict)
            success = True
        except:
            success = False
        if success:
            self.close()
        return
    
    def cancel(self):
        self.close()
        return