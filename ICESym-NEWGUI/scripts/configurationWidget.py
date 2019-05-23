#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 20:37:19 2018

@author: etekken
"""

import os
from PyQt5 import QtCore, QtGui, QtWidgets
from configurationWidget_ui import Ui_ConfigureWidget
from utils import convert_string, INSTALL_PATH

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
    
    
JSON_CONFIGURATION_KEYS = ["dtheta_rpm","filein_state","calc_engine_data","Courant","heat_flow","R_gas",\
                           "rpms","filesave_spd","filesave_state","ncycles","folder_name","ga","viscous_flow",\
                           "nsave","get_state","nappend","engine_type","nstroke","ig_order"]

CB_GET_VALUE = ["nstroke"]

PARSED = {}
EXTRAS = {}

class configurationWidget(QtWidgets.QWidget):
    def __init__(self, current_configuration = None, case_name = 'default_case'):
        QtWidgets.QWidget.__init__(self)
        self.ui_cw = Ui_ConfigureWidget()
        self.ui_cw.setupUi(self)
        self.set_restrictions()
        self.current_configuration = current_configuration
        self.set_parameters()
        self.ui_cw.case_name.setText(case_name)
        return

    def set_restrictions(self):
        """
        set the initial restrictions on lineEdits
        """
        self.ui_cw.ncycles.setValidator(QtGui.QIntValidator())
        self.ui_cw.dtheta_rpm.setValidator(QtGui.QDoubleValidator(0, 360, 3))
        self.ui_cw.Courant.setValidator(QtGui.QDoubleValidator(0, 1000, 3))
        self.ui_cw.ga.setValidator(QtGui.QDoubleValidator(0, 1000, 3))
        self.ui_cw.viscous_flow.setValidator(QtGui.QDoubleValidator(0, 1000, 3))
        self.ui_cw.heat_flow.setValidator(QtGui.QDoubleValidator(0, 1000, 3))
        self.ui_cw.nappend.setValidator(QtGui.QDoubleValidator(0, 1000, 3))
        self.ui_cw.nsave.setValidator(QtGui.QIntValidator())
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
                print dict_t[ii]
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
                if ikey=='rpms':
                    it.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable) 
                elif ikey=='ig_order': 
                    it.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsDragEnabled|QtCore.Qt.ItemIsDropEnabled)
                table.setItem(current_row, icol, it)

        if ikey=='rpms':
            self.ui_cw.nof_rpms.setValue(table.rowCount())
        elif ikey=='ig_order':
            self.ui_cw.ig_order.setDragEnabled(True)
            self.ui_cw.ig_order.setAcceptDrops(True)
            self.ui_cw.ig_order.viewport().setAcceptDrops(True)
            self.ui_cw.ig_order.setDragDropOverwriteMode(False)
            self.ui_cw.ig_order.setDropIndicatorShown(True)
            self.ui_cw.ig_order.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection) 
            self.ui_cw.ig_order.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.ui_cw.ig_order.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)   
            self.ui_cw.ig_order.dropEvent = self.dropEvent
        return

    def get_table(self, ikey, table, current_dict):
        ncols = table.columnCount()
        nrows = table.rowCount()
        current_dict[ikey] = []
        for irow in range(0,nrows):
            for icol in range(0,ncols):
                text = table.item(irow,icol).text()
                text = convert_string(text)
                current_dict[ikey].append(text)
        return
    
    def edit_ig_order(self, cyl_index, delete = True):
        if delete:
            self.current_configuration['ig_order'].remove(cyl_index)
            # restar uno a los demas
            for index,icyl in enumerate(self.current_configuration['ig_order']):
                    icyl = icyl if icyl<cyl_index else icyl-1
                    self.current_configuration['ig_order'][index] = icyl
        else:
            self.current_configuration['ig_order'].append(cyl_index)
        self.set_table('ig_order', self.ui_cw.ig_order, self.current_configuration)
        return

    def edit_rpms(self, nof_rpms):
        self.ui_cw.rpms.setRowCount(nof_rpms)
        return
    
    def enable_load_button(self, cur_index):
        if cur_index==1:
            self.ui_cw.toolButton_get_state.setEnabled(True)
            if 'filein_state' in self.current_configuration.keys():
                filein = self.current_configuration['filein_state']
                filein = os.path.basename(os.path.normpath(filein))
                self.ui_cw.filein_state.setText(filein if filein!='' and filein!='None' else 'None')
        else:
            self.ui_cw.toolButton_get_state.setEnabled(False)
            self.ui_cw.filein_state.setText('None')
            self.current_configuration['filein_state'] = 'None'
        return
    
    def load_state_file(self):
        dialog = QtGui.QFileDialog(self)
        dialog.setWindowTitle('Select a File')
        dialog.setDirectory(INSTALL_PATH)
        if dialog.exec_():
            # TODO: ver como validar esto...
            filename = dialog.selectedFiles()[0]
            self.current_configuration['filein_state'] = filename
            filename = os.path.basename(os.path.normpath(filename))
            self.ui_cw.filein_state.setText(filename)
        return
        

    def set_parameters(self):
        """
        puts the loaded configuration parameters in the items
        """
        if not self.current_configuration:
            return
        
        groupboxes = [self.ui_cw.groupBox,self.ui_cw.groupBox_2,self.ui_cw.groupBox_3]
        
        for igroupbox in groupboxes: # por cada groupbox
            
            for ikey in self.current_configuration.keys():
                if ikey not in JSON_CONFIGURATION_KEYS:
                    continue
                item_le  = igroupbox.findChild(QtWidgets.QLineEdit,ikey)
                item_cob = igroupbox.findChild(QtWidgets.QComboBox,ikey)
                item_chb = igroupbox.findChild(QtWidgets.QCheckBox,ikey)
                item_sb  = igroupbox.findChild(QtWidgets.QSpinBox,ikey)
                item_dsb = igroupbox.findChild(QtWidgets.QDoubleSpinBox,ikey)
                item_tw  = igroupbox.findChild(QtWidgets.QTableWidget,ikey)
                if item_le is not None:
                    parsed_ikey = self.parse_data(self.current_configuration,ikey)
                    # intentar convertir la cadena a string (si es un numero), de otro
                    # modo es directamente un string (por ejemplo, label)
                    try:
                        text = str(parsed_ikey)
                    except:
                        text = parsed_ikey
                    item_le.setText(text)
                elif item_cob is not None:
                    parsed_ikey = self.parse_data(self.current_configuration,ikey)
                    assert(type(parsed_ikey) == int or type(parsed_ikey) == float)
                    item_cob.setCurrentIndex(parsed_ikey)
                    # si no hay indice valido, seguramente no es el indice sino
                    # el texto del cb el que hay que buscar y setear
                    if item_cob.currentIndex()==-1:
                        try:
                            item_cob.setCurrentIndex(item_cob.findText(str(parsed_ikey)))
                        except:
                            continue
                elif item_sb is not None:
                    assert(type(self.current_configuration[ikey]) == int)
                    item_sb.setValue(self.current_configuration[ikey])
                elif item_dsb is not None:
                    parsed_ikey = self.parse_data(self.current_configuration,ikey)
                    assert(type(parsed_ikey) == float)
                    item_dsb.setValue(parsed_ikey)
                elif item_tw is not None:
                    assert(type(self.current_configuration[ikey]) == list)
                    self.set_table(ikey, item_tw, self.current_configuration)
                elif item_chb is not None:
                    assert(type(self.current_configuration[ikey]) == int)
                    item_chb.setChecked(self.current_configuration[ikey])
                else:
                    continue
        return

    def save_data(self):
        """
        puts the loaded configuration parameters in the dict and save it
        """
        if not self.current_configuration:
            return
        
        groupboxes = [self.ui_cw.groupBox,self.ui_cw.groupBox_2,self.ui_cw.groupBox_3]
        
        for igroupbox in groupboxes: # por cada groupbox
            
            for ikey in self.current_configuration.keys():
                
                if ikey not in JSON_CONFIGURATION_KEYS:
                    continue
                item_le  = igroupbox.findChild(QtWidgets.QLineEdit,ikey)
                item_cob = igroupbox.findChild(QtWidgets.QComboBox,ikey)
                item_chb = igroupbox.findChild(QtWidgets.QCheckBox,ikey)
                item_sb  = igroupbox.findChild(QtWidgets.QSpinBox,ikey)
                item_dsb = igroupbox.findChild(QtWidgets.QDoubleSpinBox,ikey)
                item_tw  = igroupbox.findChild(QtWidgets.QTableWidget,ikey)
                if item_le is not None:
                    text = item_le.text()
                    # intentar convertir la cadena a float (si es un numero), de otro
                    # modo es directamente un string (por ejemplo, label)
                    text = convert_string(text)
                    self.current_configuration[ikey] = text
                elif item_cob is not None:
                    if ikey in CB_GET_VALUE:
                        self.current_configuration[ikey] = convert_string(item_cob.currentText())
                    else:
                        self.current_configuration[ikey] = item_cob.currentIndex()
                elif item_sb is not None:
                    self.current_configuration[ikey] = item_sb.value()
                elif item_dsb is not None:
                    self.current_configuration[ikey] = item_dsb.value()
                elif item_tw is not None:
                    self.get_table(ikey, item_tw, self.current_configuration)
                elif item_chb is not None:
                    self.current_configuration[ikey] = 1 if item_chb.isChecked() else 0
                else:
                    continue
        
        return self.current_configuration

# ------- Basura necesaria para el drag and drop de la tabla ig_order

    def position(self, pos, rect, index):
        r = QtGui.QAbstractItemView.OnViewport
        margin = 2
        if pos.y() - rect.top() < margin:
            r = QtGui.QAbstractItemView.AboveItem
        elif rect.bottom() - pos.y() < margin:
            r = QtGui.QAbstractItemView.BelowItem
        elif rect.contains(pos, True):
            r = QtGui.QAbstractItemView.OnItem

        if r == QtGui.QAbstractItemView.OnItem and not (self.ui_cw.ig_order.model().flags(index) & QtCore.Qt.ItemIsDropEnabled):
            r = QtGui.QAbstractItemView.AboveItem if pos.y() < rect.center().y() else QtGui.QAbstractItemView.BelowItem

        return r

    def getSelectedRowsFast(self):
        selRows = []
        for item in self.ui_cw.ig_order.selectedItems():
            if item.row() not in selRows:
                selRows.append(item.row())
        return selRows

    def droppingOnItself(self, event, index):
        dropAction = event.dropAction()

        if self.ui_cw.ig_order.dragDropMode() == QtGui.QAbstractItemView.InternalMove:
            dropAction = QtCore.Qt.MoveAction

        if event.source() == self.ui_cw.ig_order and event.possibleActions() & QtCore.Qt.MoveAction and dropAction == QtCore.Qt.MoveAction:
            selectedIndexes = self.ui_cw.ig_order.selectedIndexes()
            child = index
            while child.isValid() and child != self.ui_cw.ig_order.rootIndex():
                if child in selectedIndexes:
                    return True
                child = child.parent()

        return False

    def dropOn(self, event):
        if event.isAccepted():
            return False, None, None, None

        index = QtCore.QModelIndex()
        row = -1
        col = -1

        if self.ui_cw.ig_order.viewport().rect().contains(event.pos()):
            index = self.ui_cw.ig_order.indexAt(event.pos())
            if not index.isValid() or not self.ui_cw.ig_order.visualRect(index).contains(event.pos()):
                index = self.ui_cw.ig_order.rootIndex()

        if self.ui_cw.ig_order.model().supportedDropActions() & event.dropAction():
            if index != self.ui_cw.ig_order.rootIndex():
                dropIndicatorPosition = self.position(event.pos(), self.ui_cw.ig_order.visualRect(index), index)

                if dropIndicatorPosition == QtGui.QAbstractItemView.AboveItem:
                    row = index.row()
                    col = index.column()
                elif dropIndicatorPosition == QtGui.QAbstractItemView.BelowItem:
                    row = index.row() + 1
                    col = index.column()
                else:
                    row = index.row()
                    col = index.column()

            if not self.droppingOnItself(event, index):
                return True, row, col, index

        return False, None, None, None

    def dropEvent(self, event):
        if event.source() == self.ui_cw.ig_order and (event.dropAction() == QtCore.Qt.MoveAction or  self.ui_cw.ig_order.dragDropMode() == QtGui.QAbstractItemView.InternalMove):
            success, row, col, topIndex = self.dropOn(event)
            if success:
                selRows =  self.getSelectedRowsFast()
                top = selRows[0]
                dropRow = row
                if dropRow == -1:
                    dropRow =  self.ui_cw.ig_order.rowCount()
                offset = dropRow - top
                for i, row in enumerate(selRows):
                    r = row + offset
                    if r >  self.ui_cw.ig_order.rowCount() or r < 0:
                        r = 0
                    self.ui_cw.ig_order.insertRow(r)
                selRows =  self.getSelectedRowsFast()
                top = selRows[0]
                offset = dropRow - top                
                for i, row in enumerate(selRows):
                    r = row + offset
                    if r >  self.ui_cw.ig_order.rowCount() or r < 0:
                        r = 0
                    for j in range(self.ui_cw.ig_order.columnCount()):
                        source = QtWidgets.QTableWidgetItem(self.ui_cw.ig_order.item(row, j))
                        self.ui_cw.ig_order.setItem(r, j, source)
                event.accept()
        return