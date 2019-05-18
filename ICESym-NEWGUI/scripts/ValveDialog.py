"""
Created on Wed May  9 11:03:06 2018
@author: etekken
"""

import json
import os
import numpy as np
import pyqtgraph as pg
from PyQt5 import QtGui, QtWidgets, QtCore
from valveDialog_ui import Ui_ValveDialog
from utils import set_plot, check_if_float, show_message, convert_string

JSON_VALVE_KEYS = ['Lvmax', 'angle_VC', 'label', 'Nval', 'Dv', \
                   'type_dat', 'typeVal', 'angle_V0', 'valve_model', 'Cd', 'Lv', \
                   'ncyl', 'tube', 'histo']

class ValveDialog(QtWidgets.QDialog):
    """
    class to manage the valve atributes. If current_valve is None, we are creating
    a new one. On the other hand, we are modifying an old one.
    """
    def __init__(self, current_dir, current_valve_scene_item = None, item_index = 0, parent = None):
        QtWidgets.QDialog.__init__(self)
        self.ui_vd = Ui_ValveDialog()
        self.ui_vd.setupUi(self)
        self.setFixedSize(403, 460)
        self.set_restrictions()
        self.current_dict = None # default valve dictionary
        self.current_dir = current_dir
        self.setWindowTitle( self.windowTitle() + " " + str(item_index) )
        if not current_valve_scene_item:
            self.load_default()
        else:
            self.current_dict = current_valve_scene_item.object
            self.set_parameters()
        self.parent = parent
        self.current_scene_item = current_valve_scene_item
#        self.save_valve_template()
        return

    def set_restrictions(self):
        """
        set the initial restrictions on lineEdits
        """
        self.ui_vd.opening_angle_lineEdit.setValidator(QtGui.QDoubleValidator(0, 360, 3))
        self.ui_vd.closing_angle_lineEdit.setValidator(QtGui.QDoubleValidator(0, 360, 3))
        self.ui_vd.diameter_lineEdit.setValidator(QtGui.QDoubleValidator(0, 100, 3))
        self.ui_vd.max_valve_lift_lineEdit.setValidator(QtGui.QDoubleValidator(0, 100, 3))
        return

    def save_valve_template(self):
        """
        save a defined user template
        """
        filename = ("./templates/valve_default.json")
        with open(filename, 'w') as openedfile:
            json.dump(self.current_dict, openedfile)
        return

    def load_default(self):
        """
        load the default valve template
        """
        filename = self.current_dir + "templates/valve_default.json"
        if not os.path.isfile(filename):
            show_message("Cannot find default valve configuration file")
            return
        with open(filename) as openedfile:
            try:
                self.current_dict = json.load(openedfile)
            except ValueError, error:
                print 'JSON object issue: %s'%error
        if self.current_dict is None:
            show_message("An error ocurred with the json default valve archive")
            return
        if not self.check_json_keys():
            return
        self.set_parameters()
        return

    def check_json_keys(self):
        """
        check if the loaded json include all the obligatory keys
        """
        for ikey in JSON_VALVE_KEYS:
            if ikey not in self.current_dict.keys():
                show_message("Wrong number of keys in json default valve archive")
                return False
        return True

    def set_parameters(self):
        """
        puts the loaded valve parameters in the items
        """
        if not self.current_dict:
            return
        self.ui_vd.name_lineEdit.setText(self.current_dict['label'])
        self.ui_vd.type_comboBox.setCurrentIndex(0 if self.current_dict['typeVal']=='int' else 1)
        self.ui_vd.number_spinBox.setValue(int(self.current_dict['Nval']))
        self.ui_vd.opening_angle_lineEdit.setText(str(\
                              np.rad2deg(float(self.current_dict['angle_V0']))))
        self.ui_vd.closing_angle_lineEdit.setText(str(\
                              np.rad2deg(float(self.current_dict['angle_VC']))))
        self.ui_vd.diameter_lineEdit.setText(str(\
                                         float(1e3*self.current_dict['Dv'])))
        self.ui_vd.model_comboBox.setCurrentIndex(self.current_dict['valve_model'])
        self.ui_vd.lift_comboBox.setCurrentIndex(int(self.current_dict['type_dat'])+1)
        self.ui_vd.max_valve_lift_lineEdit.setText(str(\
                                         float(1e3*self.current_dict['Lvmax'])))
        self.ui_vd.histo.setChecked(True if self.current_dict['histo'] else False)
        if 'Lv' in self.current_dict.keys():
            self.set_table('Lv')
        if 'Cd' in self.current_dict.keys():
            self.set_table('Cd')
        return

    def plot_table(self):
        """
        plot the current angle-lift or lift-dc table
        """
        key = ''
        if self.sender().objectName()=="plot_lift_cd_pushButton":
            key = 'Cd'
        elif self.sender().objectName()=="plot_angle_lift_pushButton":
            key = 'Lv'            
        assert(key!='')

        if key not in self.current_dict.keys() or self.current_dict[key] == []:
            show_message("Not valid data")
            return
        pg.setConfigOptions(background=None)
        pg.setConfigOptions(foreground='k')
        xdata = [item[0] for item in self.current_dict[key]]
        ydata = [item[1] for item in self.current_dict[key]]
        plot = pg.PlotWidget()
        set_plot(plot, "Angle [deg]" if key=='Lv' else "Lift[mm]", "Lift [mm]" \
                 if key=='Lv' else "Discharge Coefficient", "")
        dialog = QtWidgets.QDialog()
        dialog.setLayout(QtWidgets.QHBoxLayout())
        dialog.layout().addWidget(plot)
        dialog.setWindowTitle("Angle vs. Lift Plot" if key=='Lv' else \
                              "Lift vs. Discharge Coefficient")
        plot.plot(xdata, ydata, pen={'color': 'r', 'width': 1})
        dialog.exec_()
        return

    def load_table(self):
        """
        loads a defined archive with values for the
        angle-lift or the lift-dc tables
        """
        key = ''
        if self.sender().objectName()=="load_lift_cd_pushButton":
            key = 'Cd'
        elif self.sender().objectName()=="load_angle_lift_pushButton":
            key = 'Lv'
            
        assert(key!='')
        
        if key=='Lv' and self.ui_vd.lift_comboBox.currentIndex() != 2:
            show_message("You must select User Defined", 2)
            return

        dialog = QtWidgets.QFileDialog(self)
        dialog.setNameFilter(" (*.txt)")
        dialog.setWindowTitle('Select a File to Open')
        dialog.setDirectory("./loads")
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

        if len(self.current_dict[key].shape) != 2 or self.current_dict[key].shape[1] != 2:
            show_message("Wrong type of data")
            self.current_dict.pop(key, None)
            return
        self.set_table(key, True)
        return

    def set_table(self, key, loaded=False):
        """
        set the values of the angle-lift or lift-dc table
        """
        table = None
        if key=='Lv':
            table = self.ui_vd.angle_lift_tableWidget
            factor1 = 1.0
            factor2 = 1e3 if not loaded else 1.0
        elif key=='Cd':
            table = self.ui_vd.lift_cd_tableWidget
            factor1 = 1e3 if not loaded else 1.0
            factor2 = 1.0

        assert(table!=None)

        table.clearContents()
        table.setRowCount(0)
        for ituple in self.current_dict[key]:
            val1 = ituple[0]
            val2 = ituple[1]
            current_row = table.rowCount()
            table.setRowCount(current_row+1)
            it0 = QtWidgets.QTableWidgetItem(str(val1*factor1))
            it0.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable) 
            table.setItem(current_row, 0, it0)
            it1 = QtWidgets.QTableWidgetItem(str(val2*factor2))
            it1.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable) 
            table.setItem(current_row, 1, it1)
        return

    def lift_changed(self, cur_str):
        """
        check if the lift is user-defined or not
        """
        if cur_str != "User Defined":
            self.ui_vd.angle_lift_tableWidget.setEnabled(False)
            self.ui_vd.max_valve_lift_lineEdit.setEnabled(True)
        else:
            self.ui_vd.angle_lift_tableWidget.setEnabled(True)
            self.ui_vd.max_valve_lift_lineEdit.setEnabled(False)
        return

    def save_current_configuration(self):
        """
        save in a json file the current valve configuration
        """
        file_dialog = QtGui.QFileDialog()
        file_dialog.setWindowTitle('Save valve template')
        file_dialog.setDirectory(self.current_dir)
        file_dialog.setNameFilter('json files (*.json)')
        file_dialog.setDefaultSuffix('json')
        filename = ''
        if file_dialog.exec_() == QtGui.QFileDialog.Accepted:
            filename = file_dialog.selectedFiles()[0]
        else:
            show_message("Operation cancelled", 1)
            return

        if filename == '' or not '.json' in filename:
            show_message("Wrong file name. Has a .json extension?")
            return
        try:
            with open(filename, 'w') as openedfile:
                json.dump(self.current_dict, openedfile)
            show_message("The file has been successfully saved", 1)
        except:
            show_message("An error as ocurred trying to save the .json file")
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
                text_parsed = text_parsed*factor[icol]
                if ncols==1:
                    current_dict[ikey].append(text_parsed)
                else:
                    dict_row.append(text_parsed)
            if dict_row!=[]:
                current_dict[ikey].append(dict_row)
        return

    def accept(self):
        # verificar campos
        to_verify = []
        to_verify.append( check_if_float( self.ui_vd.opening_angle_lineEdit.text(),'opening angle') )
        to_verify.append( check_if_float( self.ui_vd.closing_angle_lineEdit.text(),'closing angle') )
        to_verify.append( check_if_float( self.ui_vd.diameter_lineEdit.text(),'diameter') )
        to_verify.append( check_if_float( self.ui_vd.max_valve_lift_lineEdit.text(),'maximum valve lift') )
        if False in to_verify:
            return

        # lineEdits
        self.current_dict['label']          = str(self.ui_vd.name_lineEdit.text())
        self.current_dict['Nval']           = self.ui_vd.number_spinBox.value()
        self.current_dict['valve_model']    = self.ui_vd.model_comboBox.currentIndex()
        self.current_dict['typeVal']        = 'int' if self.ui_vd.type_comboBox.currentIndex()==0 else 'exh'
        self.current_dict['type_dat']       = self.ui_vd.lift_comboBox.currentIndex()-1
        self.current_dict['angle_V0']       = np.deg2rad(float(self.ui_vd.opening_angle_lineEdit.text()))
        self.current_dict['angle_VC']       = np.deg2rad(float(self.ui_vd.closing_angle_lineEdit.text()))
        self.current_dict['Lvmax']          = float(self.ui_vd.max_valve_lift_lineEdit.text())*1e-3
        self.current_dict['Dv']             = float(self.ui_vd.diameter_lineEdit.text())*1e-3
        self.current_dict['histo']          = int(self.ui_vd.histo.isChecked())

        success = False
        # tablas
        try:
            self.get_table('Cd', self.ui_vd.lift_cd_tableWidget, self.current_dict, [1e-3,1.0])
            if self.ui_vd.lift_comboBox.currentIndex()==2: # user-defined
                self.get_table('Lv', self.ui_vd.angle_lift_tableWidget, self.current_dict, [1.0,1e-3])
            else:
                self.current_dict['Lv'] = []
            success = True
        except:
            success = False
        if success:
            self.parent.check_valve_type(self.current_scene_item)
            self.close()
        return
    
    def cancel(self):
        self.close()
        return