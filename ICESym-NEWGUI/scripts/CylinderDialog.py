#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 28 20:51:13 2018

@author: etekken
"""

import math, os, json
import numpy as np
import pyqtgraph as pg
from PyQt5 import QtWidgets, QtCore, QtGui
from cylinderDialog_ui import Ui_CylinderDialog
from utils import M2MM, SQM2SQCM, RAD2DEG, MM2M, CCM2CM, SQCM2SQM, DEG2RAD, CM2CCM,\
                  DEFAULT_DVP, INSTALL_PATH, convert_string, set_plot, check_if_float,\
                  show_message, LOADS_PATH

JSON_CYLINDER_KEYS = ["crank_radius", "type_ig", "label", "full_implicit", "ndof", "model_ht",\
                      "factor_ht", "piston_area", "ownState", "mass_C", "nnod", "scavenge_type",\
                      "twall", "state_ini", "nve", "head_chamber_area", "type_temperature",\
                      "rod_length", "species_model", "nvi", "delta_ca", "extras", "histo",\
                      "Vol_clearance", "Bore", "scavenge", "y", "hvap_fuel", "Q_fuel",
                      "combustion", "a_wiebe", "dtheta_comb", "combustion_model", "m_wiebe",\
                      "theta_ig_0", "pulse", "m_inj", "dtheta_inj", "theta_inj_ini", "theta_id",\
                      "ignition_delay_model","T_fuel", "phi", "intake_valves", "exhaust_valves",'state_ini_0']

PARSED = {}
PARSED['crank_radius']          = [['crank_radius',M2MM*2],[]]
PARSED['Bore']                  = [['Bore',M2MM],[]]
PARSED['rod_length']            = [['rod_length',M2MM],[]]
PARSED['piston_area']           = [['piston_area',SQM2SQCM],[]]
PARSED['delta_ca']              = [['delta_ca',RAD2DEG],[]]
PARSED['Vol_clearance']         = [[0.5*math.pi,'Bore','Bore','crank_radius'],['Vol_clearance']]
PARSED['head_chamber_area']     = [['head_chamber_area',SQM2SQCM],[]]
PARSED['model_ht']              = [['model_ht'],[]]
PARSED['hvap_fuel']             = [['hvap_fuel',MM2M],[]]
PARSED['Q_fuel']                = [['Q_fuel',CCM2CM],[]]
PARSED['theta_ig_0']            = [['theta_ig_0',RAD2DEG],[]]
PARSED['dtheta_comb']           = [['dtheta_comb',RAD2DEG],[]]
PARSED['m_inj']                 = [['m_inj',M2MM],[]]
PARSED['dtheta_inj']            = [['dtheta_inj',RAD2DEG],[]]
PARSED['theta_inj_ini']         = [['theta_inj_ini',RAD2DEG],[]]
PARSED['theta_id']              = [['theta_id',RAD2DEG],[]]
PARSED['pulse']                 = [['pulse'],[]]
PARSED['ignition_delay_model']  = [['ignition_delay_model'],[]]

DEPARSED = {}
DEPARSED['crank_radius']          = [['crank_radius',MM2M/2.0],[]]
DEPARSED['Bore']                  = [['Bore',MM2M],[]]
DEPARSED['rod_length']            = [['rod_length',MM2M],[]]
DEPARSED['piston_area']           = [['piston_area',SQCM2SQM],[]]
DEPARSED['delta_ca']              = [['delta_ca',DEG2RAD],[]]
DEPARSED['Vol_clearance']         = [['Vol_clearance'],[0.5*math.pi,'Bore','Bore','crank_radius']] #TODO
DEPARSED['head_chamber_area']     = [['head_chamber_area',SQCM2SQM],[]]
DEPARSED['model_ht']              = [['model_ht'],[]]
DEPARSED['hvap_fuel']             = [['hvap_fuel',M2MM],[]]
DEPARSED['Q_fuel']                = [['Q_fuel',CM2CCM],[]]
DEPARSED['theta_ig_0']            = [['theta_ig_0',DEG2RAD],[]]
DEPARSED['dtheta_comb']           = [['dtheta_comb',DEG2RAD],[]]
DEPARSED['m_inj']                 = [['m_inj',MM2M],[]]
DEPARSED['dtheta_inj']            = [['dtheta_inj',DEG2RAD],[]]
DEPARSED['theta_inj_ini']         = [['theta_inj_ini',DEG2RAD],[]]
DEPARSED['theta_id']              = [['theta_id',DEG2RAD],[]]
DEPARSED['pulse']                 = [['pulse'],[]]
DEPARSED['ignition_delay_model']  = [['ignition_delay_model'],[]]

CB_GET_VALUE = {}

EXTRAS = {}
EXTRAS['Vol_clearance'] = 1.0
EXTRAS['model_ht']      = -1
EXTRAS['pulse']         = -1

DEEXTRAS = {}
DEEXTRAS['Vol_clearance'] = -1.0
DEEXTRAS['model_ht']      = 1
DEEXTRAS['pulse']         = 1

HOMO_TEMP_TABLE     = ['Wall']
NOHOMO_TEMP_TABLE   = ['Piston','Intake','Exhaust','Liners']

DICT_KEYS = ['fuel','injection','combustion']

class CylinderDialog(QtWidgets.QDialog):
    """
    class to manage the cylinder atributes. If current_cylinder is None, we are 
    creating a new one. On the other hand, we are modifying an old one.
    """
    def __init__(self, current_cylinder = None, item_index = 0):
        QtWidgets.QDialog.__init__(self)
        self.ui_cd = Ui_CylinderDialog()
        self.ui_cd.setupUi(self)
        self.setBaseSize(654, 873)
        self.set_restrictions()
        self.current_dict = None # default cylinder dictionary
        self.setWindowTitle( self.windowTitle() + " " + str(item_index) )
        if not current_cylinder:
            self.load_default()
        else:
            self.current_dict = current_cylinder
        self.set_parameters()
        self.change_ignition_delay_model( self.ui_cd.ignition_delay_model.currentIndex() )
        self.change_pulse( self.ui_cd.pulse.currentIndex() )
        self.combustion_model_selection( self.ui_cd.combustion_model.currentIndex() )
        self.change_ignition_type( self.ui_cd.type_ig.currentIndex() )
        return
    
    def set_restrictions(self):
        self.ui_cd.Bore.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.crank_radius.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.rod_length.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.Vol_clearance.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.head_chamber_area.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.piston_area.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.delta_ca.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.hvap_fuel.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.m_inj.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.dtheta_inj.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.T_fuel.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.theta_inj_ini.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.theta_id.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.theta_ig_0.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.dtheta_comb.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.phi.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.a_wiebe.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_cd.m_wiebe.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        return
    
    def load_default(self):
        """
        load the default cylinder template
        """
        filename = os.path.join(INSTALL_PATH,"templates","cylinder_default.json")
        if not os.path.isfile(filename):
            show_message("Cannot find default cylinder configuration file")
            return
        with open(filename) as openedfile:
            try:
                self.current_dict = json.load(openedfile)
            except ValueError as error:
                show_message('JSON object issue: %s'%error)
                return
        if self.current_dict is None:
            show_message("An error ocurred with the json default cylinder archive")
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
        for ikey in JSON_CYLINDER_KEYS:
            if ikey not in self.current_dict.keys() in (self.current_dict['fuel'].keys() or self.current_dict['combustion'].keys() or self.current_dict['injection'].keys()):
                show_message("Wrong number of keys in json default cylinder archive")
                return False
        return True

    def parse_data(self, dict_t, ikey, dict_to_use=PARSED, dict_extras=EXTRAS):
        """
        parses the dict data to gui values
        """
        if ikey not in dict_to_use.keys():
            return dict_t[ikey]

        if ikey=='Vol_clearance' and dict_to_use==DEPARSED:
            val = (0.5*math.pi*(float(self.ui_cd.Bore.text())*MM2M)**2*float(self.ui_cd.crank_radius.text())/2.0*MM2M)
            val = val/(float(self.ui_cd.Vol_clearance.text())-1.0)
            return val
        
        newvalue = 1.0
        # aquellos que se multiplican
        for ii in dict_to_use[ikey][0]:
            if (type(ii)==str):
                newvalue = newvalue*dict_t[ii]
            else:
                newvalue = newvalue*ii
        # aquellos que se dividen
        for ii in dict_to_use[ikey][1]:
            if (type(ii)==str):
                newvalue = newvalue/dict_t[ii]
            else:
                newvalue = newvalue/ii
        if ikey in EXTRAS.keys():
            newvalue = newvalue + dict_extras[ikey]
        return newvalue
    
    def load_table(self):
        """
        loads a defined archive with values for the multiple tables
        """
        key = ''
        if self.sender().objectName()=="load_mfdot_pushButton":
            key = 'mfdot_array'
            current_dict = self.current_dict['injection']
            table = self.ui_cd.mfdot_array
        elif self.sender().objectName()=="load_xbdot_pushButton":
            key = 'xbdot_array'
            current_dict = self.current_dict['combustion']
            table = self.ui_cd.xbdot_array
        elif self.sender().objectName()=="load_state_ini_pushButton":
            key = 'state_ini'
            current_dict = self.current_dict
            table = self.ui_cd.state_ini
            
        assert(key!='')
        
        if key=='mfdot_array' and self.ui_cd.pulse.currentIndex() != 2:
            show_message("You must select a Pulse User Defined", 2)
            return

        dialog = QtWidgets.QFileDialog(self)
        dialog.setNameFilter(" (*.txt)")
        dialog.setWindowTitle('Select a File to Open')
        dialog.setDirectory(LOADS_PATH)
        if dialog.exec_():
            filename = dialog.selectedFiles()[0]
        else:
            show_message("Nothing selected", 2)
            return

        if key in current_dict.keys():
            msg = "You will lose the current values. Do you want to continue?"
            reply = show_message(msg,4,QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No)
            if reply == QtWidgets.QMessageBox.Yes:
                current_dict[key] = np.loadtxt(filename)
            else:
                show_message("Operation cancelled", 1)
                return
        else:
            current_dict[key] = np.loadtxt(filename)
            
        if key=='state_ini':
            if self.current_dict[key].shape[0] != self.ui_cd.nnod.value():
                show_message("The number of Inital States must match the Number of Nodes")
                self.current_dict[key] = [[0.0,0.0,0.0]]
                return
        
        current_dict[key] = list(current_dict[key]) # para que no sea numpy array
        self.set_table(key, table, current_dict)
        return
    
    def plot_table(self):
        """
        plot the current table
        """
        key = ''
        if self.sender().objectName()=="plot_mfdot_pushButton":
            key = 'mfdot_array'
            current_dict = self.current_dict['injection']
            xlabel = "Angle [deg]"
            ylabel = "Value [kg/s]"
            title = "Angle vs. Value"
        elif self.sender().objectName()=="plot_xbdot_pushButton":
            key = 'xbdot_array'
            current_dict = self.current_dict['combustion']
            xlabel = "Angle [deg]"
            ylabel = "Value [1/s]"
            title = "Angle vs. Value"
        assert(key!='')

        if key not in current_dict or current_dict[key] == []:
            show_message("Not valid data")
            return
        pg.setConfigOptions(background=None)
        pg.setConfigOptions(foreground='k')
        try:
            xdata = [item[0] for item in current_dict[key]]
            ydata = [item[1] for item in current_dict[key]]
            plot = pg.PlotWidget()
            set_plot(plot, xlabel, ylabel, "")
            dialog = QtWidgets.QDialog()
            dialog.setLayout(QtWidgets.QHBoxLayout())
            dialog.layout().addWidget(plot)
            dialog.setWindowTitle(title)
            plot.plot(xdata, ydata, pen={'color': 'r', 'width': 1})
            dialog.exec_()
        except:
            show_message("Cannot plot the current values")
        return

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

        if ikey=='twall':
            temp_type = self.current_dict['type_temperature']
            table.setVerticalHeaderLabels(HOMO_TEMP_TABLE if \
                                          temp_type==0 else NOHOMO_TEMP_TABLE)
        return

    def change_temp_type(self, index):
        if index==1 and len(self.current_dict['twall'])<4:
            self.current_dict['twall'] = [0.0,0.0,0.0,0.0]
        elif index==0:
            self.current_dict['twall'] = [0.0]
        self.current_dict['type_temperature'] = index
        self.set_table('twall', self.ui_cd.twall, self.current_dict)
        return

    def set_scavenge(self, index):
        self.ui_cd.scavenge_type.setEnabled(False if index==0 else True)
        if index==1 and 'scavenge_type' not in self.current_dict.keys():
            self.current_dict['scavenge_type'] = 0
        return
    
    def combustion_model_selection(self, mode):
        self.ui_cd.xbdot_array.setEnabled(False)
        self.ui_cd.a_wiebe.setEnabled(False)
        self.ui_cd.m_wiebe.setEnabled(False)
        if mode==0: # user-defined
            self.ui_cd.xbdot_array.setEnabled(True)
        elif mode==1 or mode==2: # wiebe or wiebe-3
            self.ui_cd.a_wiebe.setEnabled(True)
            self.ui_cd.m_wiebe.setEnabled(True)
            
        return
    
    def change_ignition_delay_model(self, mode):
        self.ui_cd.theta_id.setEnabled(False)
        if (mode==2): # user-defined
            self.ui_cd.theta_id.setEnabled(True)
        return
    
    def change_ignition_type(self, mode):
        self.ui_cd.diesel_groupBox.setEnabled(False)
        self.ui_cd.combustion_model.model().item(2).setEnabled(False)
        self.ui_cd.combustion_model.model().item(3).setEnabled(False)
        self.ui_cd.phi.setEnabled(False)
        if mode==0:
            self.ui_cd.phi.setEnabled(True)
        if mode==1:
            self.ui_cd.diesel_groupBox.setEnabled(True)
            self.ui_cd.combustion_model.model().item(2).setEnabled(True)
            self.ui_cd.combustion_model.model().item(3).setEnabled(True)
            
        return
    
    def change_pulse(self, mode):
        self.ui_cd.mfdot_array.setEnabled(False)
        if (mode==2): # user-defined
            self.ui_cd.mfdot_array.setEnabled(True)
        return
    
    def match_list(self, new_number, cur_list, new_item):
        while new_number!=len(cur_list):
            if new_number<len(cur_list):
                cur_list.pop()
            else:
                cur_list.append(new_item)
        return
    
    def modify_nnodes(self, ncd):
        self.match_list(ncd-1,self.current_dict['state_ini'],DEFAULT_DVP)
        self.set_table('state_ini',self.ui_cd.state_ini,self.current_dict)
        return

    def set_parameters(self):
        """
        puts the loaded cyilinder parameters in the items
        """
        if not self.current_dict:
            return

        # El state_ini tiene el inicial para el nodo 0 del tanque, luego el resto
        # son de tubos. Deben separarse en dos tablas pero el diseño original contemplaba
        # todo en una sola lista. Por ello las divido aca, creando este array "nuevo"
        self.current_dict['state_ini_0'] = []
        self.current_dict['state_ini_0'].append(self.current_dict['state_ini'][0])
        self.current_dict['state_ini'].pop(0)

        tabs = (0,4,5)

        for itab in range(6): # por cada tab del widget

            current_dict = self.current_dict if itab in tabs else \
                                            self.current_dict[DICT_KEYS[itab-1]]

            for ikey in current_dict.keys():

                if ikey not in JSON_CYLINDER_KEYS:
                    continue
                item_le  = self.ui_cd.tabWidget.widget(itab).findChild(QtWidgets.QLineEdit,ikey)
                item_cob = self.ui_cd.tabWidget.widget(itab).findChild(QtWidgets.QComboBox,ikey)
                item_chb = self.ui_cd.tabWidget.widget(itab).findChild(QtWidgets.QCheckBox,ikey)
                item_sb  = self.ui_cd.tabWidget.widget(itab).findChild(QtWidgets.QSpinBox,ikey)
                item_dsb = self.ui_cd.tabWidget.widget(itab).findChild(QtWidgets.QDoubleSpinBox,ikey)
                item_tw  = self.ui_cd.tabWidget.widget(itab).findChild(QtWidgets.QTableWidget,ikey)
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
                    typecb = type(current_dict[ikey])
                    assert(typecb == float or typecb==int or typecb==list)
                    if typecb==list:
                        item_chb.setChecked(True if 0 in current_dict[ikey] else False)
                    else:
                        item_chb.setChecked(current_dict[ikey])
                else:
                    continue
        return
    
    def get_table(self, ikey, table, current_dict, factor=1.0):
        if not table.isEnabled():
            return True
        ncols = table.columnCount()
        nrows = table.rowCount()
        current_dict[ikey] = [] if ikey!='state_ini' else current_dict[ikey]
        for irow in range(0, nrows):
            dict_row = []
            for icol in range(0, ncols):
                if not table.item(irow,icol):
                    current_dict[ikey] = []
                    show_message('Please provide all the fields for the %s table'%ikey)
                    return False
                if not check_if_float(table.item(irow,icol).text(), ikey):
                    show_message('Please provide valid values for the %s table'%ikey)
                    current_dict[ikey] = []
                    return False
                text = table.item(irow,icol).text()
                text_parsed = convert_string(text)
                text_parsed = text_parsed*factor
                if ncols==1:
                    current_dict[ikey].append(text_parsed)
                else:
                    dict_row.append(text_parsed)
            if dict_row!=[]:
                current_dict[ikey].append(dict_row)
        return True
    
    def verify_fields(self):
        to_verify = []
        to_verify.append( check_if_float( self.ui_cd.Bore.text(),'Bore') )
        to_verify.append( check_if_float( self.ui_cd.crank_radius.text(),'Stroke') )
        to_verify.append( check_if_float( self.ui_cd.rod_length.text(),'Con-rod Length') )
        to_verify.append( check_if_float( self.ui_cd.Vol_clearance.text(),'Compression Ratio') )
        to_verify.append( check_if_float( self.ui_cd.head_chamber_area.text(),'Head Chamber Area') )
        to_verify.append( check_if_float( self.ui_cd.piston_area.text(),'Piston Area') )
        to_verify.append( check_if_float( self.ui_cd.delta_ca.text(),'Piston Area') )
        to_verify.append( check_if_float( self.ui_cd.hvap_fuel.text(),'Hvap Fuel') )
        if self.current_dict['type_ig']==1:
            to_verify.append( check_if_float( self.ui_cd.m_inj.text(),'Mass Injected p/ Cycle') )
            to_verify.append( check_if_float( self.ui_cd.dtheta_inj.text(),'Duration of Injection') )
            to_verify.append( check_if_float( self.ui_cd.T_fuel.text(),'Fuel Temperature') )
            to_verify.append( check_if_float( self.ui_cd.theta_inj_ini.text(),'Start of Injection') )
            to_verify.append( check_if_float( self.ui_cd.theta_id.text(),'Ignition Delay') )
        to_verify.append( check_if_float( self.ui_cd.theta_ig_0.text(),'Start of Combustion') )
        to_verify.append( check_if_float( self.ui_cd.dtheta_comb.text(),'Duration of Combustion') )
        to_verify.append( check_if_float( self.ui_cd.phi.text(),'Equivalence Ratio') )
        to_verify.append( check_if_float( self.ui_cd.a_wiebe.text(),'Comb. Efficiency Parameter') )
        to_verify.append( check_if_float( self.ui_cd.m_wiebe.text(),'Shape Parameter') )
        if False in to_verify:
            return False
        return True
    
    def set_histogram(self):
        c = 1
        for itype_valves in ('intake_valves','exhaust_valves'):
            for ivalve in self.current_dict[itype_valves]:
                if ivalve['histo']==1:
                    self.current_dict['histo'].append(c)
                c+=1
        return
    
    def accept(self):
        """
        puts the loaded configuration parameters in the dict and save it
        """
        if not self.current_dict:
            return

        if not self.verify_fields():
            return

        self.current_dict['histo'] = []
        tabs = (0,4,5)
        
        self.current_dict['state_ini'] = []
        ret = self.get_table('state_ini_0', self.ui_cd.state_ini_0, self.current_dict)
        self.current_dict['state_ini'].append(self.current_dict['state_ini_0'][0])
        del self.current_dict['state_ini_0']

        for itab in range(6): # por cada tab del widget

            current_dict = self.current_dict if itab in tabs else \
                                            self.current_dict[DICT_KEYS[itab-1]]

            for ikey in current_dict.keys():

                if ikey not in JSON_CYLINDER_KEYS:
                    continue
                item_le  = self.ui_cd.tabWidget.widget(itab).findChild(QtWidgets.QLineEdit,ikey)
                item_cob = self.ui_cd.tabWidget.widget(itab).findChild(QtWidgets.QComboBox,ikey)
                item_chb = self.ui_cd.tabWidget.widget(itab).findChild(QtWidgets.QCheckBox,ikey)
                item_sb  = self.ui_cd.tabWidget.widget(itab).findChild(QtWidgets.QSpinBox,ikey)
                item_dsb = self.ui_cd.tabWidget.widget(itab).findChild(QtWidgets.QDoubleSpinBox,ikey)
                item_tw  = self.ui_cd.tabWidget.widget(itab).findChild(QtWidgets.QTableWidget,ikey)
                if item_le is not None:
                    text = item_le.text()
                    text = convert_string(text)
                    current_dict[ikey] = text
                    current_dict[ikey] = self.parse_data(current_dict,ikey,DEPARSED,DEEXTRAS)
                elif item_cob is not None:
                    if ikey in CB_GET_VALUE:
                        current_dict[ikey] = convert_string(item_cob.currentText())
                    else:
                        current_dict[ikey] = item_cob.currentIndex()
                    current_dict[ikey] = self.parse_data(current_dict,ikey,DEPARSED,DEEXTRAS)
                elif item_sb is not None:
                    current_dict[ikey] = item_sb.value()
                    current_dict[ikey] = self.parse_data(current_dict,ikey,DEPARSED,DEEXTRAS)
                elif item_dsb is not None:
                    current_dict[ikey] = item_dsb.value()
                    current_dict[ikey] = self.parse_data(current_dict,ikey,DEPARSED,DEEXTRAS)
                elif item_tw is not None:
                    ret = self.get_table(ikey, item_tw, current_dict)
                    if not ret:
                        return
                elif item_chb is not None:
                    if type(current_dict[ikey])==list:
                        current_dict[ikey] = [0] if item_chb.isChecked() else []
                    else:
                        current_dict[ikey] = 1 if item_chb.isChecked() else 0
                        current_dict[ikey] = self.parse_data(current_dict,ikey,DEPARSED,DEEXTRAS)
                else:
                    continue

        self.set_histogram()
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