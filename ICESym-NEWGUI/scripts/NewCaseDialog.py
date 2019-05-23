#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 03:50:15 2019

@author: etekken
"""

import os, copy, math
from PyQt5 import QtWidgets, QtCore, QtGui
from newCaseDialog_ui import Ui_NewCaseDialog
from utils import show_message, load_templates, save_data_aux, check_if_float,\
                  load_cylinder_template, MM2M, DEFAULT_DVP, INSTALL_PATH, CASES_PATH
from configurationWidget import configurationWidget
from SceneItem import SceneItem
from TubeDialog import configure_default_tube
from ValveDialog import configure_default_valve

class NewCaseDialog(QtWidgets.QDialog):
    def __init__(self):
        QtWidgets.QDialog.__init__(self)
        self.ui_ncd = Ui_NewCaseDialog()
        self.ui_ncd.setupUi(self)
        self.setBaseSize(400, 490)
        self.case_name  = None
        self.case_dir   = None
        self.case_type  = None # 1 abierto, 2 nuevo blanco, 3 nuevo wizard
        self.set_wizard_restrictions()
        return
    
    def set_wizard_restrictions(self):
        self.ui_ncd.Bore.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_ncd.crank_radius.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_ncd.rod_length.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_ncd.Vol_clearance.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        self.ui_ncd.vol_constant.setValidator(QtGui.QDoubleValidator(0, 100000, 3))
        return

    def check_new_case(self, state):
        self.ui_ncd.blank_case_radioButton.setChecked(True)
        self.ui_ncd.blank_case_radioButton.setEnabled(state)
        self.ui_ncd.wizard_radioButton.setEnabled(state)
        self.ui_ncd.case_name.clear()
        self.case_name  = None
        self.case_dir   = None
        self.case_type  = None
        self.ui_ncd.case_name.setReadOnly(not state)
        self.ui_ncd.open_toolButton.setEnabled(not state)
        self.ui_ncd.wizard_groupBox.setEnabled(False)
        return

    def enable_wizard(self, state):
        self.ui_ncd.wizard_groupBox.setEnabled(state)
        self.case_type  = 3
        return
    
    def create_case_from_wizard(self):
        default_dict = load_templates()
        objects = {}
        objects['Valves']       = []
        objects['Tubes']        = []
        objects['Atmospheres']  = []
        objects['Junctions']    = []
        objects['Cylinders']    = []
        objects['Tanks']        = []
        X_SIZE = 66
        Y_SIZE = 66

        ncyls = self.ui_ncd.ncyls.value()
        for icyl in range(0,ncyls):
            default_dict['Configurations']['ig_order'].append(icyl)
        # Cambiar los nombres de los folders segun el elegido por el usario
        default_dict['Configurations']['folder_name']    = '%s_folder'%self.case_name
        default_dict['Configurations']['filesave_state'] = '%s_state'%self.case_name
        default_dict['Configurations']['filesave_spd']   = '%s_species'%self.case_name
        default_dict['Configurations']['nstroke']        = int(self.ui_ncd.nstroke.currentText())
        cw = configurationWidget(default_dict['Configurations'], self.case_name)

        INITIAL_Y = Y_SIZE*(ncyls-1)

        # creacion de los objetos
        iobject = copy.deepcopy(default_dict['Atmospheres'])
        position = QtCore.QPoint(0,INITIAL_Y)
        item_atmosphere = SceneItem('Atmospheres', position, iobject)
        objects['Atmospheres'].append(item_atmosphere)

        iobject = copy.deepcopy(default_dict['Tubes'])
        iobject['tright']   = 'tank'
        iobject['tleft']    = 'atmosphere'
        iobject['nright']   = 0
        iobject['nleft']    = 0
        iobject['label']    = 'Tube_0'
        configure_default_tube(iobject)
        position = QtCore.QPoint(X_SIZE*2,INITIAL_Y)
        item_tube = SceneItem('Tubes', position, iobject)
        objects['Tubes'].append(item_tube)

        iobject = copy.deepcopy(default_dict['Tanks'])
        iobject['int2tube'] = [0]
        iobject['nnod']     = ncyls+2
        iobject['label']    = 'Tank_0'
        # un state_ini ya viene en el template, nnod: 1 + 1 tube + ncyls
        for isi in range(ncyls+1):
            iobject['state_ini'].append(DEFAULT_DVP)
            iobject['Cd_ports'].append(0.8)
        position = QtCore.QPoint(X_SIZE*4,INITIAL_Y)
        item_tank = SceneItem('Tanks', position, iobject)
        objects['Tanks'].append(item_tank)
        
        iobject = copy.deepcopy(default_dict['Junctions'])
        iobject['nnod']     = ncyls+1
        iobject['label']    = 'Junction_0'
        # nnod: 1 tube + ncyls
        for isi in range(ncyls+1):
            iobject['state_ini'].append(DEFAULT_DVP)
        position = QtCore.QPoint(X_SIZE*16,INITIAL_Y)
        item_junc = SceneItem('Junctions', position, iobject)
        objects['Junctions'].append(item_junc)

        # TODO: controlar esto!!
        for icyl in range(0,ncyls):
            iobject = copy.deepcopy(default_dict['Valves'])
            iobject['ncyl']     = icyl
            iobject['typeVal']  = 'int'
            iobject['label']    = 'Intake_Valve_%s'%icyl
            configure_default_valve(iobject,default_dict['Configurations']['nstroke'])
            position = QtCore.QPoint(X_SIZE*8,Y_SIZE*icyl*2)
            item_valve = SceneItem('Valves', position, iobject)
            objects['Valves'].append(item_valve)

        # un tubo de entrada, valvula, cilindro, valvula y tubo de salida por cada cilindro
        for icyl in range(0,ncyls):
            iobject = copy.deepcopy(default_dict['Tubes'])
            iobject['tright']   = 'cylinder'
            iobject['tleft']    = 'tank'
            iobject['nright']   = icyl
            iobject['nleft']    = 0
            iobject['label']    = 'Tube_%s'%len(objects['Tubes'])
            configure_default_tube(iobject)
            position = QtCore.QPoint(X_SIZE*6,Y_SIZE*icyl*2)
            item_tube = SceneItem('Tubes', position, iobject)
            objects['Tubes'].append(item_tube)
            index_tube = objects['Tubes'].index(item_tube)
            objects['Tanks'][0].object['exh2tube'].append(index_tube)
            objects['Valves'][icyl].object['tube'] = index_tube
            
            if self.ui_ncd.default_cylinder.isChecked():
                lct = load_cylinder_template(default_dict['Configurations']['nstroke'],\
                                                    int(self.ui_ncd.type_ig.currentIndex()))
                if lct != {}:
                    iobject = lct
                else:
                    return False
            else:
                iobject = copy.deepcopy(default_dict['Cylinders'])
                # atributos de los objetos (todos iguales)
                iobject['type_ig']          = int(self.ui_ncd.type_ig.currentIndex())
                iobject['Bore']             = float(self.ui_ncd.Bore.text())*MM2M
                iobject['crank_radius']     = float(self.ui_ncd.crank_radius.text())*MM2M/2.0
                iobject['rod_length']       = float(self.ui_ncd.rod_length.text())*MM2M
                val = (0.5*math.pi*(float(self.ui_ncd.Bore.text())*MM2M)**2*float(self.ui_ncd.crank_radius.text())/2.0*MM2M)
                val = val/(float(self.ui_ncd.Vol_clearance.text())-1.0)
                iobject['Vol_clearance']    = val
                iobject['head_chamber_area'] = (math.pi*iobject['Bore']**2/4.0)*1.2
                iobject['piston_area']      = (math.pi*iobject['Bore']**2/4.0)*1.1

            iobject['label']            = 'Cylinder_%s'%icyl
            # state ini son 3 siempre (valvula int, cilindro, valvula ext)
            iobject['state_ini']        = [[1.1769, 101330.0, 300.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0]]
            iobject['nnod']             = 3
            
            vol_constant = float(self.ui_ncd.vol_constant.text())
            objects['Tanks'][0].object['Volume'] = (vol_constant*math.pi*iobject['Bore']**2/4.0)*iobject['crank_radius']

            position = QtCore.QPoint(X_SIZE*10,Y_SIZE*icyl*2)
            item_cyl = SceneItem('Cylinders', position, iobject)
            objects['Cylinders'].append(item_cyl)

            iobject = copy.deepcopy(default_dict['Valves'])
            iobject['ncyl']     = icyl
            iobject['typeVal']  = 'exh'
            iobject['label']    = 'Exhaust_Valve_%s'%icyl
            configure_default_valve(iobject,default_dict['Configurations']['nstroke'])
            position = QtCore.QPoint(X_SIZE*12,Y_SIZE*icyl*2)
            item_valve = SceneItem('Valves', position, iobject)
            objects['Valves'].append(item_valve)

            iobject = copy.deepcopy(default_dict['Tubes'])
            iobject['tright']   = 'junction'
            iobject['tleft']    = 'cylinder'
            iobject['nright']   = 0
            iobject['nleft']    = icyl
            iobject['label']    = 'Tube_%s'%len(objects['Tubes'])
            configure_default_tube(iobject)
            position = QtCore.QPoint(X_SIZE*14,Y_SIZE*icyl*2)
            item_tube = SceneItem('Tubes', position, iobject)
            objects['Junctions'][0].object['type_end'].append(1)
            objects['Tubes'].append(item_tube)
            index_tube = objects['Tubes'].index(item_tube)
            objects['Junctions'][0].object['node2tube'].append(index_tube)
            
            objects['Valves'][-1].object['tube'] = index_tube

        iobject = copy.deepcopy(default_dict['Tubes'])
        iobject['tright']   = 'atmosphere'
        iobject['tleft']    = 'junction'
        iobject['nright']   = 1
        iobject['nleft']    = 0
        iobject['label']    = 'Tube_%s'%len(objects['Tubes'])
        configure_default_tube(iobject)
        position = QtCore.QPoint(X_SIZE*18,INITIAL_Y)
        item_tube = SceneItem('Tubes', position, iobject)
        objects['Tubes'].append(item_tube)
        index_tube = objects['Tubes'].index(item_tube)
        objects['Junctions'][0].object['type_end'].append(-1)
        objects['Junctions'][0].object['node2tube'].append(index_tube)

        iobject = copy.deepcopy(default_dict['Atmospheres'])
        position = QtCore.QPoint(X_SIZE*20,INITIAL_Y)
        item_atmosphere = SceneItem('Atmospheres', position, iobject)
        objects['Atmospheres'].append(item_atmosphere)
        
        for ikey in objects.keys():
            for iobject in objects[ikey]:
                iobject.position = iobject.position+QtCore.QPoint(1.5,1.5)
                iobject.pixmap.setPos(iobject.position)
            

        save_data_aux(cw, objects, self.case_dir, self.case_name, filename=None, wizard=True)
        return True
    
    def check_wizard_attributes(self):
        to_verify = []
        to_verify.append( check_if_float( self.ui_ncd.Bore.text(),'Bore') )
        to_verify.append( check_if_float( self.ui_ncd.crank_radius.text(),'Stroke') )
        to_verify.append( check_if_float( self.ui_ncd.rod_length.text(),'Con-rod Length') )
        to_verify.append( check_if_float( self.ui_ncd.Vol_clearance.text(),'Compression Ratio') )
        to_verify.append( check_if_float( self.ui_ncd.vol_constant.text(),'Constant Tank Volume') )
        if False in to_verify:
            return False
        return True

    def accept_d(self):
        if not len(self.ui_ncd.case_name.text()):
            if self.ui_ncd.new_case_checkBox.isChecked():
                show_message('Please, provide a name for the case')
            else:
                show_message('Please, open a case')
            return
        
        if self.ui_ncd.new_case_checkBox.isChecked():
            self.case_name  = str(self.ui_ncd.case_name.text())
            self.case_dir   = CASES_PATH
            if self.ui_ncd.blank_case_radioButton.isChecked():
                self.case_type  = 2
            else:
                if not self.check_wizard_attributes():
                    return
                if not self.create_case_from_wizard():
                    show_message('An error has ocurred creating the case. Please contact the developer.')
                    return
                self.case_type  = 3
        
        if not self.case_name or not self.case_dir or not self.case_type:
            show_message('An error has ocurred. Please contact the developer.')
            return
        self.accept()
        return

    def cancel(self):
        self.close()
        return

    def open_case(self):
        dialog = QtWidgets.QFileDialog(self)
        dialog.setNameFilter("Python Files (*.py)")
        dialog.setWindowTitle('Open an ICESym-GUI Case')
        dialog.setDirectory(CASES_PATH)
        if dialog.exec_():
            filename = dialog.selectedFiles()[0]
            with open(filename, "r") as f:
                line = f.readline()
                if line=="#### ---- ####\n":
                    pathName =  os.path.dirname(filename)
                    moduleName =  os.path.basename(filename).replace('.py','')
                    self.case_name  = moduleName
                    self.case_dir   = pathName
                    self.case_type  = 1
                    self.ui_ncd.case_name.setText(self.case_name)
                else:
                    show_message('Select a valid case to open')
        return
    
    def use_default_cylinder(self,state):
        self.ui_ncd.Bore.setEnabled(not state)
        self.ui_ncd.crank_radius.setEnabled(not state)
        self.ui_ncd.rod_length.setEnabled(not state)
        self.ui_ncd.Vol_clearance.setEnabled(not state)
        return