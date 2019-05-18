#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 23:11:30 2018

@author: etekken
"""

import os
import json
from PyQt5 import QtGui, QtWidgets
from atmosphereDialog_ui import Ui_AtmosphereDialog
from utils import show_message

JSON_ATMOSPHERE_KEYS = ['nnod','ndof','state_ini']

class AtmosphereDialog(QtWidgets.QDialog):
    """
    class to manage the atmosphere atributes. If current_atmosphere is None, we are creating
    a new one. On the other hand, we are modifying an old one.
    """
    def __init__(self, current_dir, current_atmosphere = None, item_index = 0):
        QtWidgets.QDialog.__init__(self)
        self.ui_ad = Ui_AtmosphereDialog()
        self.ui_ad.setupUi(self)
        self.setFixedSize(300, 260)
        self.set_restrictions()
        self.current_dict = None
        self.current_dir = current_dir
        self.setWindowTitle( self.windowTitle() + " " + str(item_index) )
        if not current_atmosphere:
            self.load_default()
        else:
            self.current_dict = current_atmosphere
            self.set_parameters()
#        self.save_valve_template()
        return
    
    def set_restrictions(self):
        """
        set the initial restrictions on lineEdits
        """
        self.ui_ad.density_lineEdit.setValidator(QtGui.QDoubleValidator(0, 360, 3))
        self.ui_ad.velocity_lineEdit.setValidator(QtGui.QDoubleValidator(0, 360, 3))
        self.ui_ad.pressure_lineEdit.setValidator(QtGui.QDoubleValidator(0, 100, 3))
        return
    
    def set_parameters(self):
        """
        puts the loaded atmosphere parameters in the items
        """
        if not self.current_dict:
            return
        self.ui_ad.density_lineEdit.setText(str(self.current_dict['state_ini'][0]))
        self.ui_ad.velocity_lineEdit.setText(str(self.current_dict['state_ini'][1]))
        self.ui_ad.pressure_lineEdit.setText(str(self.current_dict['state_ini'][2]))
        return
    
    def load_default(self):
        """
        load the default atmosphere template
        """
        filename = self.current_dir + "templates/atmosphere_default.json"
        if not os.path.isfile(filename):
            show_message("Cannot find default atmosphere configuration file")
            return
        with open(filename) as openedfile:
            try:
                self.current_dict = json.load(openedfile)
            except ValueError, error:
                print 'JSON object issue: %s'%error
        if self.current_dict is None:
            show_message("An error ocurred with the json default atmosphere archive")
            return
        if not self.check_json_keys():
            return
        self.set_parameters()
        return
    
    def check_json_keys(self):
        """
        check if the loaded json include all the obligatory keys
        """
        for ikey in JSON_ATMOSPHERE_KEYS:
            if ikey not in self.current_dict.keys():
                show_message("Wrong number of keys in json default atmosphere archive")
                return False
        return True
    
    def accept(self):
        success = False
        try:        
            self.current_dict['state_ini'][0] = float( self.ui_ad.density_lineEdit.text() )
            self.current_dict['state_ini'][1] = float( self.ui_ad.velocity_lineEdit.text() )
            self.current_dict['state_ini'][2] = float( self.ui_ad.pressure_lineEdit.text() )
            success = True
        except:
            show_message('Please check that all the fields has correct values')
        if success:
            self.close()
        return
    
    def cancel(self):
        self.close()
        return