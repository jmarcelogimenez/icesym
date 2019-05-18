#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 20:20:06 2019

@author: etekken
"""

import os, copy
from PyQt5 import QtCore, QtGui, QtWidgets
from plotTypeOneWidget_ui import Ui_PlotTypeOneWidget
from utils import show_message, PLOT_ARGUMENTS
from units import UNITS, CONVERSIONS
from numpy import loadtxt, take, array, append, trapz
from GeneralAttributes import GeneralAttributes


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

COMPONENTS_DICT = {}
COMPONENTS_DICT['Cylinders'] = 'cyl'
COMPONENTS_DICT['Tubes'] = 'tube'
COMPONENTS_DICT['Tanks'] = 'tank'
COMPONENTS_DICT['Junctions'] = 'junc'
COMPONENTS_DICT['Globals'] = 'global'

LISTNDOFA   = ['Density', 'Velocity', 'Pressure']
LISTNDOFB   = ['Density', 'Pressure','Temperature']
CYLEXTRAS_DICT = {}
CYLEXTRAS   = ['Convective Heat-Transfer Coeff','Radiactive Heat-Transfer Coeff',\
               'Convective Heat-Transfer Rate','Radiactive Heat-Transfer Rate','Burned Mass Fraction',\
               'Burned Mass Fraction Rate','Mass Flow Rate trought Intake Port','Mass Flow Rate trought Exhaust Port',\
               'Volume','Mass of Fuel','Mass of Air', 'Mass of Residual Gas','Total Heat-Transfer Rate',\
               'Fuel Chemical Energy Release', 'Torque']
CYLEXTRAS_DICT[0] = copy.deepcopy(CYLEXTRAS)
CYLEXTRAS_DICT[1] = copy.deepcopy(CYLEXTRAS)
CYLEXTRAS_DICT[2] = copy.deepcopy(CYLEXTRAS)
CYLEXTRAS_DICT[3] = copy.deepcopy(CYLEXTRAS)
CYLEXTRAS_DICT[2].remove('Volume')
CYLEXTRAS_DICT[2].remove('Torque')
#CYLEXTRAS_DICT[3].remove('Volume')

TANKEXTRAS  = ['Mass Flow Rate','Enthalpy Flow Rate','Mass','Convective Heat-Transfer Rate']

GLOBALS = ['Power Indicated','Power Effective','Torque Indicated','Torque Effective',\
           'IMEP per Cylinder','IMEP Global','FMEP per Cylinder','FMEP Global', \
           'BMEP per Cylinder','BMEP Global','SFC Indicated','SFC Effective','Mechanical Efficiency',\
           'Volumetric Efficiency per Cylinder','Volumetric Efficiency Global',\
           'Fuel Conversion Efficiency Indicated','Fuel Conversion Efficiency Effective']

CYLEXTRA_LINES = {}
CYLEXTRA_LINES['Convective Heat-Transfer Coeff']        = [1,0]
CYLEXTRA_LINES['Radiactive Heat-Transfer Coeff']        = [1,1]
CYLEXTRA_LINES['Convective Heat-Transfer Rate']         = [1,2]
CYLEXTRA_LINES['Radiactive Heat-Transfer Rate']         = [1,3]
CYLEXTRA_LINES['Burned Mass Fraction']                  = [2,0]
CYLEXTRA_LINES['Burned Mass Fraction Rate']             = [2,1]
CYLEXTRA_LINES['Mass Flow Rate trought Intake Port']    = [3,0]
CYLEXTRA_LINES['Mass Flow Rate trought Exhaust Port']   = [3,1]
CYLEXTRA_LINES['Volume']                                = [3,2]
# La linea [3,3] es para la variable Vdot, no se grafica. (Ver codigo fortran)
CYLEXTRA_LINES['Mass of Fuel']                          = [3,4]
CYLEXTRA_LINES['Mass of Air']                           = [3,5]
CYLEXTRA_LINES['Mass of Residual Gas']                  = [3,6]
CYLEXTRA_LINES['Total Heat-Transfer Rate']              = [3,7]
CYLEXTRA_LINES['Fuel Chemical Energy Release']          = [3,8]
CYLEXTRA_LINES['Torque']                                = [3,9]

TANKEXTRA_LINES = {}
TANKEXTRA_LINES['Mass Flow Rate']                   =  0
TANKEXTRA_LINES['Enthalpy Flow Rate']               =  1
TANKEXTRA_LINES['Mass']                             = -2
TANKEXTRA_LINES['Convective Heat-Transfer Rate']    = -1

class PlotTypeOneWidget(QtWidgets.QWidget):
    def __init__(self, plot_function, current_test_dir, current_configuration, current_objects, plot_type, get_oa, set_oa):
        #ptype: plot type. 0 angle, 1 time, 2 RPM, 3 Cycle
        QtWidgets.QWidget.__init__(self)
        self.ui = Ui_PlotTypeOneWidget()
        self.ui.setupUi(self)
        self.current_configuration = current_configuration
        self.current_objects = current_objects
        self.current_test_dir = current_test_dir
        self.plot_function = plot_function
        self.plot_type = plot_type
        # Funciones para obtener y setear los archivos abiertos, existe una 
        # unica instancia compartida por todos los PlotTypeOneWidget y la tiene
        # el padre postProcessWidget
        self.get_open_archives = get_oa
        self.set_open_archives = set_oa
        self.set_restrictions()
        self.set_rpms_and_cycles()
        self.choose_component('Cylinders')        
        return
    
    def set_restrictions(self):
        if self.plot_type in (0,1,3):
            self.ui.component.removeItem(4) # Globals
        elif self.plot_type==2:
            self.ga = GeneralAttributes(self.current_configuration, self.read_normal_txt, \
                                        self.read_extras_txt, self.current_objects['Cylinders'], \
                                        self.current_objects['Atmospheres'][0].object['state_ini'][0], self.current_test_dir)
        if self.plot_type in (1,3):
            self.ui.cycles.setEnabled(False)
        if self.plot_type == 2:
            self.ui.rpms.setEnabled(False)

        self.xlabel             = PLOT_ARGUMENTS[self.plot_type]['xlabel']
        self.xunits             = PLOT_ARGUMENTS[self.plot_type]['xunits']
        self.not_check_cycle    = PLOT_ARGUMENTS[self.plot_type]['ncheck_cycle']
        self.normal_x_var       = PLOT_ARGUMENTS[self.plot_type]['normal_x_var']
        self.extras_x_var       = PLOT_ARGUMENTS[self.plot_type]['extras_x_var']

        self.ui.title.setText(PLOT_ARGUMENTS[self.plot_type]['title'])
        self.ui.legend.setText(PLOT_ARGUMENTS[self.plot_type]['legend'])
        return

    def set_rpms_and_cycles(self):
        self.ui.rpms.setSpacing(3)
        self.ui.cycles.setSpacing(3)
        for irpm in self.current_configuration['rpms']:
            rpm_folder = self.current_test_dir + "/RPM_%s"%irpm
            if not os.path.isdir(rpm_folder):
                show_message("There is no folder for RPM %s. Maybe the simulation is incomplete"%irpm, 2)
                continue
            it = QtWidgets.QListWidgetItem(str(irpm))
            it.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsUserCheckable)
            it.setCheckState(0)
            self.ui.rpms.addItem(it)
            
        last_rpm = self.ui.rpms.item(self.ui.rpms.count()-1)
        if last_rpm:
            last_rpm.setCheckState(2)
        for icycle in range(0,self.current_configuration['ncycles']):
            it = QtWidgets.QListWidgetItem(str(icycle+1))
            it.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsUserCheckable)
            it.setCheckState(0)
            self.ui.cycles.addItem(it)
        last_cycle = self.ui.cycles.item(self.ui.cycles.count()-1)
        if last_cycle:
            last_cycle.setCheckState(2)
        return
    
    def get_list_items(self, list_widget):
        selected = []
        for index_list in range(0,list_widget.count()):
            it = list_widget.item(index_list)
            if it.checkState()==2:
                selected.append(str(it.text()))
        return selected
    
    def clear_comboboxes(self, indexs):
        self.ui.element.clear() if 0 in indexs else None
        self.ui.node.clear() if 1 in indexs else None
        self.ui.variable.clear() if 2 in indexs else None
        self.ui.units.clear() if 3 in indexs else None
        return
    
    def choose_component(self, component):        
        if component=='Globals':
            self.clear_comboboxes([0,1,2,3])
            self.ui.variable.addItems(GLOBALS)
            return 

        self.clear_comboboxes([0,1,2,3])
        selected_rpms = self.get_list_items(self.ui.rpms)
        added_elements = []
        for irpm in selected_rpms:
            rpm_folder = self.current_test_dir + "/RPM_%s"%irpm
            if not os.path.isdir(rpm_folder):
                show_message("There is no folder for RPM %s"%irpm)
                return
            else:
                archives = [f for f in os.listdir(rpm_folder) if (COMPONENTS_DICT[component] \
                            in f and '.txt' in f and 'extras' not in f)]
                archives.sort() # Para que me aparezcan en orden numerico
                for index,iarchive in enumerate(archives):
                    iarchive = iarchive.replace('.txt','')
                    new_element = component[0:-1]+" "+iarchive[-1]
                    if new_element not in added_elements:
                        self.ui.element.addItem(new_element)
                        added_elements.append(new_element)

        self.choose_element(component[0:-1]+" "+str(0))
        return
    
    def choose_element(self, element):
        if not self.ui.element.count():
            return
        self.clear_comboboxes([1,2,3])
        index_element = int(element[-1])
        component = element[0:-2]+'s'
        if index_element>=len(self.current_objects[component]):
            show_message("Error trying to find %s"%element)
            self.clear_comboboxes([1,2,3])
            return
        self.current_index_element = index_element
        iobject = self.current_objects[component][self.current_index_element]
        if 'histo' not in iobject.object.keys():
            show_message("Error trying to find key histo in %s"%element)
            self.clear_comboboxes([1,2,3])
            return
        for ihisto in iobject.object['histo']:
            self.ui.node.addItem(str(ihisto))
        self.choose_variable(element)
        return
    
    def choose_variable(self, element):
        if not self.ui.element.count():
            return
        self.clear_comboboxes([2,3])
        index_element = int(element[-1])
        component = element[0:-2]+'s'
        selected_rpms = self.get_list_items(self.ui.rpms)
        rpm_folder = self.current_test_dir + "/RPM_%s"%selected_rpms[0] # TODO ver bien como hacer esto (verificar que los extras esten en todas las rpm)
        archive_extra = [f for f in os.listdir(rpm_folder) if (COMPONENTS_DICT[component] \
                            in f and '.txt' in f and 'extras' in f and str(index_element) in f)]
        if component in ('Cylinders','Tanks'):
            self.ui.variable.addItems(LISTNDOFB)
        if component in ('Tubes','Junctions'):
            self.ui.variable.addItems(LISTNDOFA)
        if archive_extra!=[]:
            self.ui.variable.addItems(CYLEXTRAS_DICT[self.plot_type] if component=='Cylinders' else\
                                          TANKEXTRAS if component=='Tanks' else None)
        return

    def change_variable(self, variable):
        if variable == '':
            return
        self.ui.units.clear()
        self.ui.units.addItems(UNITS[variable])
        return

    def read_normal_txt(self, archive, node, icycle, variable_index, unit):
        open_archives = self.get_open_archives()
        # Nodo, ciclo, angulo, tiempo .. variables
        if archive not in open_archives.keys():
            A = loadtxt(archive)
            self.set_open_archives(archive,A)
        else:
            A = open_archives[archive]

        scale = CONVERSIONS[unit]
        A_node_filtered = A[A[:,0] == node]
        A_node_and_cycle_filtered = A_node_filtered if self.not_check_cycle \
        else A_node_filtered[A_node_filtered[:,1] == int(icycle)]
        data = take(A_node_and_cycle_filtered, [self.normal_x_var,4+variable_index], axis=1)
        for idata in data:
            idata[1] = idata[1]*scale
        return data

    def loadextratxt(self, archive, offset):
        with open(archive, "r") as f:
            all_data = [map(float,x.split()) for x in f.readlines()]
        A = array([])
        A = append(A,all_data)
        return A

    def read_extras_txt(self, archive, icycle, variable, component, unit):
        offset = 4 if component=='Cylinders' else 2
        variable_line = CYLEXTRA_LINES[variable][0] if component=='Cylinders' else 1
        variable_col = CYLEXTRA_LINES[variable][1] if component=='Cylinders' else TANKEXTRA_LINES[variable]
        scale = CONVERSIONS[unit]
        
        open_archives = self.get_open_archives()
        if archive not in open_archives.keys():
            A = self.loadextratxt(archive, offset)
            self.set_open_archives(archive,A)
        else:
            A = open_archives[archive]
            
        # El extras de cilindro repite un patron cada cuatro lineas.
        # Primer linea: ciclo, angulo, tiempo
        # Segunda, tercera, cuarta linea: varias variables
        # El extras de tanque repite un patron cada 2 lineas.
        # Primer linea: ciclo, angulo, tiempo, ntubos
        # Segunda linea: Mass Flow Rate (ntubos), Entalphy Flow Rate (ntubos), luego dos mas

        if component=='Tanks':
            ntubes = A[0][3]
            TANKEXTRA_LINES['Enthalpy Flow Rate'] = int(ntubes)
            variable_col = TANKEXTRA_LINES[variable]

        ndata = A.shape[0]
        data = [[ A[i][self.extras_x_var], A[i+variable_line][variable_col]*scale ] \
                for i in range(0,ndata,offset) if (int(A[i][0])==int(icycle) or self.not_check_cycle)]
        return data
    
    def verify_data(self, rpms, cycles):
        for ilist in (rpms,cycles):
            if ilist==[]:
                show_message('Please, select at least one RPM and one cycle')
                return False
        cbs = ['element','node']
        for icb in cbs:
            cb = self.ui.__getattribute__(icb)
            if cb.count()==0:
                show_message('There are no items for %s'%(icb.capitalize()))
                return False
        return True
    
    def archive_to_open(self, variable, rpm_folder, component):
        if variable in LISTNDOFA or variable in LISTNDOFB:
            archive = rpm_folder + "/" + COMPONENTS_DICT[component] + "_" + str(self.current_index_element) + ".txt"
            extras = False
        elif variable in CYLEXTRAS_DICT[self.plot_type] or variable in TANKEXTRAS:
            archive = rpm_folder + "/" + COMPONENTS_DICT[component] + "_extras_" + str(self.current_index_element) + ".txt"
            extras = True
        return (archive,extras)

    def trapz_data(self, data):
        value = (data[0][1]+data[-1][1])/2.0
        y = [el[1] for el in data]
        x = [el[0] for el in data]
        y.insert(0,value)
        x.insert(0,0.0)
        max_angle = 720.0 if max(x)>370.0 else 360.0
#        max_angle = 720.0 if x[-1]>370.0 else 360.0 # TODO: bug gui vieja
        x.append(max_angle)
        y.append(value)
        res = trapz(y,x)
        res = res/max_angle
        return res

    def get_plot_attributes(self):
        plot_attributes = {}
        plot_attributes['component'] = str(self.ui.component.currentText())
        plot_attributes['variable'] = str(self.ui.variable.currentText())
        plot_attributes['selected_cycles'] = self.get_list_items(self.ui.cycles) if self.plot_type!=3\
                        else [int(icycle) for icycle in range(1,int(self.current_configuration['ncycles'])+1)]
        plot_attributes['label'] = str(self.ui.legend.text())
        plot_attributes['selected_rpms'] = self.get_list_items(self.ui.rpms) if self.plot_type!=2\
                        else [int(irpm) for irpm in self.current_configuration['rpms']]
        plot_attributes['variable_index'] = int(self.ui.variable.currentIndex())
        plot_attributes['title'] = str(self.ui.title.text())
        plot_attributes['figure_number'] = self.ui.figure_number.currentIndex()-1
        plot_attributes['unit'] = str(self.ui.units.currentText())
        return plot_attributes

    def plot_angle_or_time(self, plot_attributes):
        datas = []
        legends = []
        plot_attributes['node'] = int(self.ui.node.currentText())
        for irpm in plot_attributes['selected_rpms']:
            rpm_folder = self.current_test_dir + "/RPM_%s"%irpm
            for icycle in plot_attributes['selected_cycles']:
                (archive, extras) = self.archive_to_open(plot_attributes['variable'],rpm_folder,plot_attributes['component'])
                try:
                    if extras:
                        data = self.read_extras_txt(archive,icycle,plot_attributes['variable'],plot_attributes['component'],plot_attributes['unit'])
                    else:
                        data = self.read_normal_txt(archive,plot_attributes['node'],icycle,plot_attributes['variable_index'],plot_attributes['unit'])
                    datas.append(data)
                    legends.append(plot_attributes['label']+"_RPM_"+str(irpm)+"_Cycle_"+str(icycle))
                except:
                    show_message('Error opening archive %s. Cannot plot this selections.'%archive)
                    return ([],[])
        return (datas,legends)

    def plot_cycle(self, plot_attributes):
        datas = []
        legends = []
        plot_attributes['node'] = int(self.ui.node.currentText())
        for irpm in plot_attributes['selected_rpms']:
            data_irpm = []
            rpm_folder = self.current_test_dir + "/RPM_%s"%irpm
            for icycle in plot_attributes['selected_cycles']:
                (archive, extras) = self.archive_to_open(plot_attributes['variable'],rpm_folder,plot_attributes['component'])
                try:
                    if extras:
                        data = self.read_extras_txt(archive,icycle,plot_attributes['variable'],plot_attributes['component'],plot_attributes['unit'])
                    else:
                        data = self.read_normal_txt(archive,plot_attributes['node'],icycle,plot_attributes['variable_index'],plot_attributes['unit'])
                    res = self.trapz_data(data)
                    data = [icycle,res]
                    data_irpm.append(data)
                except:
                    show_message('Error opening archive %s. Cannot plot this selections.'%archive)
                    return ([],[])
            datas.append(data_irpm)
            legends.append(plot_attributes['label']+"_RPM_"+str(irpm))
        return (datas,legends)

    def plot_rpm(self, plot_attributes):
        datas = []
        legends = []
        if plot_attributes['component']=='Globals':
            for icycle in plot_attributes['selected_cycles']:
                [data,legend] = self.ga.return_calculated_variable(plot_attributes['variable'],int(icycle),\
                                                                    plot_attributes['label'],plot_attributes['unit'])
                datas.append(data)
                legends.append(legend)
            # flat the lists
            datas = [y for x in datas for y in x]
            legends = [y for x in legends for y in x]
        else:
            plot_attributes['node'] = int(self.ui.node.currentText())
            for icycle in plot_attributes['selected_cycles']:
                data_icycle = []
                for irpm in plot_attributes['selected_rpms']:
                    rpm_folder = self.current_test_dir + "/RPM_%s"%irpm
                    (archive, extras) = self.archive_to_open(plot_attributes['variable'],rpm_folder,plot_attributes['component'])
                    try:
                        if extras:
                            data = self.read_extras_txt(archive,icycle,plot_attributes['variable'],plot_attributes['component'],plot_attributes['unit'])
                        else:
                            data = self.read_normal_txt(archive,plot_attributes['node'],icycle,plot_attributes['variable_index'],plot_attributes['unit'])
                        res = self.trapz_data(data)
                        data = [irpm,res]
                        data_icycle.append(data)
                    except:
                        show_message('Error opening archive %s. Cannot plot this selections.'%archive)
                        return ([],[])
                datas.append(data_icycle)
                legends.append(plot_attributes['label']+"_Cycle_"+str(icycle))
        return (datas,legends)

    def prepare_plot(self):
        self.ui.plot_pushButton.setEnabled(False)
        QtWidgets.QApplication.processEvents()
        try:
            plot_attributes = self.get_plot_attributes()
            datas = []
            legends = []
            if self.plot_type<2:
                (datas,legends) = self.plot_angle_or_time(plot_attributes)
            elif self.plot_type==2:
                (datas,legends) = self.plot_rpm(plot_attributes)
            elif self.plot_type==3:
                (datas,legends) = self.plot_cycle(plot_attributes)
    
            n_plots = self.plot_function(datas, plot_attributes['title'], legends, self.xlabel, plot_attributes['variable'], \
                                         self.xunits, plot_attributes['unit'], plot_attributes['figure_number'], self.plot_type)
            if plot_attributes['figure_number']==-1:
                self.ui.figure_number.addItem('Figure '+str(n_plots-1))
        except:
            show_message('Error trying to plot the current selection')
        self.ui.plot_pushButton.setEnabled(True)
        return