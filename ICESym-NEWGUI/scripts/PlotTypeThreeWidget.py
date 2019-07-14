#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 20:20:06 2019

@author: etekken
"""

import os
from PyQt5 import QtCore, QtGui, QtWidgets
from plotTypeThreeWidget_ui import Ui_PlotTypeThreeWidget
from utils import show_message
from numpy import loadtxt, take, array, append
from units import UNITS, CONVERSIONS
from exception_handling import handle_exception

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
COMPONENTS_DICT['Tanks'] = 'tank'

LISTNDOFA   = ['Density', 'Velocity', 'Pressure']
LISTNDOFB   = ['Density', 'Pressure','Temperature']
CYLEXTRAS_DICT = {}
CYLEXTRAS   = ['Convective Heat-Transfer Coeff','Radiactive Heat-Transfer Coeff',\
               'Convective Heat-Transfer Rate','Radiactive Heat-Transfer Rate','Burned Mass Fraction',\
               'Burned Mass Fraction Rate','Mass Flow Rate trought Intake Port','Mass Flow Rate trought Exhaust Port',\
               'Volume','Mass of Fuel','Mass of Air', 'Mass of Residual Gas','Total Heat-Transfer Rate',\
               'Fuel Chemical Energy Release', 'Torque']

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

class PlotTypeThreeWidget(QtWidgets.QWidget):
    def __init__(self, plot_function, current_run_dir, run_attributes, current_objects, plot_type, get_oa, set_oa):
        QtWidgets.QWidget.__init__(self)
        self.ui = Ui_PlotTypeThreeWidget()
        self.ui.setupUi(self)
        self.current_run_dir = current_run_dir
        self.plot_function = plot_function
        self.plot_type = plot_type
        self.change_attributes(run_attributes, current_objects)
        # Funciones para obtener y setear los archivos abiertos, existe una 
        # unica instancia compartida por todos los PlotTypeOneWidget y la tiene
        # el padre postProcessWidget
        self.get_open_archives = get_oa
        self.set_open_archives = set_oa
        self.set_restrictions()
        return

    def change_attributes(self, run_attributes, current_objects):
        self.run_attributes = run_attributes
        self.current_objects = current_objects
        self.set_rpms_and_cycles()
        self.choose_component('Cylinders')
        return

    def set_restrictions(self):
        self.not_check_cycle = False
        self.ui.title.setText('Free Plot')
        self.ui.legend.setText('Free Legend')
        return

    def set_rpms_and_cycles(self):
        self.ui.rpms.setSpacing(3)
        self.ui.cycles.setSpacing(3)
        self.ui.rpms.clear()
        self.ui.cycles.clear()
        for irpm in self.run_attributes['rpms']:
            it = QtWidgets.QListWidgetItem(str(irpm))
            it.setFlags(QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsUserCheckable)
            it.setCheckState(0)
            self.ui.rpms.addItem(it)
            
        last_rpm = self.ui.rpms.item(self.ui.rpms.count()-1)
        if last_rpm:
            last_rpm.setCheckState(2)
        for icycle in range(0,self.run_attributes['ncycles']):
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
        self.ui.x_variable.clear() if 2 in indexs else None
        self.ui.y_variable.clear() if 2 in indexs else None
        self.ui.x_units.clear() if 3 in indexs else None
        self.ui.y_units.clear() if 3 in indexs else None
        return

    def choose_component(self, component):
        self.clear_comboboxes([0,1,2,3])
        selected_rpms = self.get_list_items(self.ui.rpms)
        added_elements = []
        for irpm in selected_rpms:
            rpm_folder = os.path.join(self.current_run_dir,"RPM_%s"%irpm)
            if not os.path.isdir(rpm_folder):
                show_message("There is no folder for RPM %s"%irpm)
                return
            else:
                archives = [f for f in os.listdir(rpm_folder) if (COMPONENTS_DICT[component] \
                            in f and '.txt' in f and 'extras' not in f)]
                archives.sort() # Para que me aparezcan en orden numerico
                for index,iarchive in enumerate(archives):
                    iarchive = iarchive.replace('.txt','')
                    new_element = self.current_objects[component][int(iarchive[-1])].object['label']
                    if new_element not in added_elements:
                        self.ui.element.addItem(new_element)
                        added_elements.append(new_element)

        self.choose_element(self.ui.element.currentText())
        return

    def choose_variable(self, element):
        if not self.ui.element.count():
            return
        self.clear_comboboxes([2,3])
        index_element = self.ui.element.findText(element,QtCore.Qt.MatchFixedString)
        component = self.ui.component.currentText()
        selected_rpms = self.get_list_items(self.ui.rpms)
        rpm_folder = os.path.join(self.current_run_dir,"RPM_%s"%selected_rpms[0])
        archive_extra = [f for f in os.listdir(rpm_folder) if (COMPONENTS_DICT[component] \
                            in f and '.txt' in f and 'extras' in f and str(index_element) in f)]
        if component in ('Cylinders','Tanks'):
            self.ui.x_variable.addItems(LISTNDOFB)
            self.ui.y_variable.addItems(LISTNDOFB)
        if component in ('Tubes','Junctions'):
            self.ui.x_variable.addItems(LISTNDOFA)
            self.ui.y_variable.addItems(LISTNDOFA)
        if archive_extra!=[]:
            self.ui.x_variable.addItems(CYLEXTRAS if component=='Cylinders' else\
                                          TANKEXTRAS if component=='Tanks' else None)
            self.ui.y_variable.addItems(CYLEXTRAS if component=='Cylinders' else\
                                          TANKEXTRAS if component=='Tanks' else None)
        return

    def choose_element(self, element):
        if not self.ui.element.count():
            return
        self.clear_comboboxes([1,2,3])
        index_element = self.ui.element.findText(element,QtCore.Qt.MatchFixedString)
        component = self.ui.component.currentText()
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
        data = take(A_node_and_cycle_filtered, 4+variable_index, axis=1)
        data = [i*scale for i in data]
        
        try:
            assert(len(data) != 0)
        except:
            handle_exception('Cannot find data in %s archive for %s cycle'%(archive,icycle))
        return data

    def loadextratxt(self, archive, offset):
        with open(archive, "r") as f:
            all_data = [list(map(float,x.split())) for x in f.readlines()]
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
        data = [ A[i+variable_line][variable_col]*scale \
                for i in range(0,ndata,offset) if (int(A[i][0])==int(icycle) or self.not_check_cycle)]
        
        try:
            assert(len(data) != 0)
        except:
            handle_exception('Cannot find data in %s archive for %s cycle'%(archive,icycle))
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

    def change_x_variable(self, variable):
        if variable == '':
            return
        self.ui.x_units.clear()
        self.ui.x_units.addItems(UNITS[variable])
        return

    def change_y_variable(self, variable):
        if variable == '':
            return
        self.ui.y_units.clear()
        self.ui.y_units.addItems(UNITS[variable])
        return

    def archive_to_open(self, variable, rpm_folder, component):
        if variable in LISTNDOFA or variable in LISTNDOFB:
            archive = os.path.join(rpm_folder,COMPONENTS_DICT[component]+"_"+str(self.current_index_element)+".txt")
            extras = False
        elif variable in CYLEXTRAS or variable in TANKEXTRAS:
            archive = os.path.join(rpm_folder,COMPONENTS_DICT[component]+"_extras_"+str(self.current_index_element)+".txt")
            extras = True
        return (archive,extras)

    def get_plot_attributes(self):
        plot_attributes = {}
        plot_attributes['component'] = str(self.ui.component.currentText())
        plot_attributes['variable'] = [str(self.ui.x_variable.currentText()),str(self.ui.y_variable.currentText())]
        plot_attributes['selected_cycles'] = self.get_list_items(self.ui.cycles) if self.plot_type!=3\
                        else [int(icycle) for icycle in range(1,int(self.run_attributes['ncycles'])+1)]
        plot_attributes['label'] = str(self.ui.legend.text())
        plot_attributes['selected_rpms'] = self.get_list_items(self.ui.rpms) if self.plot_type!=2\
                        else [int(irpm) for irpm in self.run_attributes['rpms']]
        plot_attributes['variable_index'] = [int(self.ui.x_variable.currentIndex()),int(self.ui.y_variable.currentIndex())]
        plot_attributes['title'] = str(self.ui.title.text())
        plot_attributes['figure_number'] = self.ui.figure_number.currentIndex()-1
        plot_attributes['units'] = [str(self.ui.x_units.currentText()),str(self.ui.y_units.currentText())]
        return plot_attributes

    def plot_free(self, plot_attributes):
        datas = []
        legends = []
        plot_attributes['node'] = int(self.ui.node.currentText())
        for irpm in plot_attributes['selected_rpms']:
            rpm_folder = os.path.join(self.current_run_dir,"RPM_%s"%irpm)
            for icycle in plot_attributes['selected_cycles']:
                (archive_x, extras_x) = self.archive_to_open(plot_attributes['variable'][0],rpm_folder,plot_attributes['component'])
                (archive_y, extras_y) = self.archive_to_open(plot_attributes['variable'][1],rpm_folder,plot_attributes['component'])
                archives = [archive_x,archive_y]
                extras   = [extras_x,extras_y]
                data_p   = []
                for index,iarchive in enumerate(archives):
                    try:
                        if extras[index]:
                            data = self.read_extras_txt(iarchive,icycle,plot_attributes['variable'][index],plot_attributes['component'],plot_attributes['units'][index])
                        else:
                            data = self.read_normal_txt(iarchive,plot_attributes['node'],icycle,plot_attributes['variable_index'][index],plot_attributes['units'][index])
                        data_p.append(data)
                        from exception_handling import CURRENT_EXCEPTION
                        assert(not CURRENT_EXCEPTION)
                    except:
                        handle_exception('Error opening archive %s. Cannot plot this selections.'%iarchive)
                        return
                data_len = min(len(data_p[0]),len(data_p[1])) # sometimes arrays differs in length
                data_p = [[data_p[0][i],data_p[1][i]] for i in range(data_len)]
                datas.append(data_p)
                legends.append(plot_attributes['label']+"_RPM_"+str(irpm)+"_Cycle_"+str(icycle))
        return (datas,legends)

    def prepare_plot(self, _plot_attributes=None):
        self.ui.plot_pushButton.setEnabled(False)
        QtWidgets.QApplication.processEvents()
        try:
            plot_attributes = _plot_attributes if _plot_attributes else self.get_plot_attributes()
            datas = []
            legends = []
            (datas,legends) = self.plot_free(plot_attributes)
            n_plots = self.plot_function(datas, plot_attributes['title'], legends, plot_attributes['variable'][0], plot_attributes['variable'][1],\
                                         plot_attributes['units'][0], plot_attributes['units'][1], plot_attributes['figure_number'], self.plot_type)
            if plot_attributes['figure_number']==-1:
                self.ui.figure_number.addItem('Figure '+str(n_plots-1))
            from exception_handling import CURRENT_EXCEPTION
            assert(not CURRENT_EXCEPTION)
        except:
            handle_exception('Error trying to plot the current selection')
        finally:
            self.ui.plot_pushButton.setEnabled(True)
        return