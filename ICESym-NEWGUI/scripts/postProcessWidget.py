#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 00:08:25 2019

@author: etekken
"""

import os, sys, h5py
from PyQt5 import QtCore, QtGui, QtWidgets
from postProcessWidget_ui import Ui_PostProcessWidget
from PlotTypeOneWidget import PlotTypeOneWidget
from PlotTypeTwoWidget import PlotTypeTwoWidget
from PlotTypeThreeWidget import PlotTypeThreeWidget
import pyqtgraph as pg
import numpy as np
from utils import set_plot, show_message, check_two_equals, SelectionDialog, PLOT_ARGUMENTS

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
    
NTYPE_PLOTS = 6

class postProcessWidget(QtWidgets.QWidget):
    def __init__(self, current_dir, current_configuration, current_objects):
        QtWidgets.QWidget.__init__(self)
        self.ui_ppw = Ui_PostProcessWidget()
        self.ui_ppw.setupUi(self)
        self.current_dir = current_dir
        self.current_configuration = current_configuration
        self.current_objects = current_objects
        self.current_test_dir = 'runs/' + self.current_configuration['folder_name']
        self.apw = None
        self.tpw = None
        self.rpw = None
        self.cpw = None
        self.spw = None
        self.fpw = None
        self.colours        = ['r', 'b', 'k', 'g', 'y', 'c']
        self.plots          = {}
        self.legends        = {}
        self.colour_plots   = {}
        self.open_archives  = {}
        self.enable_ppw()
        return

    def enable_ppw(self):
        if not os.path.isdir(self.current_test_dir):
            self.ui_ppw.tabWidget_plots.setEnabled(False)
        else:
            if self.load_header_file() and self.set_plot_widgets():
                self.ui_ppw.tabWidget_plots.setEnabled(True)
        return

    def load_header_file(self):
        filename = self.current_test_dir+'/header.py'
        if not os.path.isfile(filename):
            show_message('There is an error with the header.py file needed to configure the postprocess')
            return False
        pathName =  os.path.dirname(filename)
        moduleName =  os.path.basename(filename).replace('.py','')
        sys.path.append(str(pathName))
        self.configData = __import__(str(moduleName))
        
        # TODO: final_times = rpm/60.0 * ciclos * strokespc/2.0
        
        # Atributo que solo existe en header.py (?) se lo paso para que pueda
        # usarlo la funcion getMasses de GeneralAttributes
        for index,icylinder in enumerate(self.configData.Cylinders):
            self.current_objects['Cylinders'][index].object['angleClose'] = icylinder['angleClose']
            self.current_objects['Cylinders'][index].object['Q_fuel'] = icylinder['Q_fuel']
            
        self.current_configuration['rpms'] = self.configData.rpms
        self.current_configuration['final_times'] = self.configData.final_times
        return True

    def set_plot_widgets(self):
        if not self.apw:
            self.apw = PlotTypeOneWidget(self.plot, self.current_test_dir, self.current_configuration, \
                                         self.current_objects, 0, self.get_open_archives, self.set_open_archives)
            self.ui_ppw.widget_angle_layout.addWidget(self.apw)
            self.apw.setAutoFillBackground(True)
        if not self.tpw:
            self.tpw = PlotTypeOneWidget(self.plot, self.current_test_dir, self.current_configuration, \
                                         self.current_objects, 1, self.get_open_archives, self.set_open_archives)
            self.ui_ppw.widget_time_layout.addWidget(self.tpw)
            self.tpw.setAutoFillBackground(True)
        if not self.rpw:
            self.rpw = PlotTypeOneWidget(self.plot, self.current_test_dir, self.current_configuration, \
                                         self.current_objects, 2, self.get_open_archives, self.set_open_archives)
            self.ui_ppw.widget_rpms_layout.addWidget(self.rpw)
            self.rpw.setAutoFillBackground(True)
        if not self.cpw:
            self.cpw = PlotTypeOneWidget(self.plot, self.current_test_dir, self.current_configuration, \
                                         self.current_objects, 3, self.get_open_archives, self.set_open_archives)
            self.ui_ppw.widget_cycles_layout.addWidget(self.cpw)
            self.cpw.setAutoFillBackground(True)
        if not self.spw:
            self.spw = PlotTypeTwoWidget(self.plot, self.current_test_dir, self.current_configuration, \
                                         self.current_objects, 4, self.get_open_archives, self.set_open_archives)
            self.ui_ppw.widget_space_layout.addWidget(self.spw)
            self.spw.setAutoFillBackground(True)
        if not self.fpw:
            self.fpw = PlotTypeThreeWidget(self.plot, self.current_test_dir, self.current_configuration, \
                                         self.current_objects, 5, self.get_open_archives, self.set_open_archives)
            self.ui_ppw.widget_free_layout.addWidget(self.fpw)
            self.fpw.setAutoFillBackground(True)

        for ip in range(0,NTYPE_PLOTS):
            self.plots[ip]           = []
            self.legends[ip]         = []
            self.colour_plots[ip]    = []
        return True

    def advance_colour(self, current_plot, current_color_dict):
        current_color = current_color_dict[current_plot]
        index = self.colours.index(current_color)
        index = index+1 if index<len(self.colours)-1 else -1
        current_color_dict[current_plot] = self.colours[index]
        return

    def get_open_archives(self):
        return self.open_archives

    def set_open_archives(self,archive,data):
        self.open_archives[archive] = data
        return

    def choose_widgets(self, plot_type):
        if plot_type==0:
            tabWidget   = self.ui_ppw.tabWidget_figures_angles
            plotWidget  = self.apw
        elif plot_type==1:
            tabWidget   = self.ui_ppw.tabWidget_figures_time
            plotWidget  = self.tpw
        elif plot_type==2:
            tabWidget   = self.ui_ppw.tabWidget_figures_rpms
            plotWidget  = self.rpw
        elif plot_type==3:
            tabWidget   = self.ui_ppw.tabWidget_figures_cycles
            plotWidget  = self.cpw
        elif plot_type==4:
            tabWidget   = self.ui_ppw.tabWidget_figures_space
            plotWidget  = self.spw
        elif plot_type==5:
            tabWidget   = self.ui_ppw.tabWidget_figures_free
            plotWidget  = self.fpw
        return (tabWidget,plotWidget)

    def remove_figure(self, remove_plot, plot_type):
        try:
            (tabWidget,plotWidget) = self.choose_widgets(plot_type)
            for index,iplot in enumerate(self.plots[plot_type]):
                if remove_plot == iplot:
                    self.plots[plot_type].remove(iplot)
                    self.colour_plots[plot_type].pop(index)
                    self.legends[plot_type].pop(index)
                    tabWidget.removeTab(index)
                    del iplot
                    break
                
            # acomodar numeracion de tabWidgets que quedan y eliminar del combobox de figuras
            for indextab in range(index,tabWidget.count()):
                tabWidget.setTabText(indextab, 'Figure '+str(int(tabWidget.tabText(indextab)[-1])-1))
            plotWidget.ui.figure_number.removeItem(plotWidget.ui.figure_number.count()-1)
        except:
            show_message('An error has occurred. Cannot delete this figure')
        return
    
    def remove_curve(self, remove_plot, plot_type):
        try:
            (tabWidget,plotWidget) = self.choose_widgets(plot_type)
            for index,iplot in enumerate(self.plots[plot_type]):
                if remove_plot == iplot:
                    plot_items = self.plots[plot_type][index].listDataItems()
                    break
            legends = []
            for ilegenditem in self.legends[plot_type][index].items:
                legends.append( ilegenditem[1].text )
    
            if len(plot_items):
                dialog = QtWidgets.QDialog()
                SDialog = SelectionDialog()
                SDialog.setupUi(dialog, legends)
                ret = dialog.exec_()
                ret = True
                if ret:
                    ncurves = SDialog.listWidget.count()
                    er = 0
                    for row in range(ncurves):
                        item = SDialog.listWidget.item(row)
                        checkState = item.checkState()
                        if checkState == 2:
                            self.plots[plot_type][index].removeItem(self.plots[plot_type][index].listDataItems()[row-er])
                            self.legends[plot_type][index].removeItem(self.legends[plot_type][index].items[row-er][1].text)
                            er+=1
                    show_message('Curve(s) successfully erased!',1)
                else:
                    show_message('Operation cancelled',1)
            else:
                show_message('This plot has not curves',1)
        except:
            show_message('An error has occurred. Cannot delete this curves')
        return

    def plot(self, current_data, title, legend_texts, xlabel, ylabel, xunits, yunits, figure_number, plot_type):
        
        # no permitir dos leyendas iguales en un misma plot, puesto que luego para borrar
        # una curva pyqtgraph usa el nombre
        legends_to_check = legend_texts
        if figure_number>=0:
            for item in self.legends[plot_type][figure_number].items:
                legends_to_check.append(item[1].text)
        if check_two_equals(legend_texts):
            show_message('There are two legends with the same name. Please, select another legend name for this selections')
            return
        
        pg.setConfigOptions(background='w')
        pg.setConfigOptions(foreground='k')
        added = False
        (tabWidget,plotWidget) = self.choose_widgets(plot_type)

        for index,idata in enumerate(current_data):

            xdata = []
            ydata = []

            for iidata in idata:
                xdata.append(iidata[0])
                ydata.append(iidata[1])

            if figure_number==-1 and not added:
                new_plot = pg.PlotWidget()
                set_plot(new_plot, xlabel, ylabel, title, xunits, yunits)
                new_plot.getPlotItem().getViewBox().menu.addAction("Delete Curve",  lambda: self.remove_curve(new_plot,plot_type))
                new_plot.getPlotItem().getViewBox().menu.addAction("Delete Figure", lambda: self.remove_figure(new_plot,plot_type))                
                tab = QtWidgets.QWidget()
                tabWidget.addTab(tab,'Figure '+str(len(self.plots[plot_type])))
                tab = tabWidget.widget(len(self.plots[plot_type]))
                self.plots[plot_type].append(new_plot)
                self.colour_plots[plot_type].append(self.colours[0])
                tab.setLayout(QtWidgets.QHBoxLayout())
                tab.layout().addWidget(new_plot)
                legend = pg.LegendItem(offset=(0,1))
                legend.setParentItem(self.plots[plot_type][-1].getPlotItem())
                self.legends[plot_type].append(legend)
                added = True

            it = self.plots[plot_type][figure_number].plot(xdata, ydata,\
                                 pen={'color': self.colour_plots[plot_type][figure_number], 'width': 1})
            self.advance_colour(len(self.plots[plot_type])-1,self.colour_plots[plot_type])
            self.legends[plot_type][figure_number].addItem(it,legend_texts[index])
        return len(self.plots[plot_type])

    def save_postpro(self):
        name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Post Process As', "./", "Hierarchical Data Format Files (*.hdf)")
        filename = name[ 0 ]
        filename = filename+'.hdf' if filename.find('.hdf')==-1 else filename
        try:
            with h5py.File(filename, 'w') as outfile:
                for itypeplot in range(0,NTYPE_PLOTS):
                    for iplot in range(len(self.plots[itypeplot])):
                        # cada curva del plot, en formato (x,y)
                        plot_items = self.plots[itypeplot][iplot].getPlotItem().listDataItems()
                        plot_data = None
                        legends = []
                        for index,idataitem in enumerate(plot_items):
                            data = idataitem.getData()
                            data = [ [data[0][i],data[1][i]] for i in range(len(data[0])) ]
                            data = np.array(data)
                            plot_data = np.concatenate((plot_data, data), axis=0) if plot_data is not None else data
                            legends.append( self.legends[itypeplot][iplot].items[index][1].text )
                        plot_dataset = outfile.create_dataset(PLOT_ARGUMENTS[itypeplot]['title']+' '+str(iplot), data=plot_data)
                        plot_dataset.attrs['title']     = str(self.plots[itypeplot][iplot].getPlotItem().titleLabel.text)
                        plot_dataset.attrs['xlabel']    = str(self.plots[itypeplot][iplot].getPlotItem().getAxis('bottom').labelText)
                        plot_dataset.attrs['xunits']    = str(self.plots[itypeplot][iplot].getPlotItem().getAxis('bottom').labelUnits)
                        plot_dataset.attrs['ylabel']    = str(self.plots[itypeplot][iplot].getPlotItem().getAxis('left').labelText)
                        plot_dataset.attrs['yunits']    = str(self.plots[itypeplot][iplot].getPlotItem().getAxis('left').labelUnits)
                        plot_dataset.attrs['legends']   = legends
                        plot_dataset.attrs['nplots']    = len(plot_items)
            show_message('Post Process successfully saved!',1)
        except:
            show_message('Error saving the archive %s'%filename)
        return

    def exists_plots(self):
        for itp in range(0,NTYPE_PLOTS):
            if len(self.plots[itp]):
                return True
        return False

    def load_postpro(self):
        dialog = QtWidgets.QFileDialog(self)
        dialog.setNameFilter("Hierarchical Data Format Files (*.hdf)")
        dialog.setWindowTitle('Open a Post Process File')
        dialog.setDirectory(self.current_dir)
        if dialog.exec_():
            filename = dialog.selectedFiles()[0]            
            if self.exists_plots():
                msg = "Do you want to preserve all the current figures?"
                reply = show_message(msg,4,QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No)
                if reply == QtWidgets.QMessageBox.No:
                    for itp in range(0,NTYPE_PLOTS):
                        for ip in self.plots[itp]:
                            self.remove_figure(ip,itp)
            try:
                with h5py.File(filename, 'r') as openfile:
                    for ikey in openfile.keys():
                        plot_dataset    = openfile[ikey]
                        prange          = plot_dataset.shape[0]/plot_dataset.attrs['nplots']
                        datas           = []
                        legends         = []
                        for icurve in range(plot_dataset.attrs['nplots']):
                            iplot = []
                            for ip in range(prange):
                                iplot.append( plot_dataset[ip+(icurve*prange)] )
                            datas.append(iplot)
                            legends.append( plot_dataset.attrs['legends'][icurve] )
                        title = plot_dataset.attrs['title']
                        plot_type = 0 if 'Angle' in ikey else 1 if 'Time' in ikey else 2 if \
                                        'RPM' in ikey else 3 if 'Cycle' in ikey else 4 if 'Space' in ikey else 5
                        self.plot(datas, title, legends, plot_dataset.attrs['xlabel'],\
                                  plot_dataset.attrs['ylabel'], plot_dataset.attrs['xunits'], plot_dataset.attrs['yunits'], -1, plot_type)
                        (tabWidget,plotWidget) = self.choose_widgets(plot_type)
                        plotWidget.ui.figure_number.addItem('Figure '+str(len(self.plots[plot_type])-1))
                show_message('Post Process successfully loaded!',1)
            except:
                show_message('Error opening the archive %s'%filename)
        return
        