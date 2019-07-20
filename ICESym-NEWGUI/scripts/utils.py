#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 11 19:54:58 2018

@author: etekken
"""

import math, os, json
from sys import platform
from PyQt5 import QtGui, QtWidgets, QtCore

# conversion constants
RAD2DEG     = 180.0/math.pi
DEG2RAD     = math.pi/180.0
M2MM        = 1e3
SQM2SQCM    = 1e4
CM2CCM      = 1e6
MM2M        = 1e-3
CCM2CM      = 1e-6
SQCM2SQM    = 1e-4

# images paths
VALVE_EXH_PATH  = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','valve_exhaust_color.svg')
VALVE_INT_PATH  = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','valve_intake_color.svg')
ATMOS_PATH      = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','atmosphere_color.svg')
TUBE_PATH       = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','tube_color.svg')
TANK_PATH       = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','tank_color.svg')
JUNCTION_PATH   = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','Junction_color.svg')
CYLINDER_PATH   = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','cylinder_color.svg')
VALVE_EXH_PATH_NC  = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','valve_exhaust_thumb.svg')
VALVE_INT_PATH_NC  = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','valve_intake_thumb.svg')
ATMOS_PATH_NC      = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','atmosphere_thumb.svg')
TUBE_PATH_NC       = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','tube_thumb.svg')
TANK_PATH_NC       = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','tank_thumb.svg')
JUNCTION_PATH_NC   = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','Junction_thumb.svg')
CYLINDER_PATH_NC   = os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','cylinder_thumb.svg')

# paths, for all the cases the same
SIMULATOR_PATH = os.path.join(os.environ["ICESYM_INST_DIR"],'simulator')
RUNS_PATH      = os.path.join(os.environ["ICESYM_INST_DIR"],'runs')
CASES_PATH     = os.path.join(os.environ["ICESYM_INST_DIR"],'cases')
LOADS_PATH     = os.path.join(os.environ["ICESYM_INST_DIR"],'loads')
DOCS_PATH      = os.path.join(os.environ["ICESYM_INST_DIR"],'docs')
TEMPLATES_PATH = os.path.join(os.environ["ICESYM_INST_DIR"],'templates')
INSTALL_PATH   = os.environ["ICESYM_INST_DIR"]

ICON_PATHS = {}
ICON_PATHS['Valves']         = VALVE_INT_PATH
ICON_PATHS['Atmospheres']    = ATMOS_PATH
ICON_PATHS['Tubes']          = TUBE_PATH
ICON_PATHS['Tanks']          = TANK_PATH
ICON_PATHS['Junctions']      = JUNCTION_PATH
ICON_PATHS['Cylinders']      = CYLINDER_PATH
ICON_PATHS['Valve_exh']      = VALVE_EXH_PATH
ICON_PATHS['Valve_int']      = VALVE_INT_PATH

ICON_PATHS_NC = {}
ICON_PATHS_NC['Valves']         = VALVE_INT_PATH_NC
ICON_PATHS_NC['Atmospheres']    = ATMOS_PATH_NC
ICON_PATHS_NC['Tubes']          = TUBE_PATH_NC
ICON_PATHS_NC['Tanks']          = TANK_PATH_NC
ICON_PATHS_NC['Junctions']      = JUNCTION_PATH_NC
ICON_PATHS_NC['Cylinders']      = CYLINDER_PATH_NC
ICON_PATHS_NC['Valve_exh']      = VALVE_EXH_PATH_NC
ICON_PATHS_NC['Valve_int']      = VALVE_INT_PATH_NC

# default density, velocity and pressure
DEFAULT_DVP = [1.1769, 0.1, 101330.0]

PLOT_ARGUMENTS = {}
PLOT_ARGUMENTS[0] = {}
PLOT_ARGUMENTS[1] = {}
PLOT_ARGUMENTS[2] = {}
PLOT_ARGUMENTS[3] = {}
PLOT_ARGUMENTS[4] = {}
PLOT_ARGUMENTS[5] = {}
PLOT_ARGUMENTS[0]['title']          = 'Angle Plot'
PLOT_ARGUMENTS[0]['legend']         = 'Angle Legend'
PLOT_ARGUMENTS[0]['xlabel']         = 'Angle'
PLOT_ARGUMENTS[0]['xunits']         = 'deg'
PLOT_ARGUMENTS[0]['ncheck_cycle']   = False
PLOT_ARGUMENTS[0]['normal_x_var']   = 2
PLOT_ARGUMENTS[0]['extras_x_var']   = 1
PLOT_ARGUMENTS[1]['title']          = 'Time Plot'
PLOT_ARGUMENTS[1]['legend']         = 'Time Legend'
PLOT_ARGUMENTS[1]['xlabel']         = 'Time'
PLOT_ARGUMENTS[1]['xunits']         = 'sec'
PLOT_ARGUMENTS[1]['ncheck_cycle']   = True
PLOT_ARGUMENTS[1]['normal_x_var']   = 3
PLOT_ARGUMENTS[1]['extras_x_var']   = 2
PLOT_ARGUMENTS[2]['title']          = 'RPM Plot'
PLOT_ARGUMENTS[2]['legend']         = 'RPM Legend'
PLOT_ARGUMENTS[2]['xlabel']         = 'RPM'
PLOT_ARGUMENTS[2]['xunits']         = ''
PLOT_ARGUMENTS[2]['ncheck_cycle']   = False
PLOT_ARGUMENTS[2]['normal_x_var']   = 2
PLOT_ARGUMENTS[2]['extras_x_var']   = 1
PLOT_ARGUMENTS[3]['title']          = 'Cycle Plot'
PLOT_ARGUMENTS[3]['legend']         = 'Cycle Legend'
PLOT_ARGUMENTS[3]['xlabel']         = 'Cycle'
PLOT_ARGUMENTS[3]['xunits']         = ''
PLOT_ARGUMENTS[3]['ncheck_cycle']   = False
PLOT_ARGUMENTS[3]['normal_x_var']   = 2
PLOT_ARGUMENTS[3]['extras_x_var']   = 1
PLOT_ARGUMENTS[4]['title']          = 'Space Plot'
PLOT_ARGUMENTS[5]['title']          = 'Free Plot'

# default plots to graph when the simulation finish
DEFAULT_PLOTS = []
DEFAULT_PLOTS.append([5,'Cylinders',0,0,['Volume','Pressure'],['m^3','Pa'],[-1],[-1],'PV_Cyl0','PV Plot Cylinder 0',-1])
DEFAULT_PLOTS.append([0,'Cylinders',0,-1,'Temperature','K',[-1],[-1],'Temperature','Cylinder Temperatures',-1])
DEFAULT_PLOTS.append([0,'Cylinders',0,0,'Mass Flow Rate trought Intake Port','kg/s',[-1],[0],'Mass Flow Rate trought Intake Port','Mass Flow Rate (Intake)',-1])
DEFAULT_PLOTS.append([0,'Cylinders',0,0,'Mass Flow Rate trought Intake Port','kg/s',[-1],[-1],'Mass Flow Rate trought Intake Port','Mass Flow Rate (Intake)',0])
DEFAULT_PLOTS.append([0,'Cylinders',0,0,'Mass Flow Rate trought Exhaust Port','kg/s',[-1],[0],'Mass Flow Rate trought Exhaust Port','Mass Flow Rate (Exhaust)',-1])
DEFAULT_PLOTS.append([0,'Cylinders',0,0,'Mass Flow Rate trought Exhaust Port','kg/s',[-1],[-1],'Mass Flow Rate trought Exhaust Port','Mass Flow Rate (Exhaust)',0])
DEFAULT_PLOTS.append([2,'Globals',0,0,'Power Indicated','W',[-1],[],'Power Indicated','Power',-1])
DEFAULT_PLOTS.append([2,'Globals',0,0,'Power Effective','W',[-1],[],'Power Effective','Power',0])
DEFAULT_PLOTS.append([2,'Globals',0,0,'Torque Indicated','N.m',[-1],[],'Torque Indicated','Torque',-1])
DEFAULT_PLOTS.append([2,'Globals',0,0,'Torque Effective','N.m',[-1],[],'Torque Effective','Torque',0])

# tab's information
INDEX_TAB_CONFIGURATION = 0
INDEX_TAB_MODELING      = 1
INDEX_TAB_RUN           = 2
INDEX_TAB_POSTPROCESS   = 3
TAB_INFORMATION = {}
TAB_INFORMATION[INDEX_TAB_CONFIGURATION]    = {}
TAB_INFORMATION[INDEX_TAB_MODELING]         = {}
TAB_INFORMATION[INDEX_TAB_RUN]              = {}
TAB_INFORMATION[INDEX_TAB_POSTPROCESS]      = {}
TAB_INFORMATION[INDEX_TAB_CONFIGURATION]['name']    = 'Configuration'
TAB_INFORMATION[INDEX_TAB_MODELING]['name']         = 'Modeling'
TAB_INFORMATION[INDEX_TAB_RUN]['name']              = 'Run'
TAB_INFORMATION[INDEX_TAB_POSTPROCESS]['name']      = 'Post Process'

# curve's format information
CURVE_COLORS = {}
CURVE_COLORS['Red']      = 'r'
CURVE_COLORS['Blue']     = 'b'
CURVE_COLORS['Black']    = 'k'
CURVE_COLORS['Green']    = 'g'
CURVE_COLORS['Yellow']   = 'y'
CURVE_COLORS['Cyan']     = 'c'
CURVE_COLORS['Magenta']  = 'm'
CURVE_LINE_FORMATS = {}
CURVE_LINE_FORMATS['Solid']         = QtCore.Qt.SolidLine
CURVE_LINE_FORMATS['Dash']          = QtCore.Qt.DashLine
CURVE_LINE_FORMATS['Dot']           = QtCore.Qt.DotLine
CURVE_LINE_FORMATS['Dash Dot']      = QtCore.Qt.DashDotLine
CURVE_LINE_FORMATS['Dash Dot Dot']  = QtCore.Qt.DashDotDotLine

def set_plot(plot, leg_bottom, leg_left, title, units_x = '', units_y = ''):
    """
    set the parameters to the given plot
    """
    labelStyle = {'color': '#000', 'type': 'Ubuntu', 'size': '9pt'}
    plot.setTitle(title, **labelStyle)
    labelStyle = {'color': '#000', 'font-type': 'Ubuntu', 'font-size': '9pt'}
    plot.setLabel('bottom', leg_bottom, units_x, **labelStyle)
    plot.setLabel('left', leg_left, units_y, **labelStyle)
    bottomaxis = plot.getAxis('bottom')
    leftaxis = plot.getAxis('left')
    leftaxis.enableAutoSIPrefix(False)
    bottomaxis.enableAutoSIPrefix(False)
    font = QtGui.QFont("Ubuntu", 9)
    bottomaxis.setFont(font)
    leftaxis.setFont(font)
    leftaxis.tickFont = font
    bottomaxis.tickFont = font
    plot.showGrid(True, True)
    return

def remove_values_from_list(the_list, val):
    """
    Remueve valores en los elementos de una lista dada.

    Keyword arguments:
    the_list -- la lista a procesar
    val -- el valor a eliminar en los elementos de la lista
    Return values:
    the_list -- la lista procesada
    """
    the_list = [element.replace(val, '') for element in the_list]
    return the_list

def convert_string(text):
    """
    castea al string al tipo de dato numerico correspondiente, de otro modo 
    devuelve el string sin modificar
    """
    try:
        text = int(text)
    except:
        try:
            text = float(text)
        except:
            None
    return text

def show_errors(currentFolder):
    """
    Chequear si hubo algun error e imprimirlo en un mensaje
    """
    filename = os.path.join(currentFolder,'error.log')
    #Si fue creado y tiene algo escrito lo muestro y lo muevo
    if os.path.isfile(filename) and os.path.getsize(filename) > 0:
        with open(filename, 'r') as log:
            content = log.readlines()
            while '\n' in content:
                content.remove('\n')
            content = ''.join(content)
            if 'Killed\n'==content: # Para que al matar no salte error
                return True

            show_message(content)
            log.close()
            return True
    elif os.path.isfile(filename): #Si fue creado pero esta vacio lo borro
        if 'linux' in platform:
            command = 'rm %s '%filename
        elif 'win' in platform:
            command = 'del /f %s '%filename
        os.system(command)
    return False

def check_if_float(text, argument):
    try:
        float(text)
        return True
    except:
        show_message('Problem with argument %s'%argument)
        return False
    
def get_messages_font():
    font = QtGui.QFont()
    font.setFamily("Ubuntu")
    font.setPointSize(10)
    return font
    
def show_message(content, message_type = 3, buttons = QtWidgets.QMessageBox.Ok):
    # 1 Information, 2 Warning, 3 Critical, 4 Question
    font = get_messages_font()
    mtitle = 'Information' if message_type==1 else 'Warning' if message_type==2 else 'Error' if message_type==3 else 'Question'
    w = QtWidgets.QMessageBox(message_type,mtitle,content,buttons)
    w.setFont(font)
    w.setWindowIcon(QtGui.QIcon(os.path.join(os.environ["ICESYM_INST_DIR"],'images','newicons','icon.png')))
    QtWidgets.QApplication.processEvents()
    reply = w.exec_()
    return reply

def load_templates():
    """
    load the default templates
    """
    default_dict = {}
    default_dict['Valves']         = None
    default_dict['Tubes']          = None
    default_dict['Atmospheres']    = None
    default_dict['Junctions']      = None
    default_dict['Cylinders']      = None
    default_dict['Tanks']          = None
    default_dict['Configurations'] = None
    items = ['valve','atmosphere','cylinder','junction','tube','tank','configuration']
    for item in items:
        filename = os.path.join(INSTALL_PATH,"templates","%s_default.json"%item)
        if not os.path.isfile(filename):
            show_message("Cannot find default %s template file"%item)
            return {}
        parsed_name = item.capitalize()+'s'
        with open(filename) as openedfile:
            try:
                default_dict[parsed_name] = json.load(openedfile)
            except ValueError as error:
                show_message('JSON object issue: %s'%error)
                return {}
        if default_dict[parsed_name] is None:
            show_message("An error ocurred with the %s json template archive"%item)
            return {}
    return default_dict

def load_cylinder_template(nstroke, type_ig):
    type_ig_s = 'SI' if type_ig==0 else 'CI'
    filename = os.path.join(INSTALL_PATH,"templates","cylinder_default_%sstroke_%s.json"%(nstroke,type_ig_s))
    if not os.path.isfile(filename):
        show_message("Cannot find %s"%filename)
        return {}
    with open(filename) as openedfile:
        try:
            default_dict = json.load(openedfile)
        except ValueError as error:
            show_message('JSON object issue: %s'%error)
            return {}
    if default_dict is None:
        show_message("An error ocurred with the %s json template archive"%filename)
        return {}
    return default_dict

def explode_atribute(atribute):
    line = ''
    if not(isinstance(atribute,list)):
        try:
            if isinstance(atribute,str) or isinstance(atribute,unicode):
                line = line + "'" + atribute + "'"		
            else:
                line = line + str(atribute)
        except:
            if isinstance(atribute,str):
                line = line + "'" + atribute + "'"		
            else:
                line = line + str(atribute)
        return line
    else:	
        line = line + "["
        for l in range(len(atribute)):
            line = line + explode_atribute(atribute[l])			
            if l+1 < len(atribute):
                line = line + ", "
        line = line + "]"
    return line
    
def get_lines(current_dict, object_name, iobject):
    lines = object_name + str(iobject) + " = dict()\n"
    
    subdicts = []
    for ikey in current_dict:
        if type(current_dict[ikey])==dict:
            subdicts.append(ikey)
            continue
        lines += object_name + str(iobject) + "['" + ikey + "'] = " + explode_atribute(current_dict[ikey]) + "\n"
    
    for isubdict in subdicts:
        lines += init_object_lines(isubdict)
        lines += get_lines(current_dict[isubdict], isubdict, iobject)
        lines += object_name + str(iobject) + "['" + isubdict + "'] = " + isubdict + str(iobject) +"\n\n"
        lines += end_object_lines(isubdict)
        
    if object_name=='Valves' and current_dict['ncyl']!=-1:
        valve_type = "exhaust_valves" if current_dict['typeVal']=='exh' else "intake_valves"
        lines += "Cylinders" + str(current_dict['ncyl']) + "['" + valve_type + "'].append(" + object_name +  str(iobject) + ")\n\n"
    if object_name in ['Cylinders','Valves','Junctions','Tubes','Tanks','Atmospheres']:
        lines += "\n" + object_name + ".append(" + object_name +  str(iobject) + ")\n\n"
    elif object_name in ['Simulator']:
        lines += "\nSimulator = "	 + object_name +  str(iobject) + "\n\n"
    return lines
    
def init_object_lines(object_name):
    lines = "\n#--------- Inicializacion de " + object_name + "\n\n"
    if object_name in ['Cylinders','Valves','Junctions','Tubes','Tanks','Atmospheres']:
        lines += object_name + " = []\n\n"
    return lines

def end_object_lines(object_name):
    lines = "\n#--------- FIN Inicializacion de " + object_name + "\n\n"
    return lines

def save_data_aux(cw, objects, case_dir, case_name, filename = None, wizard = False):
    if not filename:
        filename = os.path.join(case_dir,case_name+".py")

    with open(filename, "w") as f:

        lines = ["#### ---- ####\n", "# Archivo generado por SimulatorGUI\n", "# CIMEC - Santa Fe - Argentina \n", 
                 "# Adecuado para levantar desde Interfaz Grafica \n", "# O para correr desde consola mediante $python main.py \n",
                 "#### ---- ####\n\n"]
        f.writelines(lines)

        # Simulator configuration
        current_dict = cw.save_data()
        lines = get_lines(current_dict, 'Simulator', 0)
        f.writelines(lines)

        # Para guardar, tengo que limpiar estos dos dicts. Luego los recupero porque
        # se usan para guardar los histogramas de los cilindros
        cyl_i = {}
        cyl_e = {}
        for icylinder in objects['Cylinders']:
            cyl_i[icylinder] = icylinder.object['intake_valves']
            cyl_e[icylinder] = icylinder.object['exhaust_valves']
            icylinder.object['intake_valves']   = []
            icylinder.object['exhaust_valves']  = []

        # Rest of the things
        object_names = ['Cylinders','Valves','Tubes','Tanks','Junctions','Atmospheres']
        for object_name in object_names:
            lines = init_object_lines(object_name)
            f.writelines(lines)
            for index, iobject in enumerate(objects[object_name]):
                current_dict = iobject.object
                current_dict['position'] = (int(iobject.pixmap.pos().x()),int(iobject.pixmap.pos().y()))
                lines = get_lines(current_dict, object_name, index)
                f.writelines(lines)
            lines = end_object_lines(object_name)
            f.writelines(lines)

        lines = "kargs = {'Simulator':Simulator, 'Cylinders':Cylinders, 'Junctions':Junctions, 'Tubes':Tubes, 'Tanks':Tanks, 'Atmospheres':Atmospheres}"
        f.writelines(lines)
    
    if not wizard:
        show_message('The file has been successfully saved', 1)

    for icylinder in objects['Cylinders']:
        icylinder.object['intake_valves']   = cyl_i[icylinder]
        icylinder.object['exhaust_valves']  = cyl_e[icylinder]
    return

def save_current_configuration_aux(parent, dict_to_save):
    """
    Save a defined user configuration
    """
    name = QtWidgets.QFileDialog.getSaveFileName(parent, 'Save Current Configuration', TEMPLATES_PATH, "json Files (*.json)")
    if name[0] != '':
        try:
            filename = name[0]
            if '.json' not in filename:
                filename = filename+'.json'
            with open(filename, 'w') as openedfile:
                json.dump(dict_to_save, openedfile)
            show_message('Configuration successfully saved!',1)
        except:
            show_message('There was an error trying to save this configuration.')
    return

def load_configuration_aux(parent, ctype):
    """
    Load a defined user configuration
    """
    dialog = QtGui.QFileDialog(parent)
    dialog.setNameFilter("json Files (*.json)")
    dialog.setWindowTitle('Load a %s configuration'%ctype)
    dialog.setDirectory(TEMPLATES_PATH)
    new_configuration = None
    if dialog.exec_():
        try:
            filename = dialog.selectedFiles()[0]
            with open(filename, "r") as openedfile:
                new_configuration = json.load(openedfile)
            return (True,new_configuration)
        except:
            show_message('There was an error trying to load this configuration.')
    return (False,new_configuration)


def check_two_equals(nlist):
    """
    Given a list, checks if 2 items are equal in it
    """
    n = 0
    for i in range(n,len(nlist)):
        for j in range(n+1,len(nlist)):
            if nlist[i]==nlist[j]:
                return True
        n+=1
        if n==len(nlist):
            return False
    return False

def dir_size(start, follow_links = True, start_depth = 0):
    """
    Given a directory, returns its size in bytes
    """
    from exception_handling import handle_exception
    # Get a list of all names of files and subdirectories in directory start
    try: 
        dir_list = os.listdir(start)
        from exception_handling import CURRENT_EXCEPTION
        assert(not CURRENT_EXCEPTION)
    except:
        handle_exception('Cannot list directory %s'%start)
        return

    total = 0
    for item in dir_list:
        # Get statistics on each item--file and subdirectory--of start
        path = os.path.join(start, item)
        try: 
            stats = os.stat(path)
            from exception_handling import CURRENT_EXCEPTION
            assert(not CURRENT_EXCEPTION)
        except:
            handle_exception('Cannot list directory %s'%path)
            return
        # The size in bytes is in the seventh item of the stats tuple, so:
        total += stats[6]
        # recursive descent if warranted
        if os.path.isdir(path) and (follow_links or not os.path.islink(path)):
            bytes = dir_size(path, follow_links, start_depth+1)
            total += bytes
    return total

class SelectionDialog(object):
    def setupUi(self, Dialog, items):
        Dialog.setObjectName("Dialog")
        Dialog.setFixedSize(200, 188)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.listWidget = QtWidgets.QListWidget(Dialog)
        self.listWidget.setObjectName("listWidget")
        for i in items:
            item = QtWidgets.QListWidgetItem()
            item.setCheckState(QtCore.Qt.Unchecked)
            item.setText(i)
            self.listWidget.addItem(item)
        self.verticalLayout.addWidget(self.listWidget)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)
        self.retranslateUi(Dialog)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle("Select boundaries")
        __sortingEnabled = self.listWidget.isSortingEnabled()
        self.listWidget.setSortingEnabled(False)
        self.listWidget.setSortingEnabled(__sortingEnabled)