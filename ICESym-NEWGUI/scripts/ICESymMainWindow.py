#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 16:45:43 2018

@author: etekken
"""

import os, json, copy, sys
from PyQt5 import QtCore, QtGui, QtWidgets
from ICESymMainWindow_ui import Ui_ICESymMainWindow
from configurationWidget import configurationWidget
from postProcessWidget import postProcessWidget
from ValveDialog import ValveDialog, configure_default_valve
from AtmosphereDialog import AtmosphereDialog
from TubeDialog import TubeDialog, configure_default_tube
from CylinderDialog import CylinderDialog
from TankDialog import TankDialog
from JunctionDialog import JunctionDialog
from PyQt5.QtSvg import QGraphicsSvgItem
from LogTabWidget import LogTabWidget
from SceneItem import SceneItem
from UsageDialog import UsageDialog
from utils import show_message, load_templates, save_data_aux, ICON_PATHS,\
                  ICON_PATHS_NC, DEFAULT_DVP, INSTALL_PATH, CASES_PATH,\
                  INDEX_TAB_MODELING, INDEX_TAB_RUN, INDEX_TAB_POSTPROCESS,\
                  TAB_INFORMATION

TREE_POSITION = {} 
TREE_POSITION['Cylinders']      = 0
TREE_POSITION['Tubes']          = 1
TREE_POSITION['Junctions']      = 2
TREE_POSITION['Tanks']          = 3
TREE_POSITION['Valves']         = 4
TREE_POSITION['Atmospheres']    = 5

OFFSET_CONECTION = 0

CONNECTION_RULES = {}
CONNECTION_RULES['Atmospheres']  = ['Tubes']
CONNECTION_RULES['Tubes']        = ['Atmospheres','Junctions','Tanks','Valves']
CONNECTION_RULES['Junctions']    = ['Tubes']
CONNECTION_RULES['Cylinders']    = ['Valves']
CONNECTION_RULES['Tanks']        = ['Tubes']
CONNECTION_RULES['Valves']       = ['Cylinders','Tubes']

class ICESymMainWindow(QtWidgets.QMainWindow):
    def __init__(self, case_name = None, case_dir = None, case_type = None):
        QtWidgets.QMainWindow.__init__(self)
        self.ui = Ui_ICESymMainWindow()
        self.ui.setupUi(self)
        self.setGeometry(0, 0, 1350, 760)
        # ubica la interfaz en el centro de la pantalla
        self.centerOnScreen()

        self.ui.tab_configuration.setAutoFillBackground(True)
        self.ui.tab_modeling.setAutoFillBackground(True)

        self.scene = QtWidgets.QGraphicsScene(self.ui.canvas_widget)
        self.view = QtWidgets.QGraphicsView(self.scene)
        self.view.setSceneRect(0, 0, 3000, 1500)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.view)
        self.ui.canvas_widget.setLayout(layout)
        self.view.setMouseTracking(True)
        self.view.mousePressEvent = self.mouse_press
        self.view.mouseDoubleClickEvent = self.double_click
        self.view.mouseMoveEvent = self.mouse_movement
        self.view.mouseReleaseEvent = self.mouse_release
        self.resizeEventOld = self.view.resizeEvent
        self.view.resizeEvent = self.resizeView
        self.view.keyPressEvent = self.key_press_event

        self.click_positions = []
        self.sel_rect = None   
        self.rectangle_selection_active = None
        self.item_to_copy = None
        self.current_mouse_position = QtCore.QPointF(0.0,0.0)
        self.current_selected_item = None
        self.ui.canvas_widget.show()
        self.centroids = []
        self.grid_lines = []

        # array con todos los items presentes
        self.scene_items = []
        # array con tuplas de items con conexiones
        self.scene_connections = []

        # array of ints. Every int represents the number of objects in scene 
        # (in the same order as the tree)
        self.number_of_objects = [0, 0, 0, 0, 0, 0]
        # dictionary with objects
        self.objects = {}
        self.objects['Valves']           = []
        self.objects['Tubes']            = []
        self.objects['Atmospheres']      = []
        self.objects['Junctions']        = []
        self.objects['Cylinders']        = []
        self.objects['Tanks']            = []
        
        self.configure_default = {}
        self.configure_default['Valves'] = configure_default_valve
        self.configure_default['Tubes']  = configure_default_tube

        self.old_size_view = self.view.size()
        
        # templates para cada objeto
        self.default_dict = load_templates()
        
        # widgets de cada tab
        self.cw     = None
        self.ppw    = None
        self.ltw    = None
        
        self.current_tab_widget_index = 0
        
        self.case_name  = case_name
        self.case_dir   = case_dir
        if case_type in (1,3): # caso abierto o creado con wizard (previamente guardado)
            self.open_data()
        elif case_type==2: # caso nuevo en blanco
            self.current_configuration = self.default_dict['Configurations']
            self.current_configuration['folder_name']    = '%s_folder'%case_name
            self.current_configuration['filesave_state'] = '%s_state'%case_name
            self.current_configuration['filesave_spd']   = '%s_species'%case_name
            self.set_configuration_run_and_postProcess_widgets()
            show_message('New case successfully created!',1)
        return

    def set_configuration_run_and_postProcess_widgets(self):
        if self.cw:
            self.ui.tab_configuration.layout().removeWidget(self.cw)
        if self.ppw:
            self.ui.tab_postProcess.layout().removeWidget(self.ppw)
        if self.ltw:
            self.ui.tab_run.layout().removeWidget(self.ltw)
        self.cw = configurationWidget(self.current_configuration, self.case_name)
        self.ui.tab_configuration.layout().addWidget(self.cw)
        self.cw.setAutoFillBackground(True)            
        self.ppw = postProcessWidget(self.current_configuration,self.objects)
        self.ui.tab_postProcess.layout().addWidget(self.ppw)
        self.ppw.setAutoFillBackground(True)
        self.ltw = LogTabWidget(self.case_name, self.case_dir, self.current_configuration['folder_name'],\
                                self.current_configuration['rpms'], self.save_data, self.plot_defaultPostProcess_after_run)
        self.ui.tab_run.layout().addWidget(self.ltw)
        self.ltw.setAutoFillBackground(True)
        return

    def addTreeChild(self, name, number_in_tree):
        tree_subitem = QtWidgets.QTreeWidgetItem()
        tree_subitem.setText(0, name[0:-1]+' '+str(self.number_of_objects[number_in_tree]))
        pixmap = QtGui.QPixmap()
        pixmap.load(ICON_PATHS_NC[name])
        icon = QtGui.QIcon(pixmap)
        tree_subitem.setIcon(0,icon)
        tree_item = self.ui.componentsTreeWidget.topLevelItem(number_in_tree)
        tree_item.addChild(tree_subitem)
        self.number_of_objects[number_in_tree] += 1
        return

    def removeTreeChild(self, itype, name):
        sub_item = self.ui.componentsTreeWidget.findItems(name,QtCore.Qt.MatchExactly|QtCore.Qt.MatchRecursive)
        tree_item = self.ui.componentsTreeWidget.topLevelItem(TREE_POSITION[itype])
        index_of_sub_item = tree_item.indexOfChild(sub_item[0])
        tree_item.removeChild(sub_item[0])
        assert(len(sub_item)==1)
        del sub_item[0]
        self.number_of_objects[TREE_POSITION[itype]] -= 1
        # Renumera todos los que quedaron por detras
        for tree_index in range(index_of_sub_item,tree_item.childCount()):
            tree_item_n = tree_item.child(tree_index)
            object_type = tree_item_n.text(0).split()[0]
            tree_item_n.setText(0, "%s %s"%(object_type,tree_index))
        return

    def addSceneItem(self, itype, position = QtCore.QPoint(0,0), iobject = None, new = False, copying = False):
        if new or copying:
            iobject = copy.deepcopy(self.default_dict[itype] if new else iobject)
            if 'label' in iobject.keys():
                if new:
                    new_label = '%s_%s'%(iobject['label'],len(self.objects[itype]))
                    if itype in self.configure_default.keys():
                        self.configure_default[itype](iobject,self.current_configuration['nstroke'])
                elif copying:
                    new_label = '%s_copy_0'%iobject['label']
                    c = 1
                    while not self.check_unique_label(new_label):
                        new_label = new_label = '%s_copy_%s'%(iobject['label'],c)
                        c = c+1
                iobject['label'] = new_label
        item = SceneItem(itype, position, iobject)
        self.scene_items.append(item)
        self.scene.addItem(item.pixmap)
        self.addTreeChild(itype.capitalize(), TREE_POSITION[itype])
        self.objects[itype].append(item)
        return

    def resizeView(self, event):
        for iline in self.grid_lines:
            self.scene.removeItem(iline)
        self.drawGrid( self.view.sceneRect() )
        for conection in self.scene_connections:
            self.scene.removeItem(conection[3])
            self.scene.removeItem(conection[4])
        for conection in self.scene_connections:
            self.scene.addItem(conection[3])
            self.scene.addItem(conection[4])
        
        self.resizeEventOld(event)
        for item in self.scene_items:
            pos = self.view.mapToScene(item.position)
            item.pixmap.setPos(pos)
        return

    def drawGrid(self, windowsize):        
        X_SIZE = 66
        Y_SIZE = 66
        centroids_x = []
        centroids_y = []
        self.centroids = []
        self.grid_lines = []
        pencil = QtGui.QPen(QtCore.Qt.gray, 1)
        for xi in range(X_SIZE, int(windowsize.width()), X_SIZE):
            p1 = self.view.mapToScene( QtCore.QPoint( float(xi) , 0) )
            p2 = self.view.mapToScene( QtCore.QPoint( float(xi) , windowsize.height()) )
            self.grid_lines.append( self.scene.addLine( QtCore.QLineF(p1, p2), pencil ) )
        for xi in range(int(X_SIZE/2), int(windowsize.width()), X_SIZE):
            centroids_x.append(xi)
        for yi in range(Y_SIZE, int(windowsize.height()), Y_SIZE):
            p1 = self.view.mapToScene( QtCore.QPoint( 0, float(yi) ) )
            p2 = self.view.mapToScene( QtCore.QPoint( windowsize.width(), float(yi)) )
            self.grid_lines.append( self.scene.addLine( QtCore.QLineF(p1, p2), pencil ) )
        for yi in range(int(Y_SIZE/2), int(windowsize.height()), Y_SIZE):
            centroids_y.append(yi)

        pencil = self.getLinePen(QtCore.Qt.black)
        for index_x,icentroid_x in enumerate(centroids_x):
            for index_y,icentroid_y in enumerate(centroids_y):
                pc = self.view.mapToScene( QtCore.QPoint( float(icentroid_x) , float(icentroid_y)) )
                self.centroids.append( pc )
        return
    
    def getLinePen(self, color):
        pencil = QtGui.QPen(color, 1)
        pencil.setStyle(QtCore.Qt.SolidLine)
        pencil.setCosmetic(True)
        return pencil

    def drawConection(self, item1, item2, new_connection = False):
        pencil = self.getLinePen(QtCore.Qt.black)
        pos_item1 = item1.pixmap.pos()
        pos_item2 = item2.pixmap.pos()
        l = self.scene.addLine( QtCore.QLineF( pos_item1.x() + item1.pixmap.boundingRect().width() + OFFSET_CONECTION, \
                                 pos_item1.y() + item1.pixmap.boundingRect().height()/2.0, pos_item2.x() - OFFSET_CONECTION, \
                                 pos_item2.y() + item2.pixmap.boundingRect().height()/2.0 - OFFSET_CONECTION), pencil )
        c1 = self.scene.addEllipse( QtCore.QRectF(pos_item1.x() + item1.pixmap.boundingRect().width() - 2.0, \
                                                 pos_item1.y() + item1.pixmap.boundingRect().height()/2.0 - 2.0, 4.0, 4.0), \
                                                 pencil, QtGui.QBrush(QtCore.Qt.black))
        c2 = self.scene.addEllipse( QtCore.QRectF(pos_item2.x() - 2.0, pos_item2.y() + item1.pixmap.boundingRect().height()/2.0 - 2.0, 4.0, 4.0), \
                                                 pencil, QtGui.QBrush(QtCore.Qt.black))
        self.scene_connections.append( (item1,item2,l,c1,c2) )
        if new_connection:
            self.modify_conections_items(item1,item2,new_connection)
        return
    
    def modify_conections_items(self, item1, item2, new_connection=False):
        """
        Setea los argumentos necesarios para cada item segun la conexion
        item1: item a la izquierda de la conexion
        item2: item a la derecha de la conexion
        """
        item1_index = self.objects[item1.type].index(item1)
        item2_index = self.objects[item2.type].index(item2)
        
        if item1.type=='Tubes':
            # Caso especial de valvula, hay que poner el cilindro 
            # con el cual conecta
            if item2.type=='Valves':
                item1.object['nright'] = item2.object['ncyl']
                item1.object['tright'] = 'cylinder'
            else:
                item1.object['nright'] = item2_index
                item1.object['tright'] = item2.type[0:-1].lower()
        if item1.type=='Valves':
            keydict = 'ncyl' if item2.type=='Cylinders' else 'tube'
            item1.object[keydict] = item2_index
            if item2.type=='Cylinders' and new_connection:
                item2.object['state_ini'].append(DEFAULT_DVP)
                item2.object['nnod'] = item2.object['nnod']+1
        if item1.type=='Tanks':
            item1.object['exh2tube'].append( item2_index )
            if new_connection:
                item1.object['state_ini'].append ( item2.object['state_ini'][-1] )
                item1.object['Cd_ports'].append ( 0.8 )
                item1.object['nnod'] = item1.object['nnod']+1
        if item1.type=='Junctions':
            item1.object['type_end'].append(-1)
            item1.object['node2tube'].append( item2_index )
            if new_connection:
                item1.object['state_ini'].append(DEFAULT_DVP)
                item1.object['nnod'] = item1.object['nnod']+1
        if item1.type=='Cylinders':
            item1.object['exhaust_valves'].append(item2.object)
        
        if item2.type=='Tubes':
            # Caso especial de valvula, hay que poner el cilindro 
            # con el cual conecta
            if item1.type=='Valves':
                item2.object['nleft'] = item1.object['ncyl']
                item2.object['tleft'] = 'cylinder'
            else:
                item2.object['nleft'] = item1_index
                item2.object['tleft'] = item1.type[0:-1].lower()
        if item2.type=='Valves':
            keydict = 'ncyl' if item1.type=='Cylinders' else 'tube'
            item2.object[keydict] = item1_index
            if item1.type=='Cylinders' and new_connection:
                item1.object['state_ini'].append(DEFAULT_DVP)
                item1.object['nnod'] = item1.object['nnod']+1
        if item2.type=='Tanks':
            item2.object['int2tube'].append( item1_index )
            if new_connection:
                item2.object['state_ini'].append ( item1.object['state_ini'][-1] )
                item2.object['Cd_ports'].append ( 0.8 )
                item2.object['nnod'] = item2.object['nnod']+1
        if item2.type=='Junctions':
            item2.object['type_end'].append(1)
            item2.object['node2tube'].append( item1_index )
            if new_connection:
                item2.object['state_ini'].append(DEFAULT_DVP)
                item2.object['nnod'] = item2.object['nnod']+1
        if item2.type=='Cylinders':
            item2.object['intake_valves'].append(item1.object)
        return

    def deleteConnection(self, item1, item2):
        for connection in self.scene_connections:
            if  (item1 == connection[0] or item1 == connection[1]) and \
                (item2 == connection[0] or item2 == connection[1]):
                    self.scene_connections.remove(connection)
        return

    def existsConnection(self, item1, item2):
        eval_items = (item1,item2)
        for conection_items in self.scene_connections:
            cur_items = (conection_items[0],conection_items[1])
            if set(cur_items) == set(eval_items):
                show_message("There\'s a connection between this items")
                return True
        return False

    def existsRightConnection(self, item1):
        for conection_items in self.scene_connections:
            if item1 == conection_items[0] and item1.type in ('Valves','Cylinders'):
                show_message("There\'s a right connection for this item")
                return True
        return False

    def connectionIsAllowed(self,item1,item2):
        type1 = item1.type
        type2 = item2.type

        # Caso especial 1: Si quiero conectar un tubo con una valvula intake
        # que aun no esta conectada a un cilindro, debo primero asegurarme esta conexion
        if type1=='Tubes' and type2=='Valves':
            if item2.object['typeVal']=='int' and item2.object['ncyl']==-1:
                message = 'To connect a Tube with an Intake Valve, first connect the last one with a Cylinder'
                show_message(message)
                return False
            if item2.object['typeVal']=='exh':
                message = 'The Valve must be for Intake'
                show_message(message)
                return False
        # Caso especial 2: Si quiero conectar un cilindro con una valula, esta
        # debe ser exhaust
        if type1=='Cylinders' and type2=='Valves':
            if item2.object['typeVal']=='int':
                message = 'The Valve must be for Exhaust'
                show_message(message)
                return False
        # Caso especial 3: Si quiero conectar una valvula con un cilindro, esta
        # debe ser intake
        if type1=='Valves' and type2=='Cylinders':
            if item1.object['typeVal']=='exh':
                message = 'The Valve must be for Intake'
                show_message(message)
                return False

        if type2 in CONNECTION_RULES[type1]:
            return True
        else:
            message = type1 + " can connect only with: \n\n"
            for itype in CONNECTION_RULES[type1]:
                message += "- " + itype + "\n"
            show_message(message)
            return False

    def updateCaseSetup(self, QTreeWidgetItem, index):
        if not QTreeWidgetItem:
            return
        component = QTreeWidgetItem.text(0)

        if 'Tube' in component:
            item_index = component.replace('Tube ','')
            try:
                item_index = int(item_index)
                item = self.objects['Tubes'][item_index]
                dialog = TubeDialog(item.object,item_index)
            except:
                return
        elif 'Cylinder' in component:
            item_index = component.replace('Cylinder ','')
            try:
                item_index = int(item_index)
                item = self.objects['Cylinders'][item_index]
                dialog = CylinderDialog(item.object,item_index)
            except:
                return
        elif 'Junction' in component:
            item_index = component.replace('Junction ','')
            try:
                item_index = int(item_index)
                item = self.objects['Junctions'][item_index]
                dialog = JunctionDialog(item.object,item_index)
            except:
                return
        elif 'Tank' in component:
            item_index = component.replace('Tank ','')
            try:
                item_index = int(item_index)
                item = self.objects['Tanks'][item_index]
                dialog = TankDialog(item.object,item_index)
            except:
                return
        elif 'Atmosphere' in component:
            item_index = component.replace('Atmosphere ','')
            try:
                item_index = int(item_index)
                item = self.objects['Atmospheres'][item_index]
                dialog = AtmosphereDialog(item.object,item_index)
            except:
                return
        elif 'Valve' in component:
            item_index = component.replace('Valve ','')
            try:
                item_index = int(item_index)
                item = self.objects['Valves'][item_index]
                dialog = ValveDialog(item,item_index,self)
            except:
                return
        else:
            return

        if dialog:
            dialog.exec_()
            item.object = dialog.current_dict
        return
    
    def load_templates(self):
        """
        load the default templates
        """
        items = ['valve','atmosphere','cylinder','junction','tube','tank','configuration']
        for item in items:
            filename = os.path.join(INSTALL_PATH,"templates","%s_default.json"%item)
            if not os.path.isfile(filename):
                show_message("Cannot find default %s template file"%item)
                return False
            parsed_name = item.capitalize()+'s'
            with open(filename) as openedfile:
                try:
                    self.default_dict[parsed_name] = json.load(openedfile)
                except ValueError as error:
                    show_message('JSON object issue: %s'%error)
                    return
            if self.default_dict[parsed_name] is None:
                show_message("An error ocurred with the %s json template archive"%item)
                return False
        return True
    
    def centerOnScreen(self):
        '''Centers the window on the screen.'''
        resolution = QtWidgets.QDesktopWidget().screenGeometry()
        self.move((resolution.width() / 2) - (self.frameSize().width() / 2),
                  (resolution.height() / 2) - (self.frameSize().height() / 2))
        
    def searchItem(self, itype, tobject):
        for sitem in self.objects[itype]:
            if sitem.object == tobject:
                return sitem
        return None
    
    def drawAllConnections(self):
        for itype in self.objects.keys():
            for sitem in self.objects[itype]:
                itype = sitem.type
                if itype == "Tubes":
                    sitem2type_r  = sitem.object['tright']
                    sitem2index_r = sitem.object['nright']
                    sitem2type_l  = sitem.object['tleft']
                    sitem2index_l = sitem.object['nleft']
                    if sitem2index_r != -1:
                        if sitem2type_r in ['tank','tube','junction','atmosphere']:
                            sitem2type_r = sitem2type_r.title() + 's'
                            sitem2 = self.objects[sitem2type_r][sitem2index_r]
                        else:
                            if sitem2type_r == 'cylinder': # si esta conectado a un cilindro va a una valvula de por medio
                                sitem2 = self.objects['Valves'][sitem2index_r]
                        self.drawConection(sitem, sitem2, False)
                        if sitem2index_l != -1 and sitem2type_l == 'atmosphere': # a la izquierda una atmosphere
                            sitem2 = self.objects['Atmospheres'][sitem2index_l]
                            self.drawConection(sitem2, sitem, False) # notar el invertido de items

                if itype == "Junctions":
                    lenTubes = len(sitem.object['node2tube'])
                    for i in range(0,lenTubes):
                        if sitem.object['type_end'][i] == -1: # de escape, es decir a izquierda
                            sitem2 = self.objects['Tubes'][ sitem.object['node2tube'][i] ]
                            self.drawConection(sitem, sitem2, False)

                if itype == "Tanks":
                    lenTubes = len(sitem.object['exh2tube'])
                    for i in range(0,lenTubes):
                        sitem2 = self.objects['Tubes'][ sitem.object['exh2tube'][i] ]
                        self.drawConection(sitem, sitem2, False)

                if itype == "Cylinders":
                    lenValvess = len(sitem.object['exhaust_valves'])
                    for i in range(0,lenValvess):
                        Valves = sitem.object['exhaust_valves'][i]
#                        sitem.object['intake_valves'] = []      # este key es solo para dibujar, luego jode
#                        sitem.object['exhaust_valves']  = []    # este key es solo para dibujar, luego jode
                        sitem2 = self.searchItem('Valves', Valves)
                        if sitem2 != None:
                            self.drawConection(sitem, sitem2, False)

                if itype == "Valves":
                    if sitem.object['typeVal'] == "exh" and sitem.object['tube']!=-1:
                        sitem2 = self.objects['Tubes'][ sitem.object['tube'] ]
                        self.drawConection(sitem, sitem2, False)
                    if sitem.object['typeVal'] == "int" and sitem.object['ncyl']!=-1:
                        sitem2 = self.objects['Cylinders'][ sitem.object['ncyl'] ]
                        self.drawConection(sitem, sitem2, False)
        return
    
    def deleteTanksAndJunctionsTubes(self):
        for itanks in self.objects['Tanks']:
            itanks.object['int2tube'] = []
            itanks.object['exh2tube'] = []
        for ijuncs in self.objects['Junctions']:
            ijuncs.object['node2tube'] = []
            ijuncs.object['type_end'] = []
        return
    
    def orderTanksConnections(self):
        for itanks in self.objects['Tanks']:
            itanks.object['int2tube'].sort()
            itanks.object['exh2tube'].sort()
        return
    
    def eraseItem(self, item):
        """
        Dado un item específico, lo elimina a el mismo y todas sus conexiones.
        Debe tambien cambiar la numeracion de las conexiones
        """
        self.scene.removeItem(item.pixmap)
        self.scene_items.remove(item)

        to_delete=[]
        for connection in self.scene_connections:
            item1 = connection[0]
            item2 = connection[1]
            if item==item1 or item==item2:
                to_delete.append(connection)
        for connection in to_delete:
            self.eraseConnection(connection)
            
        if item.type=='Cylinders':
            self.cw.edit_ig_order(self.objects[item.type].index(item),True)

        self.objects[item.type].remove(item)

        self.deleteTanksAndJunctionsTubes()

        first_connections   = [c for c in self.scene_connections if c[0].type==item.type\
                             or c[1].type==item.type]
        rest_of_connections = [c for c in self.scene_connections if c[0].type!=item.type\
                             and c[1].type!=item.type]
        # Primero modifico las conexiones que involucren al item, 
        # luego el resto..
        for connection in first_connections:
            self.modify_conections_items(connection[0],connection[1])
        for connection in rest_of_connections:
            self.modify_conections_items(connection[0],connection[1])

        self.orderTanksConnections()
        return
    
    def eraseConnection(self, connection):
        """
        Elimina una conexion dada
        """
        if connection in self.scene_connections:
            self.scene.removeItem(connection[2])
            self.scene.removeItem(connection[3])
            self.scene.removeItem(connection[4])
            self.scene_connections.remove(connection)
        item1 = connection[0]
        item2 = connection[1]
        if item1.type=='Junctions' or item2.type=='Junctions':
            item_junction = item1 if item1.type=='Junctions' else item2
            item_tube = item1 if item1.type=='Tubes' else item2
            tube_number = self.objects['Tubes'].index(item_tube)
            tube_index = item_junction.object['node2tube'].index(tube_number)
            item_junction.object['node2tube'].remove(tube_number)
            item_junction.object['type_end'].pop(tube_index)
            item_junction.object['state_ini'].pop(-1)
            item_junction.object['nnod'] = item_junction.object['nnod']-1
            item_junction.object['histo'] = []
        if item1.type=='Tanks' or item2.type=='Tanks':
            item_tank = item1 if item1.type=='Tanks' else item2
            item_tube = item1 if item1.type=='Tubes' else item2
            tube_number = self.objects['Tubes'].index(item_tube)
            tube_type = 'int2tube' if item_tube==item1 else 'exh2tube'
            item_tank.object[tube_type].remove(tube_number)
            item_tank.object['state_ini'].pop(-1)
            item_tank.object['Cd_ports'].pop(-1)
            item_tank.object['nnod'] = item_tank.object['nnod']-1
            item_tank.object['histo'] = []
        if item1.type=='Tubes' or item2.type=='Tubes':
            item_tube = item1 if item1.type=='Tubes' else item2
            item_tube.object['nleft' if item2==item_tube else 'nright'] = -1
            item_tube.object['tleft' if item2==item_tube else 'tright'] = "<none>"
        if item1.type=='Valves' or item2.type=='Valves':
            item_valve = item1 if item1.type=='Valves' else item2
            item_valve.object['ncyl'] = -1
            item_valve.object['tube'] = -1
            # eliminar ambas conexiones, con el tubo y con el cilindro
            for connection_v in self.scene_connections:
                item1_c = connection_v[0]
                item2_c = connection_v[1]
                if item_valve == item1_c or item_valve == item2_c:
                    self.eraseConnection(connection_v)
        if item1.type=='Cylinders' or item2.type=='Cylinders':
            item_cylinder = item1 if item1.type=='Cylinders' else item2
            item_cylinder.object['nnod'] = item_cylinder.object['nnod']-1
            item_cylinder.object['state_ini'].pop(-1)
            item_valve = item1 if item1.type=='Valves' else item2
            type_valve = 'intake_valves' if item2==item_cylinder else 'exhaust_valves'
            item_cylinder.object[type_valve].remove(item_valve.object)
        return

    def drawObjects(self, itype, dict_data):
        for idata in dict_data:
            self.addSceneItem(itype, QtCore.QPoint(idata['position'][0], idata['position'][1]), idata, False, False)
        return
    
    def cleanObjects(self):
        for ikey in self.objects.keys():
            self.objects[ikey] = []
        self.scene.clear()
        self.scene_items = []
        self.scene_connections = []
        self.number_of_objects = [0, 0, 0, 0, 0, 0]        
        for itype in TREE_POSITION.keys():
            tree_item = self.ui.componentsTreeWidget.topLevelItem(TREE_POSITION[itype])
            while tree_item.childCount():
                it = tree_item.child(0)
                tree_item.removeChild(it)
                del it
        return

    def check_valve_type(self, scene_item):
        self.scene.removeItem(scene_item.pixmap)
        if scene_item.object['typeVal']=='int':
            scene_item.pixmap = QGraphicsSvgItem(ICON_PATHS['Valve_int'])
        else:
            scene_item.pixmap = QGraphicsSvgItem(ICON_PATHS['Valve_exh'])
        scene_item.pixmap.setFlags(QtWidgets.QGraphicsItem.ItemIsSelectable)
        scene_item.pixmap.setPos(scene_item.position)        
        self.scene.addItem(scene_item.pixmap)
        return

    def find_optimum_cell_position(self, current_position):
        # Buscar que centroide se acerca mas a la posicion
        min_dist = 100
        import math
        min_centroid = self.centroids[0]
        for icentroid in self.centroids:
            dist =  math.sqrt((icentroid.x()-current_position.x())**2+(icentroid.y()-current_position.y())**2)
            if dist<min_dist:
                min_dist = dist
                min_centroid = icentroid
        return min_centroid+QtCore.QPointF(-32.0,-32.0)
    
    def check_unique_label(self, object_label, show_message_b = False):
        for item in self.scene_items:
            if 'label' in item.object.keys():
                if item.object['label'] == object_label:
                    if show_message_b:
                        show_message('There is already an object with this name. Please, select another one.')
                    return False
        return True

# Keyboard manipulation and auxiliary functions
# -----------------------------------------------------------------------------

    def key_press_event(self, event):
        if not self.check_tab(INDEX_TAB_MODELING):
            return
        if event.key() == QtCore.Qt.Key_Delete:
            for item in self.scene_items:
                if item.pixmap.isSelected():
                    item_name = item.type[0:-1] + " " + str(self.objects[item.type].index(item))
                    msg = "You will erase %s and all their connections. Do you want to continue?"%item_name
                    reply = show_message(msg,4,QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No)
                    if reply == QtWidgets.QMessageBox.Yes:
                        try:
                            self.eraseItem(item)
                            self.removeTreeChild(item.type,item_name)
                            show_message("Item erased", 1)
                        except:
                            show_message("There was an error trying to erase this item")
                    else:
                        show_message("Operation cancelled", 1)
                        return

            for connection in self.scene_connections:
                item1 = connection[0]
                item2 = connection[1]
                line = connection[2]
                if line.pen().color()==QtCore.Qt.red:
                    msg = "You will erase the connection between %s %s and %s %s. Do you want to continue?"\
                    %(item1.type[0:-1],self.objects[item1.type].index(item1),item2.type[0:-1],self.objects[item2.type].index(item2))
                    reply = show_message(msg,4,QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No)
                    if reply == QtWidgets.QMessageBox.Yes:
                        try:
                            self.eraseConnection(connection)
                            show_message("Connection erased", 1)
                        except:
                            show_message("There was an error trying to erase this connection")
                    else:
                        show_message("Operation cancelled", 1)
                        return
                    return
        if event.modifiers() == QtCore.Qt.ControlModifier:
            if event.key() == QtCore.Qt.Key_C:
                for item in self.scene_items:
                    if item.pixmap.isSelected():
                        self.item_to_copy = item
            if event.key() == QtCore.Qt.Key_V:
                if self.item_to_copy:
                    msg = "Do you want to paste the selected %s?"%self.item_to_copy.type[0:-1]
                    reply = show_message(msg,4,QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No)
                    if reply == QtWidgets.QMessageBox.Yes:
                        object_position = self.find_optimum_cell_position(self.current_mouse_position)
                        self.addSceneItem(self.item_to_copy.type, object_position,\
                                          self.item_to_copy.object, False, True)
                        self.item_to_copy = None
        return

# Mouse manipulation and auxiliary functions
# -----------------------------------------------------------------------------

    def double_click(self, event):
        if not self.check_tab(INDEX_TAB_MODELING):
            return
        for item in self.scene_items:
            if item.pixmap.isSelected():
                dialog = None
                item_index = self.objects[item.type].index(item)
                if (item.type == 'Valves'):
                    # A valve le debo pasar el SceneItem y el parent porque en
                    # caso de que cambie el tipo, debe cambiar el icono
                    dialog = ValveDialog(item,item_index,self)
                if (item.type == 'Atmospheres'):
                    dialog = AtmosphereDialog(item.object,item_index)
                if (item.type == 'Tubes'):
                    dialog = TubeDialog(item.object,item_index)
                if (item.type == 'Cylinders'):
                    dialog = CylinderDialog(item.object,item_index)
                if (item.type == 'Tanks'):
                    dialog = TankDialog(item.object,item_index)
                if (item.type == 'Junctions'):
                    dialog = JunctionDialog(item.object,item_index)
                if dialog:
                    dialog.exec_()
                    item.object = dialog.current_dict
        self.current_selected_item = None
        return

    def mouse_movement(self, event):
        if not self.check_tab(INDEX_TAB_MODELING):
            return

        # TODO: under construction
        if self.rectangle_selection_active:
            self.rect_sel_end = event.pos()
            pencil = QtGui.QPen(QtCore.Qt.black, 1)
            pencil.setStyle(QtCore.Qt.DashLine)
            pencil.setCosmetic(True)
            if self.sel_rect:
                self.scene.removeItem(self.sel_rect)
            self.sel_rect = self.scene.addRect( QtCore.QRectF(self.rect_sel_beg, self.rect_sel_end), \
                                               pencil, QtGui.QBrush(QtCore.Qt.NoBrush) )            

        # mapFromScene mapea las coordenadas del objeto, que tienen como origen el vertice
        # inferior izquierdo y el centro de coordenadas es el centro de la pantalla, en coordenadas
        # globales donde el 0,0 del objeto esta en el vertice superior izquierdo y el de la vista
        # en el vertice superior izquierdo
        # event.pos() da la posicion del click, con el origen en el vertice inferior izq
        if self.current_selected_item and self.current_selected_item.pixmap.isSelected():
            if self.current_selected_item.offset == -1:
                self.current_selected_item.offset = event.pos() - \
                self.view.mapFromScene(QtCore.QPoint(self.current_selected_item.pixmap.pos().x(),\
                                                     self.current_selected_item.pixmap.pos().y()))
            current_position = self.view.mapToScene( (event.pos() - \
                                              self.current_selected_item.offset).x(),\
                                             (event.pos() - self.current_selected_item.offset).y() )
            
            position = self.find_optimum_cell_position(current_position)
            self.current_selected_item.pixmap.setPos(position)            
            offset = event.pos() - QtCore.QPoint(self.current_selected_item.pixmap.pos().x(), \
                               self.current_selected_item.pixmap.pos().y())
            position = QtCore.QPoint((event.pos() - offset).x(), (event.pos() - offset).y())
            self.current_selected_item.position = position
            for conection in self.scene_connections:
                if self.current_selected_item == conection[0] or \
                   self.current_selected_item == conection[1]:
                    # ver como hacer mas optimo esto
                    self.deleteConnection(conection[0], conection[1])
                    self.scene.removeItem(conection[2])
                    self.scene.removeItem(conection[3])
                    self.scene.removeItem(conection[4])
                    self.drawConection(conection[0], conection[1])

        self.current_mouse_position = self.view.mapToScene(event.pos())
        return

    def mouse_release(self, event):
        if not self.check_tab(INDEX_TAB_MODELING):
            return
        if self.rectangle_selection_active:
            self.rectangle_selection_active = False
            self.scene.removeItem(self.sel_rect)
            
        self.current_selected_item = None

        # Ver si hay alguna conexion entre dos items
        if event.button() == QtCore.Qt.RightButton:
            for item1 in self.scene_items:
                if item1.pixmap.isSelected():
                    offset = self.view.mapFromScene( QtCore.QPoint(0,0) )
                    posMouse = event.pos() - offset
                    for item2 in self.scene_items:
                        if item1 == item2:
                            continue
                        item2x = item2.pixmap.x()
                        item2y = item2.pixmap.y()
                        item2w = item2.pixmap.boundingRect().width()
                        item2h = item2.pixmap.boundingRect().height()
                        # Si existe un objeto seleccionado, y el release se hace 
                        # sobre otro objeto, se dibuja una conexion entre
                        # ellos, apuntando al segundo (solo si no existe uno)
                        if self.check_if_item_is_inside(posMouse, item2x, item2y, item2w, item2h)\
                        and not self.existsRightConnection(item1)\
                        and not self.existsConnection(item1,item2)\
                        and self.connectionIsAllowed(item1,item2):
                            self.drawConection(item1, item2, True)
        return

    def check_if_item_is_inside(self, click, posx, posy, width, height):
        # Ver si un click fue adentro de un item
        return ( click.x() > posx and click.x() < posx + width and click.y() > posy and click.y() < posy + width )

    def check_if_line_is_inside(self, click, p1, p2, length):
        # Ver si un click fue adentro de una linea (uso ffs)
        dxc = click.x() - p1.x()
        dyc = click.y() - p1.y()
        dxl = p2.x() - p1.x()
        dyl = p2.y() - p1.y()
        cross = dxc*dyl - dyc*dxl
        if abs(cross)>400.0:
            return False
        alpha1 = QtCore.QPointF.dotProduct(click-p1,p2-p1)/(length**2)
        alpha2 = QtCore.QPointF.dotProduct(p2-click,p2-p1)/(length**2)
        return alpha1>0.0 and alpha1<1.0 and alpha2>0.0 and alpha2<1.0

    def mouse_press(self, event):
        if not self.check_tab(INDEX_TAB_MODELING):
            return
        offset = self.view.mapFromScene(QtCore.QPoint(0,0))
        mouseClick = event.pos() - offset

        if event.button() == QtCore.Qt.LeftButton: # Check if there is a selected item
            # TODO: under construction
            if event.modifiers() == QtCore.Qt.ControlModifier:
                self.view.setDragMode(QtWidgets.QGraphicsView.RubberBandDrag)
                self.rect_sel_beg = event.pos()
                self.rect_sel_end = event.pos()
                self.rectangle_selection_active = True

            self.current_selected_item = None
            for item in self.scene_items:
                item.pixmap.setSelected(False)
                item.offset = -1
            for connection in self.scene_connections: # Check if a connection is selected
                line = connection[2]
                pencil = self.getLinePen(QtCore.Qt.black)
                line.setPen(pencil)
                line.setSelected(False)

            for item in self.scene_items: # Check if an item is selected
                if self.check_if_item_is_inside(mouseClick, item.pixmap.x(), item.pixmap.y(),\
                                                item.pixmap.boundingRect().width(),\
                                                item.pixmap.boundingRect().height()):
                    item.pixmap.setSelected(True)
                    self.current_selected_item = item
                    return

            for connection in self.scene_connections: # Check if a connection is selected
                line = connection[2]
                p1 = line.line().p1()
                p2 = line.line().p2()
                length = line.line().length()
                if self.check_if_line_is_inside(mouseClick, p1, p2, length):
                    pencil = self.getLinePen(QtCore.Qt.red)
                    line.setPen(pencil)
                    # line.setSelected(True) # Not working
                    return
        return

# Free Actions
# -----------------------------------------------------------------------------
    
    def open_terminal(self):
        command = "xterm -e 'cd %s && /bin/bash' &"%INSTALL_PATH
        os.system(command)
        return

    def close_case(self):
        msg = "Do you want to close the current case?"
        reply = show_message(msg,4,QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No)
        if reply == QtWidgets.QMessageBox.Yes:
            self.cleanObjects()
            self.drawGrid( self.view.sceneRect() )
            self.current_configuration = self.default_dict['Configurations']
            self.case_name = 'default_case'
            self.case_dir = CASES_PATH
            self.set_configuration_run_and_postProcess_widgets()
            show_message('Case successfully closed!', 1)
        return
    
    def open_data(self):
        sys.path.append(str(self.case_dir))
        externalData = __import__(str(self.case_name))
        self.cleanObjects()
        self.drawGrid( self.view.sceneRect() )

        self.current_configuration = externalData.Simulator
        self.drawObjects('Cylinders',   externalData.Cylinders)
        self.drawObjects('Tanks',       externalData.Tanks)
        self.drawObjects('Junctions',   externalData.Junctions)
        self.drawObjects('Tubes',       externalData.Tubes)
        self.drawObjects('Valves',      externalData.Valves)
        self.drawObjects('Atmospheres', externalData.Atmospheres)
        
        self.drawAllConnections()
        self.set_configuration_run_and_postProcess_widgets()
        
        show_message('Case successfully opened!', 1)
        return

    def open_case(self):      
        dialog = QtGui.QFileDialog(self)
        dialog.setNameFilter("Python Files (*.py)")
        dialog.setWindowTitle('Open an ICESym-GUI Case')
        dialog.setDirectory(INSTALL_PATH)
        if dialog.exec_():
            filename = dialog.selectedFiles()[0]
            with open(filename, "r") as f:
                line = f.readline()
                if line=="#### ---- ####\n":                
                    pathName =  os.path.dirname(filename)
                    moduleName =  os.path.basename(filename).replace('.py','')
                    self.case_name = moduleName
                    self.case_dir = pathName
                    self.open_data()
                else:
                    show_message('Please select a valid ICESym case to open')
        else:
            show_message('Nothing to open', 1)
        return
    
    def save_as(self):
        name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File As', "./", "Python Files (*.py)")
        if name[0] != '':
            filename = name[ 0 ]
            filename = filename+'.py' if filename.find('.py')==-1 else filename
            self.save_data(filename)
        return
    
    def save_data(self, filename = None, wizard = None):
        try:
            save_data_aux(self.cw, self.objects, self.case_dir, self.case_name, filename, wizard)
            self.ltw.change_attributes(self.current_configuration['folder_name'],self.current_configuration['rpms'])
            self.ppw.change_attributes(self.current_configuration,self.objects)
        except:
            show_message('An error has occurred trying to save the file')
        return
    
    def exit(self):
        msg = "Do you want to exit ICESym?"
        reply = show_message(msg,4,QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No)
        if reply == QtWidgets.QMessageBox.Yes:
            self.close()
        return

    def show_usage(self):
        ud = UsageDialog()
        ud.exec_()
        return
    
    def plot_defaultPostProcess_after_run(self):
        self.ui.tabWidget.setCurrentIndex(INDEX_TAB_POSTPROCESS)
        self.ui.actionDefault_Post_Process.setEnabled(False)
        self.ppw.plot_defaults()
        self.ui.actionDefault_Post_Process.setEnabled(True)
        return

# Tabs check/manipulation auxiliary functions
# -----------------------------------------------------------------------------

    def check_tab(self, tab_index, show_error_message=True):
        if self.current_tab_widget_index!=tab_index:
                if show_error_message:
                    show_message('This action is enabled only in %s Tab'%TAB_INFORMATION[tab_index]['name'])
                return False
        return True

    def change_tab(self, tab_index):
        if tab_index==INDEX_TAB_POSTPROCESS:
            self.ppw.enable_ppw()
        self.current_tab_widget_index = tab_index
        return

# Modeling Tab Actions
# -----------------------------------------------------------------------------
        
    def addValve(self):
        if not self.check_tab(INDEX_TAB_MODELING):
            return
        self.addSceneItem('Valves',QtCore.QPoint(0,0),None,True,False)
        return

    def addAtmosphere(self):
        if not self.check_tab(INDEX_TAB_MODELING):
            return
        self.addSceneItem('Atmospheres',QtCore.QPoint(0,0),None,True,False)
        return

    def addTube(self):
        if not self.check_tab(INDEX_TAB_MODELING):
            return
        self.addSceneItem('Tubes',QtCore.QPoint(0,0),None,True,False)
        return

    def addTank(self):
        if not self.check_tab(INDEX_TAB_MODELING):
            return
        self.addSceneItem('Tanks',QtCore.QPoint(0,0),None,True,False)
        return

    def addJunction(self):
        if not self.check_tab(INDEX_TAB_MODELING):
            return
        self.addSceneItem('Junctions',QtCore.QPoint(0,0),None,True,False)
        return

    def addCylinder(self):
        if not self.check_tab(INDEX_TAB_MODELING):
            return
        self.addSceneItem('Cylinders',QtCore.QPoint(0,0),None,True,False)
        self.cw.edit_ig_order(len(self.objects['Cylinders'])-1,False)
        return

# Run Tab Actions
# -----------------------------------------------------------------------------

    def run_simulation(self):
        if not self.check_tab(INDEX_TAB_RUN):
            return
        self.ltw.run_simulation()
        return

    def kill_simulation(self):
        if not self.check_tab(INDEX_TAB_RUN):
            return
        self.ltw.kill_process()
        return

# Post Process Tab Actions
# -----------------------------------------------------------------------------
        
    def save_postpro(self):
        if not self.check_tab(INDEX_TAB_POSTPROCESS):
            return
        self.ppw.save_postpro()
        return

    def load_postpro(self):
        if not self.check_tab(INDEX_TAB_POSTPROCESS):
            return
        self.ppw.load_postpro()
        return

    def plot_defaults(self):
        if not self.check_tab(INDEX_TAB_POSTPROCESS):
            return
        self.ui.actionDefault_Post_Process.setEnabled(False)
        self.ppw.plot_defaults()
        self.ui.actionDefault_Post_Process.setEnabled(True)
        return