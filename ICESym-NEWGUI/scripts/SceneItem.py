#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 20:56:39 2019

@author: etekken
"""

from PyQt5 import QtCore, QtWidgets
from PyQt5.QtSvg import QGraphicsSvgItem
from utils import ICON_PATHS

OFFSET_IMAGES = {}
OFFSET_IMAGES['Valves']         = QtCore.QPoint(0,0)
OFFSET_IMAGES['Cylinders']      = QtCore.QPoint(0,0)
OFFSET_IMAGES['Tubes']          = QtCore.QPoint(0,0)
OFFSET_IMAGES['Junctions']      = QtCore.QPoint(0,0)
OFFSET_IMAGES['Tanks']          = QtCore.QPoint(0,0)
OFFSET_IMAGES['Atmospheres']    = QtCore.QPoint(0,0)

class SceneItem():
    def __init__(self, itype, position, iobject):
        self.type = itype
        self.position = position + OFFSET_IMAGES[itype]
        self.offset = QtCore.QPoint(0,0)
        self.object = iobject
        self.current_celd = [0,0]
        self.set_pixmap()
        return
    
    def set_pixmap(self):
        if self.type=='Valves':
            if self.object['typeVal']=='int':
                self.pixmap = QGraphicsSvgItem(ICON_PATHS['Valve_int'])
            else:
                self.pixmap = QGraphicsSvgItem(ICON_PATHS['Valve_exh'])
        else:
            self.pixmap = QGraphicsSvgItem(ICON_PATHS[self.type])

        self.pixmap.setFlags(QtWidgets.QGraphicsItem.ItemIsSelectable)
        self.pixmap.setPos(self.position-self.offset)
        return