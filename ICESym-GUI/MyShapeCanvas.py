import wx
import wx.lib.ogl as ogl

class MyShapeCanvas(ogl.ShapeCanvas):
    def __init__(self, home, parent = None, id = -1, pos = wx.DefaultPosition, size = wx.DefaultSize, style = wx.BORDER, name = "ShapeCanvas"):
        ogl.ShapeCanvas.__init__(self, parent, id, pos, size, style, name)
        self.home_reference = home
