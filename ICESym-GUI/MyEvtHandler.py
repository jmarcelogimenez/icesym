#clase que controla los eventos
#en prueba aun

import wx
import wx.lib.ogl as ogl

class MyEvtHandler(ogl.ShapeEvtHandler):
    def __init__(self, log, frame):
        ogl.ShapeEvtHandler.__init__(self)
        self.log = log
        self.statbarFrame = frame

    def UpdateStatusBar(self, shape):
        x, y = shape.GetX(), shape.GetY()
        width, height = shape.GetBoundingBoxMax()
        #self.statbarFrame.SetStatusText("Pos: (%d, %d)  Size: (%d, %d)" %
                                       # (x, y, width, height))


    def OnLeftClick(self, x, y, keys=0, attachment=0):
        print "in event OnLeftClick"
        shape = self.GetShape()
        canvas = shape.GetCanvas()
        dc = wx.ClientDC(canvas)
        canvas.PrepareDC(dc)

        if shape.Selected():
            shape.Select(False, dc)
            #canvas.Redraw(dc)
            canvas.Refresh(False)
        else:
            redraw = False
            shapeList = canvas.GetDiagram().GetShapeList()
            toUnselect = []
            
            index = 0
            for s in shapeList:
                if s.Selected():
                    # If we unselect it now then some of the objects in
                    # shapeList will become invalid (the control points are
                    # shapes too!) and bad things will happen...
                    toUnselect.append(s)
             
            shape.Select(True, dc)

            if toUnselect:
                for s in toUnselect:
                    s.Select(False, dc)

                ##canvas.Redraw(dc)
                canvas.Refresh(False)

        self.UpdateStatusBar(shape)

        
    def OnEndDragLeft(self, x, y, keys=0, attachment=0):
        shape = self.GetShape()
        ogl.ShapeEvtHandler.OnEndDragLeft(self, x, y, keys, attachment)
        canvas = shape.GetCanvas()
        home = canvas.home_reference
        shapeList = canvas.GetDiagram().GetShapeList()

        if not shape.Selected():
            self.OnLeftClick(x, y, keys, attachment)
        
        print "termina en x,y: ",x,y
        
        index = 0
        for s in shapeList:
            if isinstance(s,ogl.BitmapShape):
                if s == shape:
                    print index
                    break
                index = index + 1

        home.updatePosition(index,x,y)
        self.UpdateStatusBar(shape)


    def OnSizingEndDragLeft(self, pt, x, y, keys, attch):
        ogl.ShapeEvtHandler.OnSizingEndDragLeft(self, pt, x, y, keys, attch)
        self.UpdateStatusBar(self.GetShape())


    def OnMovePost(self, dc, x, y, oldX, oldY, display):
        shape = self.GetShape()
        ogl.ShapeEvtHandler.OnMovePost(self, dc, x, y, oldX, oldY, display)
        self.UpdateStatusBar(shape)
        if "wxMac" in wx.PlatformInfo:
            shape.GetCanvas().Refresh(False)

    def OnRightClick(self,  x, y, keys=0, attachment=0):
        shapeTo = self.GetShape()
        canvas = shapeTo.GetCanvas()
        home = canvas.home_reference
        shapeList = canvas.GetDiagram().GetShapeList()

        ToIndex = 0;
        FromIndex = 0;
        it = 0
        shapeFrom = shapeTo

        shapes = home.shapes
         
        if not(shapeTo.Selected()):
            for s in shapeList:
                if isinstance(s,ogl.BitmapShape):
                    if s == shapeTo:
                        ToIndex = it
                    if s.Selected():
                        FromIndex = it
                        shapeFrom = s
                    it = it+1
         
        if not(shapeTo == shapeFrom):
            home.makeConection(shapeFrom, FromIndex, shapeTo, ToIndex)

        canvas.Refresh(True)
        self.UpdateStatusBar(shapeFrom)
        

    def OnLeftDoubleClick(self,x, y, keys, attachment):
        print "in event double click"
        shape = self.GetShape()
        canvas = shape.GetCanvas()
        dc = wx.ClientDC(canvas)
        canvas.PrepareDC(dc)
        home = canvas.home_reference
        
        shapeList = canvas.GetDiagram().GetShapeList()

        if shape.Selected():
            shape.Select(False, dc)
            #canvas.Redraw(dc)
            canvas.Refresh(False)
        else:
            redraw = False
            shapeList = canvas.GetDiagram().GetShapeList()
            toUnselect = []

        index = 0
        for s in shapeList:
            if s.Selected():
                # If we unselect it now then some of the objects in
                # shapeList will become invalid (the control points are
                # shapes too!) and bad things will happen...
                toUnselect.append(s)
            if s == shape:
                print index
                break
            if isinstance(s,ogl.BitmapShape):
                index = index + 1

        home.editGeneric(index)


class MyLineEvtHandler(ogl.ShapeEvtHandler):
    
    def __init__(self, log, frame):
        ogl.ShapeEvtHandler.__init__(self)
        self.log = log
        self.statbarFrame = frame
        
    def UpdateStatusBar(self, shape):
        x, y = shape.GetX(), shape.GetY()
        width, height = shape.GetBoundingBoxMax()
        #self.statbarFrame.SetStatusText("Pos: (%d, %d)  Size: (%d, %d)" %
                                       # (x, y, width, height))


    def OnLeftClick(self, x, y, keys=0, attachment=0):
        shape = self.GetShape()
        canvas = shape.GetCanvas()
        dc = wx.ClientDC(canvas)
        canvas.PrepareDC(dc)

        if shape.Selected():
            shape.Select(False, dc)
            #canvas.Redraw(dc)
            canvas.Refresh(False)
        else:
            redraw = False
            shapeList = canvas.GetDiagram().GetShapeList()
            toUnselect = []
            
            index = 0
            for s in shapeList:
                if s.Selected():
                    # If we unselect it now then some of the objects in
                    # shapeList will become invalid (the control points are
                    # shapes too!) and bad things will happen...
                    toUnselect.append(s)
             
            shape.Select(True, dc)

            if toUnselect:
                for s in toUnselect:
                    s.Select(False, dc)

                ##canvas.Redraw(dc)
                canvas.Refresh(False)

        self.UpdateStatusBar(shape)

        
    def OnEndDragLeft(self, x, y, keys=0, attachment=0):
        pass


    def OnSizingEndDragLeft(self, pt, x, y, keys, attch):
        pass


    def OnMovePost(self, dc, x, y, oldX, oldY, display):
        pass

    def OnKeyDown(self, evt):
        print "tecla presionada"
        if self.logKeyDn:
            self.GetParent().keylog.LogKeyEvent("KeyDown", evt)
        if self.callSkip:
            evt.Skip()

    def OnRightClick(self,  x, y, keys=0, attachment=0):
        wx.MessageBox("Aca deberia permitir cancelar","Error")
        
    def OnLeftDoubleClick(self,x, y, keys, attachment):
        pass
        

keyMap = {
    wx.WXK_BACK : "WXK_BACK",
    wx.WXK_TAB : "WXK_TAB",
    wx.WXK_RETURN : "WXK_RETURN",
    wx.WXK_ESCAPE : "WXK_ESCAPE",
    wx.WXK_SPACE : "WXK_SPACE",
    wx.WXK_DELETE : "WXK_DELETE",
    wx.WXK_START : "WXK_START",
    wx.WXK_LBUTTON : "WXK_LBUTTON",
    wx.WXK_RBUTTON : "WXK_RBUTTON",
    wx.WXK_CANCEL : "WXK_CANCEL",
    wx.WXK_MBUTTON : "WXK_MBUTTON",
    wx.WXK_CLEAR : "WXK_CLEAR",
    wx.WXK_SHIFT : "WXK_SHIFT",
    wx.WXK_ALT : "WXK_ALT",
    wx.WXK_CONTROL : "WXK_CONTROL",
    wx.WXK_MENU : "WXK_MENU",
    wx.WXK_PAUSE : "WXK_PAUSE",
    wx.WXK_CAPITAL : "WXK_CAPITAL",
    #wx.WXK_PRIOR : "WXK_PRIOR",
    #wx.WXK_NEXT : "WXK_NEXT",
    wx.WXK_END : "WXK_END",
    wx.WXK_HOME : "WXK_HOME",
    wx.WXK_LEFT : "WXK_LEFT",
    wx.WXK_UP : "WXK_UP",
    wx.WXK_RIGHT : "WXK_RIGHT",
    wx.WXK_DOWN : "WXK_DOWN",
    wx.WXK_SELECT : "WXK_SELECT",
    wx.WXK_PRINT : "WXK_PRINT",
    wx.WXK_EXECUTE : "WXK_EXECUTE",
    wx.WXK_SNAPSHOT : "WXK_SNAPSHOT",
    wx.WXK_INSERT : "WXK_INSERT",
    wx.WXK_HELP : "WXK_HELP",
    wx.WXK_NUMPAD0 : "WXK_NUMPAD0",
    wx.WXK_NUMPAD1 : "WXK_NUMPAD1",
    wx.WXK_NUMPAD2 : "WXK_NUMPAD2",
    wx.WXK_NUMPAD3 : "WXK_NUMPAD3",
    wx.WXK_NUMPAD4 : "WXK_NUMPAD4",
    wx.WXK_NUMPAD5 : "WXK_NUMPAD5",
    wx.WXK_NUMPAD6 : "WXK_NUMPAD6",
    wx.WXK_NUMPAD7 : "WXK_NUMPAD7",
    wx.WXK_NUMPAD8 : "WXK_NUMPAD8",
    wx.WXK_NUMPAD9 : "WXK_NUMPAD9",
    wx.WXK_MULTIPLY : "WXK_MULTIPLY",
    wx.WXK_ADD : "WXK_ADD",
    wx.WXK_SEPARATOR : "WXK_SEPARATOR",
    wx.WXK_SUBTRACT : "WXK_SUBTRACT",
    wx.WXK_DECIMAL : "WXK_DECIMAL",
    wx.WXK_DIVIDE : "WXK_DIVIDE",
    wx.WXK_F1 : "WXK_F1",
    wx.WXK_F2 : "WXK_F2",
    wx.WXK_F3 : "WXK_F3",
    wx.WXK_F4 : "WXK_F4",
    wx.WXK_F5 : "WXK_F5",
    wx.WXK_F6 : "WXK_F6",
    wx.WXK_F7 : "WXK_F7",
    wx.WXK_F8 : "WXK_F8",
    wx.WXK_F9 : "WXK_F9",
    wx.WXK_F10 : "WXK_F10",
    wx.WXK_F11 : "WXK_F11",
    wx.WXK_F12 : "WXK_F12",
    wx.WXK_F13 : "WXK_F13",
    wx.WXK_F14 : "WXK_F14",
    wx.WXK_F15 : "WXK_F15",
    wx.WXK_F16 : "WXK_F16",
    wx.WXK_F17 : "WXK_F17",
    wx.WXK_F18 : "WXK_F18",
    wx.WXK_F19 : "WXK_F19",
    wx.WXK_F20 : "WXK_F20",
    wx.WXK_F21 : "WXK_F21",
    wx.WXK_F22 : "WXK_F22",
    wx.WXK_F23 : "WXK_F23",
    wx.WXK_F24 : "WXK_F24",
    wx.WXK_NUMLOCK : "WXK_NUMLOCK",
    wx.WXK_SCROLL : "WXK_SCROLL",
    wx.WXK_PAGEUP : "WXK_PAGEUP",
    wx.WXK_PAGEDOWN : "WXK_PAGEDOWN",
    wx.WXK_NUMPAD_SPACE : "WXK_NUMPAD_SPACE",
    wx.WXK_NUMPAD_TAB : "WXK_NUMPAD_TAB",
    wx.WXK_NUMPAD_ENTER : "WXK_NUMPAD_ENTER",
    wx.WXK_NUMPAD_F1 : "WXK_NUMPAD_F1",
    wx.WXK_NUMPAD_F2 : "WXK_NUMPAD_F2",
    wx.WXK_NUMPAD_F3 : "WXK_NUMPAD_F3",
    wx.WXK_NUMPAD_F4 : "WXK_NUMPAD_F4",
    wx.WXK_NUMPAD_HOME : "WXK_NUMPAD_HOME",
    wx.WXK_NUMPAD_LEFT : "WXK_NUMPAD_LEFT",
    wx.WXK_NUMPAD_UP : "WXK_NUMPAD_UP",
    wx.WXK_NUMPAD_RIGHT : "WXK_NUMPAD_RIGHT",
    wx.WXK_NUMPAD_DOWN : "WXK_NUMPAD_DOWN",
    #wx.WXK_NUMPAD_PRIOR : "WXK_NUMPAD_PRIOR",
    wx.WXK_NUMPAD_PAGEUP : "WXK_NUMPAD_PAGEUP",
    #wx.WXK_NUMPAD_NEXT : "WXK_NUMPAD_NEXT",
    wx.WXK_NUMPAD_PAGEDOWN : "WXK_NUMPAD_PAGEDOWN",
    wx.WXK_NUMPAD_END : "WXK_NUMPAD_END",
    wx.WXK_NUMPAD_BEGIN : "WXK_NUMPAD_BEGIN",
    wx.WXK_NUMPAD_INSERT : "WXK_NUMPAD_INSERT",
    wx.WXK_NUMPAD_DELETE : "WXK_NUMPAD_DELETE",
    wx.WXK_NUMPAD_EQUAL : "WXK_NUMPAD_EQUAL",
    wx.WXK_NUMPAD_MULTIPLY : "WXK_NUMPAD_MULTIPLY",
    wx.WXK_NUMPAD_ADD : "WXK_NUMPAD_ADD",
    wx.WXK_NUMPAD_SEPARATOR : "WXK_NUMPAD_SEPARATOR",
    wx.WXK_NUMPAD_SUBTRACT : "WXK_NUMPAD_SUBTRACT",
    wx.WXK_NUMPAD_DECIMAL : "WXK_NUMPAD_DECIMAL",
    wx.WXK_NUMPAD_DIVIDE : "WXK_NUMPAD_DIVIDE",

    wx.WXK_WINDOWS_LEFT : "WXK_WINDOWS_LEFT",
    wx.WXK_WINDOWS_RIGHT : "WXK_WINDOWS_RIGHT",
    wx.WXK_WINDOWS_MENU : "WXK_WINDOWS_MENU",

    wx.WXK_COMMAND : "WXK_COMMAND",

    wx.WXK_SPECIAL1 : "WXK_SPECIAL1",
    wx.WXK_SPECIAL2 : "WXK_SPECIAL2",
    wx.WXK_SPECIAL3 : "WXK_SPECIAL3",
    wx.WXK_SPECIAL4 : "WXK_SPECIAL4",
    wx.WXK_SPECIAL5 : "WXK_SPECIAL5",
    wx.WXK_SPECIAL6 : "WXK_SPECIAL6",
    wx.WXK_SPECIAL7 : "WXK_SPECIAL7",
    wx.WXK_SPECIAL8 : "WXK_SPECIAL8",
    wx.WXK_SPECIAL9 : "WXK_SPECIAL9",
    wx.WXK_SPECIAL10 : "WXK_SPECIAL10",
    wx.WXK_SPECIAL11 : "WXK_SPECIAL11",
    wx.WXK_SPECIAL12 : "WXK_SPECIAL12",
    wx.WXK_SPECIAL13 : "WXK_SPECIAL13",
    wx.WXK_SPECIAL14 : "WXK_SPECIAL14",
    wx.WXK_SPECIAL15 : "WXK_SPECIAL15",
    wx.WXK_SPECIAL16 : "WXK_SPECIAL16",
    wx.WXK_SPECIAL17 : "WXK_SPECIAL17",
    wx.WXK_SPECIAL18 : "WXK_SPECIAL18",
    wx.WXK_SPECIAL19 : "WXK_SPECIAL19",
    wx.WXK_SPECIAL2 : "WXK_SPECIAL2",
}
