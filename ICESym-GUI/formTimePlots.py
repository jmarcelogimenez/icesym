
import wx

class timePlots(wx.Dialog):
    def __init__(self, *args, **kwds):
        # begin wxGlade: timePlots.__init__
        kwds["style"] = wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER
        wx.Dialog.__init__(self, *args, **kwds)
        self.panel_buttons = wx.Panel(self, -1)
        self.panel_configure = wx.ScrolledWindow(self, -1, style=wx.TAB_TRAVERSAL)

        self.accept = wx.Button(self.panel_buttons, wx.ID_OK, "")
        self.cancel = wx.Button(self.panel_buttons, wx.ID_CANCEL, "")
 
        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_BUTTON, self.onAccept, self.accept)
        # end wxGlade

    def __set_properties(self):
        # begin wxGlade: timePlots.__set_properties
        self.SetTitle("Time Plots")
        self.SetSize((200, 300))
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: timePlots.__do_layout
        grid_sizer = wx.FlexGridSizer(2, 1, 0, 0)
        grid_sizer_configure = wx.FlexGridSizer(5, 2, 0, 0)

        self.panel_configure.SetSizer(grid_sizer_configure)
        
        grid_sizer_buttons = wx.FlexGridSizer(1, 2, 0, 0)
        grid_sizer_buttons.Add(self.accept, 1, wx.EXPAND, 0) 
        grid_sizer_buttons.Add(self.cancel, 1, wx.EXPAND, 0)
        self.panel_buttons.SetSizer(grid_sizer_buttons)

        grid_sizer.Add(self.panel_configure, 1, wx.EXPAND, 0)
        grid_sizer.Add(self.panel_buttons, 1, wx.EXPAND, 0)
        self.SetSizer(grid_sizer)
        self.Layout()
        # end wxGlade

    def onAccept(self, event): # wxGlade: timePlots.<event_handler>
        print "Event handler `onAccept' not implemented!"
        event.Skip()

