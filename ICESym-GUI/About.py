import wx
import wx.html

class About(wx.Dialog):
    def __init__(self, parent):
		wx.Dialog.__init__(self, parent, -1, 'About ICESym-GUI',
				          size=(450, 350) )
		html = wx.html.HtmlWindow(self)
		html.LoadPage('html/about.html')
		button = wx.Button(self, wx.ID_OK, "Close")
		sizer = wx.BoxSizer(wx.VERTICAL)
		sizer.Add(html, 1, wx.EXPAND|wx.ALL, 5)
		sizer.Add(button, 0, wx.ALIGN_CENTER|wx.ALL, 5)
		self.SetSizer(sizer)
		self.Layout()

