
import wx

class spacePlots(wx.Dialog):
    data = dict()
    labels = dict()
    listNdofA = ['density', 'velocity', 'pressure']
    listNdofB = ['density', 'pressure','temperature']
    parsedData = dict()
    def __init__(self, *args, **kwds):
		# begin wxGlade: anglePlots.__init__
		kwds["style"] = wx.DEFAULT_DIALOG_STYLE
		wx.Dialog.__init__(self, *args, **kwds)
		
		self.panel_gral = wx.Panel(self, -1)
		self.panel_buttons = wx.Panel(self, -1)
		self.panel_configure = wx.ScrolledWindow(self, -1, style=wx.TAB_TRAVERSAL)

		self.accept = wx.Button(self.panel_buttons, wx.ID_OK, "")
		self.cancel = wx.Button(self.panel_buttons, wx.ID_CANCEL, "")
 
		#self.__set_properties()
		#self.__do_layout()

		self.Bind(wx.EVT_BUTTON, self.onAccept, self.accept)
        # end wxGlade


    def __set_properties(self):
        # begin wxGlade: spacePlots.__set_properties
		self.SetTitle("Space Plots")
		self.SetSize((400, 225))
		for i in self.data:
			self.data[i].SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
			if not(i=="element"):
				self.labels[i].SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
		self.labels['text'].SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
		self.accept.SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
		self.cancel.SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: spacePlots.__do_layout
		grid_sizer = wx.FlexGridSizer(2, 1, 10, 10)
		grid_sizer_buttons = wx.FlexGridSizer(1, 2, 0, 0)
		grid_sizer_configure = wx.FlexGridSizer(4, 2, 0, 0)
		grid_sizer_element = wx.FlexGridSizer(2, 2, 0, 0)
		
		grid_sizer_buttons.Add(self.accept, 1, wx.ALIGN_CENTER_HORIZONTAL, 0) 
		grid_sizer_buttons.Add(self.cancel, 1, wx.ALIGN_CENTER_HORIZONTAL, 0)
		self.panel_buttons.SetSizer(grid_sizer_buttons)

		grid_sizer_configure.Add(self.labels['rpm'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_configure.Add(self.data['rpm'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_configure.Add(self.labels['time'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_configure.Add(self.data['time'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_configure.Add(self.labels['label'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_configure.Add(self.data['label'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_configure.Add(self.labels['figure'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_configure.Add(self.data['figure'], 0, wx.ALIGN_CENTER_VERTICAL, 0)	

		grid_sizer_element.Add(self.labels['element_number'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_element.Add(self.data['element_number'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_element.Add(self.labels['element_ndof'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_element.Add(self.data['element_ndof'], 0, wx.ALIGN_CENTER_VERTICAL, 0)

		grid_sizer_configure.Add(grid_sizer_element, 0, wx.ALIGN_CENTER_VERTICAL, 0)	
		grid_sizer_configure.Add(self.labels['text'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		self.panel_configure.SetSizer(grid_sizer_configure)

		grid_sizer.Add(self.panel_configure, 1, wx.ALIGN_CENTER_HORIZONTAL, 0)
		grid_sizer.Add(self.panel_buttons, 1, wx.ALIGN_CENTER_HORIZONTAL, 0)
		self.panel_gral.SetSizer(grid_sizer)
		self.SetSizer(grid_sizer)
		self.Layout()
        # end wxGlade

    def setParameters(self,plotData):

		self.labels['rpm'] = wx.StaticText(self.panel_configure, -1, "Choose RPM: ")
		self.data['rpm'] = wx.Choice(self.panel_configure,choices=plotData['rpm'])
		self.data['rpm'].SetStringSelection(str(plotData['rpm'][-1]))

		self.labels['time'] = wx.StaticText(self.panel_configure, -1, "Choose Time (max " + str(plotData['time'][-1]) + ") : ")
		self.data['time'] = wx.TextCtrl(self.panel_configure,-1,str(plotData['time'][-1]))
	
		self.labels['label'] = wx.StaticText(self.panel_configure, -1, "Label: ")
		self.data['label'] = wx.TextCtrl(self.panel_configure, -1, "")

		self.labels['figure'] = wx.StaticText(self.panel_configure, -1, "Choose Figure: ")
		self.data['figure'] = wx.Choice(self.panel_configure,choices=plotData['figure'])
		self.data['figure'].SetStringSelection(str(plotData['figure'][0]))

		self.labels['element_number'] = wx.StaticText(self.panel_configure, -1, "Choose Element: ")
		self.data['element_number'] = wx.Choice(self.panel_configure, -1, choices=plotData['Tubes']['labels'])
		self.labels['element_ndof'] = wx.StaticText(self.panel_configure, -1, "Choose Variable: ")
		self.data['element_ndof'] = wx.Choice(self.panel_configure, -1, choices=self.listNdofA)
		self.labels['text'] = wx.StaticText(self.panel_configure, -1, "Space Plots available only for Tubes \n with all nodes included in the story: ")

		self.infoData = plotData
		self.Bind(wx.EVT_CHOICE, self.onChangeRPM, self.data['rpm'])
		
		self.__set_properties()
		self.__do_layout()

    def onAccept(self, event): # wxGlade: spacePlots.<event_handler>
		can_out = 1		
		for d in self.data:
			if isinstance(self.data[d],wx.Choice) or isinstance(self.data[d],wx.RadioBox):
				if self.data[d].GetSelection() == -1 and self.data[d].IsEnabled():
					self.data[d].SetBackgroundColour("pink")
					can_out = 0
				else:
					self.data[d].SetBackgroundColour(wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
					self.parsedData[d] = self.data[d].GetStringSelection()
					if d in ["element_ndof", "element_number"]: #no necesito la palabra, sino el numero
						self.parsedData[d+"_str"] = self.parsedData[d]
						self.parsedData[d] = self.data[d].GetSelection()
			else:
				if self.data[d].GetValue() == '':
					self.data[d].SetBackgroundColour("pink")
					can_out = 0
				if d=="time":	
					tm = float(self.data[d].GetValue())
					rpm = self.data['rpm'].GetSelection()
					maxTime = self.infoData['time'][rpm]
					if tm < 0 or tm > maxTime:
						self.data[d].SetBackgroundColour("pink")
						can_out = 0
					else:
						self.data[d].SetBackgroundColour(wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
						self.parsedData[d] = tm
				else:
					self.data[d].SetBackgroundColour(wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
					self.parsedData[d] = self.data[d].GetValue()
	
		if can_out == 1:
			self.EndModal(wx.ID_OK)		
		else:
			wx.MessageBox("Please Complete all data or correct the error in the time", "Error")


    def onChangeRPM(self, event):
		rpm = self.data['rpm'].GetSelection()
		maxTime = self.infoData['time'][rpm]
		self.labels['time'].SetLabel("Choose Time (max " + str(maxTime) + ") : ")
		self.data['time'].SetValue(str(maxTime))
		

