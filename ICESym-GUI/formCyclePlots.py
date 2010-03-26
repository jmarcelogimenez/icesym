import string
import wx
from unitsData import *

class cyclePlots(wx.Dialog):
    data = dict()
    labels = dict()
    listNdofA = ['density', 'velocity', 'pressure']
    listNdofB = ['density', 'pressure', 'temperature']
    cylExtras= ['convective heat-transfer coeff','radiactive heat-transfer coeff','convective heat-transfer rate','radiactive heat-transfer rate','burned mass fraction','burned mass fraction rate','mass flow rate trought intake port','mass flow rate trought exhaust port','mass of fuel','mass of air', 'mass of residual gas','total heat-transfer rate','fuel chemical energy release','torque']
    tankExtras = ['mass flow rate','enthalpy flow rate','mass','heat-transfer rate']
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
		# begin wxGlade: anglePlots.__set_properties
		self.SetTitle("Cycle Plots")
		self.SetSize((275, 350))
		for i in self.data:
			self.data[i].SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
			if not(i=="element"):
				self.labels[i].SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
		self.accept.SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
		self.cancel.SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: cyclePlots.__do_layout
		grid_sizer = wx.FlexGridSizer(2, 1, 10, 10)
		grid_sizer_buttons = wx.FlexGridSizer(1, 2, 0, 0)
		grid_sizer_configure = wx.FlexGridSizer(5, 2, 0, 0)
		grid_sizer_element = wx.FlexGridSizer(3, 2, 0, 0)
		
		grid_sizer_buttons.Add(self.accept, 1, wx.ALIGN_CENTER_HORIZONTAL, 0) 
		grid_sizer_buttons.Add(self.cancel, 1, wx.ALIGN_CENTER_HORIZONTAL, 0)
		self.panel_buttons.SetSizer(grid_sizer_buttons)

		grid_sizer_configure.Add(self.labels['rpm'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_configure.Add(self.data['rpm'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_configure.Add(self.labels['label'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_configure.Add(self.data['label'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_configure.Add(self.labels['figure'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_configure.Add(self.data['figure'], 0, wx.ALIGN_CENTER_VERTICAL, 0)	
		grid_sizer_configure.Add(self.data['element'], 0, wx.ALIGN_CENTER_VERTICAL, 0)	

		grid_sizer_element.Add(self.labels['element_number'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_element.Add(self.data['element_number'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_element.Add(self.labels['element_nnod'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_element.Add(self.data['element_nnod'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_element.Add(self.labels['element_ndof'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_element.Add(self.data['element_ndof'], 0, wx.ALIGN_CENTER_VERTICAL, 0)

		grid_sizer_configure.Add(grid_sizer_element, 0, wx.ALIGN_CENTER_VERTICAL, 0)	
		grid_sizer_configure.Add(self.labels['units'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		grid_sizer_configure.Add(self.data['units'], 0, wx.ALIGN_CENTER_VERTICAL, 0)
		self.panel_configure.SetSizer(grid_sizer_configure)

		grid_sizer.Add(self.panel_configure, 1, wx.ALIGN_CENTER_HORIZONTAL, 0)
		grid_sizer.Add(self.panel_buttons, 1, wx.ALIGN_CENTER_HORIZONTAL, 0)
		self.panel_gral.SetSizer(grid_sizer)
		self.SetSizer(grid_sizer)
		self.Layout()
        # end wxGlade

    def setParameters(self,plotData):

		self.labels['rpm'] = wx.StaticText(self.panel_configure, -1, "Choose RPM: ")
		#self.data['rpm'] = wx.Choice(self.panel_configure,choices=plotData['rpm'])
		#self.data['rpm'].SetStringSelection(str(plotData['rpm'][-1]))
		self.data['rpm'] = wx.CheckListBox(self.panel_configure,choices=plotData['rpm'],size=(100,80))
		self.data['rpm'].Check(len(plotData['rpm'])-1, True)

		self.labels['label'] = wx.StaticText(self.panel_configure, -1, "Label: ")
		self.data['label'] = wx.TextCtrl(self.panel_configure, -1, "")

		self.labels['figure'] = wx.StaticText(self.panel_configure, -1, "Choose Figure: ")
		self.data['figure'] = wx.Choice(self.panel_configure,choices=plotData['figure'])
		self.data['figure'].SetStringSelection(str(plotData['figure'][0]))

		self.data['element'] = wx.RadioBox(self.panel_configure, -1, "", choices=["Globals", "Cylinders", "Tubes", "Tanks", "Junctions"], majorDimension=0, style=wx.RA_SPECIFY_ROWS)
		self.labels['element_number'] = wx.StaticText(self.panel_configure, -1, "Choose Element: ")
		self.data['element_number'] = wx.Choice(self.panel_configure, -1, choices=[])
		self.labels['element_nnod'] = wx.StaticText(self.panel_configure, -1, "Choose Node: ")
		self.data['element_nnod'] = wx.Choice(self.panel_configure, -1, choices=[])
		self.labels['element_ndof'] = wx.StaticText(self.panel_configure, -1, "Choose Variable: ")
		self.data['element_ndof'] = wx.Choice(self.panel_configure, -1, choices=[])
		self.labels['units'] = wx.StaticText(self.panel_configure, -1, "Choose Units: ")
		self.data['units'] = wx.Choice(self.panel_configure, -1, choices=[])
		self.infoData = plotData
		self.Bind(wx.EVT_RADIOBOX, self.onChangeElement, self.data['element'])
		self.Bind(wx.EVT_CHOICE, self.onChangeElementNumber, self.data['element_number'])
		self.Bind(wx.EVT_CHOICE, self.onChangeElementNode, self.data['element_nnod'])
		self.Bind(wx.EVT_CHOICE, self.onChangeElementNdof, self.data['element_ndof'])

		self.__set_properties()
		self.__do_layout()


    def onAccept(self, event): # wxGlade: anglePlots.<event_handler>
		can_out = 1		
		for d in self.data:
			if isinstance(self.data[d],wx.Choice) or isinstance(self.data[d],wx.RadioBox):
				if self.data[d].GetSelection() == -1 and self.data[d].IsEnabled():
					self.data[d].SetBackgroundColour("pink")
					can_out = 0
				else:
					self.data[d].SetBackgroundColour(wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
					self.parsedData[d] = self.data[d].GetStringSelection()
					if d == "element_nnod": #por lass dudas elige una valve (agrega intake o exhaust)
						l = string.split(str(self.parsedData[d]))
						self.parsedData[d] = l[0]
					if d in ["element_ndof", "element", "element_number"]: #no necesito la palabra, sino el numero
						self.parsedData[d+"_str"] = self.parsedData[d]
						self.parsedData[d] = self.data[d].GetSelection()
			elif isinstance(self.data[d],wx.CheckListBox):
				self.parsedData[d] = []
				for i in range(len(self.infoData[d])):
					if self.data[d].IsChecked(i):
						self.parsedData[d].append(self.infoData[d][i])
				if self.parsedData[d] == []:
					self.data[d].SetBackgroundColour("pink")
					can_out = 0
				else:
					self.data[d].SetBackgroundColour(wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
			else:
				if self.data[d].GetValue() == '':
					self.data[d].SetBackgroundColour("pink")
					can_out = 0
				else:
					self.data[d].SetBackgroundColour(wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
					self.parsedData[d] = self.data[d].GetValue()
	
		if can_out == 1:
			self.EndModal(wx.ID_OK)		
		else:
			wx.MessageBox("Please Complete all data", "Error")

    def onChangeElement(self, event):
		#Clear() limpia la lista del choice
		#AppendItems(list) agrega la lista entera
		elem = self.data['element'].GetSelection()
		if elem == 0:
			self.data['element_number'].Clear()
			self.data['element_nnod'].Clear()
			self.data['element_ndof'].Clear()
			#self.data['element_number'].AppendItems(self.infoData['Globals']['labels'])
			self.labels['element_number'].SetLabel("Global Variable: ")
			self.labels['element_nnod'].SetLabel("- - -")
			self.labels['element_ndof'].SetLabel("- - -")
			self.data['element_nnod'].Enable(False)
			self.data['element_ndof'].Enable(False)
		elif elem == 1:
			self.data['element_number'].Clear()
			self.data['element_nnod'].Clear()
			self.data['element_ndof'].Clear()
			self.data['element_number'].AppendItems(self.infoData['Cylinders']['labels'])
			self.labels['element_number'].SetLabel("Choose Element: ")
			self.labels['element_nnod'].SetLabel("Choose Node: ")
			self.labels['element_ndof'].SetLabel("Choose Variable: ")
			self.data['element_nnod'].Enable(True)
			self.data['element_ndof'].Enable(True)
		elif elem == 2:
			self.data['element_number'].Clear()
			self.data['element_nnod'].Clear()
			self.data['element_ndof'].Clear()
			self.data['element_number'].AppendItems(self.infoData['Tubes']['labels'])
			self.labels['element_number'].SetLabel("Choose Element: ")
			self.labels['element_nnod'].SetLabel("Choose Node: ")
			self.labels['element_ndof'].SetLabel("Choose Variable: ")
			self.data['element_nnod'].Enable(True)
			self.data['element_ndof'].Enable(True)
		elif elem == 3:
			self.data['element_number'].Clear()
			self.data['element_nnod'].Clear()
			self.data['element_ndof'].Clear()
			self.data['element_number'].AppendItems(self.infoData['Tanks']['labels'])
			self.labels['element_number'].SetLabel("Choose Element: ")
			self.labels['element_nnod'].SetLabel("Choose Node: ")
			self.labels['element_ndof'].SetLabel("Choose Variable: ")
			self.data['element_nnod'].Enable(True)
			self.data['element_ndof'].Enable(True)
		elif elem == 4:
			self.data['element_number'].Clear()
			self.data['element_nnod'].Clear()
			self.data['element_ndof'].Clear()
			self.data['element_number'].AppendItems(self.infoData['Junctions']['labels'])
			self.labels['element_number'].SetLabel("Choose Element: ")
			self.labels['element_nnod'].SetLabel("Choose Node: ")
			self.labels['element_ndof'].SetLabel("Choose Variable: ")
			self.data['element_nnod'].Enable(True)
			self.data['element_ndof'].Enable(True)

    def onChangeElementNumber(self, event):
		elem = self.data['element'].GetSelection()
		number = int(self.data['element_number'].GetSelection())
		if elem == 0:
			print "no globals data"
		if elem == 1:
			self.data['element_ndof'].Clear()
			self.data['element_nnod'].Clear()
			self.data['element_nnod'].AppendItems(self.infoData['Cylinders']['histos'][number] + self.infoData['Cylinders']['valves_intake'][number] + self.infoData['Cylinders']['valves_exhaust'][number] )
		elif elem == 2:
			self.data['element_nnod'].Clear()
			self.data['element_ndof'].Clear()
			self.data['element_nnod'].AppendItems(self.infoData['Tubes']['histos'][number])
		elif elem == 3:
			self.data['element_nnod'].Clear()
			self.data['element_ndof'].Clear()
			self.data['element_nnod'].AppendItems(self.infoData['Tanks']['histos'][number])
		elif elem == 4:
			self.data['element_nnod'].Clear()
			self.data['element_ndof'].Clear()
			self.data['element_nnod'].AppendItems(self.infoData['Junctions']['histos'][number])
	
    def onChangeElementNode(self, event):
		elem = self.data['element'].GetSelection()
		number = self.data['element_number'].GetSelection()
		nnod = self.data['element_nnod'].GetSelection()

		if elem == 1:
			self.data['element_ndof'].Clear()
			if not(str(nnod) in self.infoData['Cylinders']['histos'][number]):
				self.data['element_ndof'].AppendItems(self.listNdofA)
			else:
				choices = []			
				if self.infoData['Cylinders']['extras'][number] == 1:
					choices = self.listNdofB + self.cylExtras
				else:
					choices = self.listNdofB
				self.data['element_ndof'].AppendItems(choices)
		elif elem == 2:
			self.data['element_ndof'].Clear()
			self.data['element_ndof'].AppendItems(self.listNdofA)
		elif elem == 3:
			self.data['element_ndof'].Clear()
			choices = []			
			if self.infoData['Tanks']['extras'][number] == 1:
				choices = self.listNdofB + self.tankExtras
			else:
				choices = self.listNdofB
			self.data['element_ndof'].AppendItems(choices)
		elif elem == 4:
			self.data['element_ndof'].Clear()
			self.data['element_ndof'].AppendItems(self.listNdofA)

    def onChangeElementNdof(self, event):
		ndof = self.data['element_ndof'].GetStringSelection()
		self.data['units'].Clear()
		self.data['units'].AppendItems(possibleUnits[ndof])
		if len(possibleUnits[ndof]) > 0:
			self.data['units'].SetStringSelection(str(possibleUnits[ndof][0]))
		else:
			self.data['units'].Enable(False)

