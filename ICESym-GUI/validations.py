"""
	This file sets all validator for each input in each form
"""
import wx
import  string
string.digits=string.digits+'e./+-'
string_list = string.digits+','

class numberValidator(wx.PyValidator):
	
	def __init__(self):
		wx.PyValidator.__init__(self)			
	
	def Clone(self):
		"""
		Note that every validator must implement the Clone() method.
		"""
		return numberValidator()

	def validateNumber(self, text):
		number=1
		if len(text) == 0:
			#wx.MessageBox("This field must contain some text!", "Error")
			number=0
		else:
			for x in text:
				if x not in string.digits:
					number=0
					break
		return number

	def validateList(self, text):
		res=1
		if len(text) == 0:
			#wx.MessageBox("This field must contain some text!", "Error")
			res = 0
		else:
			for x in text:
				if x not in string_list:
					res = 0
					break
		return res
	
	def Validate(self, win, tipo):
		textCtrl = self.GetWindow()
		text = textCtrl.GetValue()
		
		ret = 1
		if(tipo=='number'):
			ret = self.validateNumber(text)
		if(tipo=='list'):
			ret = self.validateList(text)

		if(ret==1):
			textCtrl.SetBackgroundColour(wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
			textCtrl.Refresh()
			return True
		else:
			textCtrl.SetBackgroundColour("pink")
			textCtrl.SetFocus()
			textCtrl.Refresh()
			return False
		
	def TransferToWindow(self):
		print "en open"
		return True
	
	def TransferFromWindow(self):
		print "en close"
		return True	
