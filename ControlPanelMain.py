import wx

class ControlPanel(wx.Panel):
	def __init__(self, parent, **kwargs):
		global coord1, coord2, Data
		wx.Panel.__init__(self, parent, **kwargs)
		self.grid = wx.GridBagSizer(hgap=5, vgap=5)	#Helps automatically place the controls

		#Linecolor Control
		self.linecolorlist = ['white','red','blue','green']
		self.lctext = wx.StaticText(self, label="Line Color:")
		self.grid.Add(self.lctext, pos=(0,0), flag=wx.ALL)
		self.editlc = wx.ComboBox(self, size=(100, -1), choices=self.linecolorlist, style=wx.TE_PROCESS_ENTER)
		self.grid.Add(self.editlc, pos=(1,0), span=(1,2), flag=wx.ALL)

		# 2D Save Plot Button
		self.save2Dbutton = wx.Button(self, label="Save PES")
		self.grid.Add(self.save2Dbutton, pos=(2,0), flag=wx.ALL)

		# Sliders to set cbrange of 2D plot
		self.cminText = wx.StaticText(self, label="Color Bar Minimum:")
		self.grid.Add(self.cminText, pos=(3,0), flag=wx.ALL, span=(1,2))
		self.cmaxText = wx.StaticText(self, label="Color Bar Maximum:")
		self.grid.Add(self.cmaxText, pos=(5,0), flag=wx.ALL)
		self.cmin = self.GetGrandParent().GetParent().panel1.cmin
		self.cmax = self.GetGrandParent().GetParent().panel1.cmax
		self.cbminSlider = wx.Slider(self, -1, 0, 0, 100, size = (100,-1), style=wx.HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
		self.cbminSlider.SetTickFreq(5,1)
		self.cbmaxSlider = wx.Slider(self, -1, 100, 0, 100, size = (100,-1), style=wx.HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
		self.cbmaxSlider.SetTickFreq(5,1)
		self.SetSizer(self.grid)
		self.grid.Add(self.cbminSlider, pos=(4,0), span=(1,2), flag=wx.ALL)
		self.grid.Add(self.cbmaxSlider, pos=(6,0), span=(1,2), flag=wx.ALL)

		# Manually input the coordinates of the drawn line
		self.c1xEntry = wx.TextCtrl(self,-1,"%.2f" % coord1.x,size=(50,-1))
		self.c1yEntry = wx.TextCtrl(self,-1,"%.2f" % coord1.y,size=(50,-1))
		self.c2xEntry = wx.TextCtrl(self,-1,"%.2f" % coord2.x,size=(50,-1))
		self.c2yEntry = wx.TextCtrl(self,-1,"%.2f" % coord2.y,size=(50,-1))

		self.grid.Add(self.c1xEntry, pos = (7,0), flag = wx.ALL)
		self.grid.Add(self.c1yEntry, pos = (7,1), flag = wx.ALL)
		self.grid.Add(self.c2xEntry, pos = (8,0), flag = wx.ALL)
		self.grid.Add(self.c2yEntry, pos = (8,1), flag = wx.ALL)

		# Bind the previously defined controls to the proper functions
		self.Bind(wx.EVT_COMBOBOX, self.LineReColor, self.editlc)
		self.Bind(wx.EVT_TEXT_ENTER, self.LineReColor, self.editlc)
		self.Bind(wx.EVT_BUTTON, self.Save2DPlot, self.save2Dbutton)
		self.Bind(wx.EVT_SLIDER, self.CBRUpdate)

	def CBRUpdate(self,e):
		self.GetGrandParent().GetParent().panel1.cmax = max(Data.flatten())*(float(self.cbmaxSlider.GetValue())/100.0)
		self.GetGrandParent().GetParent().panel1.cmin = max(Data.flatten())*(float(self.cbminSlider.GetValue())/100.0)
		self.GetGrandParent().GetParent().panel1.Draw()

	def Save2DPlot(self, e): #Save the 2D PES plot to a file
		dlg = wx.FileDialog(self, "Save 2D Plot", dirname, "", "*.*", wx.SAVE)
		if dlg.ShowModal() == wx.ID_OK:
			self.filename = dlg.GetFilename()
			self.dirname = dlg.GetDirectory()
			self.plot2Dfile = os.path.join(self.dirname, self.filename)
			self.GetGrandParent().GetParent().panel1.figure.savefig(self.plot2Dfile,dpi=300)
		dlg.Destroy()

	def LineReColor(self,e): #Recolor the line with the desired color
		global linecolor
		linecolor = str(self.editlc.GetValue().strip().lower())
		self.GetGrandParent().GetParent().panel1.AddLine(coord1, coord2)