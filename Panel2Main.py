#######################################
#
# Panel containing the 1D data slice
#
######################################
import wx

class PlotPanel2(wx.Panel):
	def __init__(self, parent, **kwargs):
		wx.Panel.__init__(self, parent, **kwargs)
		self.figure = mpl.figure.Figure()
		self.axes = self.figure.add_subplot(111)
		self.canvas = FigureCanvas(self, -1, self.figure)
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		self.sizer.Add(self.canvas, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER)
		self.SetSizer(self.sizer)
		self.Fit()

	def Draw(self):
		global XSlice, YSlice, SSlice, DataSlice
		self.axes.cla()
		XSlice = np.linspace(coord1.x,coord2.x,200)
		YSlice = np.linspace(coord1.y,coord2.y,200)
		self.localdx = XSlice[1] - XSlice[0]
		self.localdy = YSlice[1] - YSlice[0]
		self.localds = np.sqrt(self.localdx**2 + self.localdy**2)
		SSlice = [self.localds*x for x in range(200)]
		DataSlice = []
		for x,y in zip(XSlice,YSlice):
			DataSlice.append(DataFit(x,y))
		self.axes.plot(SSlice,DataSlice)
		self.canvas.draw()