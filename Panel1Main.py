
################################################
#
# Panel containing the 2D plot of the FDTD data
#
################################################
import wx
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import numpy as np

class PlotPanel1(wx.Panel):
	def __init__(self, parent, **kwargs):
		wx.Panel.__init__(self, parent, **kwargs)
		self.activecoord = 1 		#Defines coordinate 1 to be "active", coordinate 2 inactive
		self.LoadData()
		self.Draw()

	def OnPress(self,e):	#On mouse click, set the active coordinate
		try:
			self.SetActive(e.xdata,e.ydata)
			print 'Press:   button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(e.button, e.x, e.y, e.xdata, e.ydata)
		except:
			pass

	def OnRelease(self,e):	#On mouse release, set the inactive coordinate
		try:
			self.SetInactive(e.xdata,e.ydata)
			self.AddLine(coord1,coord2)
			print 'Release: button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(e.button, e.x, e.y, e.xdata, e.ydata)
		except:
			pass

	def SetActive(self,x,y):	#Better name would be "primary"
		if self.activecoord == 1:
			coord1.x = x
			coord1.y = y
		else:
			coord2.x = x
			coord2.y = y

	def SetInactive(self,x,y):	# "secondary"
		if self.activecoord == 2:
			coord1.x = x
			coord1.y = y
		else:
			coord2.x = x
			coord2.y = y

	def SwitchActiveCoordinate(self): #Flip the active and inactive coordinate
		if self.activecoord == 1:
			self.activecoord = 2
		else:
			self.activecoord = 1

	def LoadData(self):	#Get the data needed for plotting
		self.X    = self.GetGrandParent().fdtd.X
		self.Y    = self.GetGrandParent().fdtd.Y
		self.maxX = max(self.X)
		self.maxY = max(self.Y)
		self.Nx   = len(self.X)
		self.Ny   = len(self.Y)

		#Set plotting variables
		self.xmin = 0.0
		self.xmax = self.maxX
		self.ymin = 0.0
		self.ymax = self.maxY
		self.cmin = 0.0
		self.cmax = int(max(self.GetGrandParent().fdtd.EMF.flatten()))
		if self.cmin == self.cmax:
			self.cmax += 0.5
		self.numxtics = 3
		self.numytics = 5
		self.numctics = 5
		self.color    = 'jet'

#		DataFit = interpolate.interp2d(self.X,self.Y,self.GetGrandParent().fdtd,kind='cubic')

	def Draw(self):	#Plot the 2D data
		try:	#If already plotted, clears to prevent unwanted artifacts
			self.figure.clf()
		except:
			pass

		self.figure = plt.Figure()							#Figure needed to plot cbar
		self.axes = self.figure.add_subplot(111)			#For the heatmap
		self.canvas = FigureCanvas(self, -1, self.figure)	#No canvas, no plot
		self.axes.cla()										#Clear if anything was left. I don't know why I need this, but it improves functionality

		#Define all parameters for the tick locations, axis limits, labels etc...
		self.xticks = np.linspace(self.xmin, self.xmax, self.numxtics)
		self.yticks = np.linspace(self.ymin, self.ymax, self.numytics)
		self.cticks = np.linspace(self.cmin,self.cmax,self.numctics)

		self.axes.set_xlim([self.xmin, self.xmax])
		self.axes.set_ylim([self.ymin, self.ymax])

		self.axes.xaxis.set_ticks(self.xticks)
		self.axes.yaxis.set_ticks(self.yticks)

		self.axes.tick_params(axis='x',which='both',bottom='on',top='off',direction='out',width=1.5)
		self.axes.tick_params(axis='y',which='both',left='on',right='off',direction='out',width=1.5)

		self.axes.set_xticklabels(self.xticks,fontweight='bold',fontsize=self.GetGrandParent().plotparams['Font']['fontsize'])
		self.axes.set_yticklabels(self.yticks,fontweight='bold',fontsize=self.GetGrandParent().plotparams['Font']['fontsize'])

		self.axes.set_xlabel("X-plane / %i nm" % int(self.GetGrandParent().fdtd.aval), fontweight='bold',fontsize=self.GetGrandParent().plotparams['Font']['fontsize'])
		self.axes.set_ylabel("Y-plane / %i nm" % int(self.GetGrandParent().fdtd.aval), fontweight='bold',fontsize=self.GetGrandParent().plotparams['Font']['fontsize'])

		t = self.axes.set_title(r'$\mathbf{|\vec{E}|^2/ |\vec E_0|^2}$', fontweight='bold',fontsize=self.GetGrandParent().plotparams['Font']['fontsize']+1)
		t.set_y(1.04) 							#Offset the title so it doesn't overlap with the plot
		im = self.axes.imshow(self.GetGrandParent().fdtd.EMF,cmap=self.color,origin='lower',alpha=0.9,extent=[self.xmin,self.xmax,self.ymin,self.ymax])
		im.set_clim(self.cmin,self.cmax) 		#Set the colorbar limits
		self.cbar = self.figure.colorbar(im, ticks = self.cticks)
		self.cbar.ax.set_yticklabels(self.cticks,fontweight='bold',fontsize=self.GetGrandParent().plotparams['Font']['fontsize'])

		self.canvas.draw() #Refreshed the canvas so that plot appears

		#Need to reset these to get the plot to diplay and resize properly
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		self.sizer.Add(self.canvas, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER)
		self.SetSizer(self.sizer)
		self.Fit()
		self.GetParent().ResetSash()
		self.cidpress = self.canvas.mpl_connect('button_press_event',self.OnPress)
		self.cidrelease = self.canvas.mpl_connect('button_release_event',self.OnRelease)

	def AddLine(self, Pt1, Pt2):	#Add line to the 2D plot to show where the slice is taken
		global linecolor
		try:
			self.line_segment.pop(0).remove()
		except:
			pass
		self.line_segment = self.axes.plot([Pt1.x,Pt2.x],[Pt1.y,Pt2.y],color=self.GetGrandParent().plotparams['LineSlice']['color'])
		self.GetGrandParent().panel2.Draw()				#Draw and update the 1D Panel to show the slice
		self.GetGrandParent().panel2.canvas.draw()
		self.canvas.draw()								#Update 2D plot

	def Transpose(self,e):	#Flip x and y dimensions of the data set, reloads and replots data
		self.GetGrandParent().fdtd.transpose()
		self.LoadData()
		self.Draw()
