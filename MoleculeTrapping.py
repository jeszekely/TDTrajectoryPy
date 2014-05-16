#!/usr/bin/env python
#
# Interactive python GUI to analyze 2D FDTD data, make plots, run calculations
# written by Joshua E. Szekely, January 2014
# v 1.0.0

import os, wx, re
import wx.stc as stc
import json
import subprocess
import pprint

import numpy as np
import scipy as sp
from scipy import interpolate
from scipy.integrate import odeint

#Import all necessary plotting modules
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

#Define a simple class for containing (x,y) coordinates
class XY():
	x = 0.0
	y = 0.0

#######################################################################################
#
# Begin Defining Global Variables that are constant across plotting frames and panels
#
########################################################################################

coord1 = XY()		#Coordinates for defining slice of data from the full 2D set
coord2 = XY()
linecolor = 'white'	#Color of the line on the 2D plot

inputfile = "/Users/joshuaszekely/Dropbox/Codes/MolecularAlignment/AgFilmX-3col.txt"		#default datafile containing the Electric field, in 3-col format
dirname = "/Users/joshuaszekely/Dropbox/Codes/MolecularAlignment/"							#Directory of datafile
fontsize = 12		#Font size used by default in plots

X = []				#Arrays to store axis data
Y = []
Data = []

XCOL = 1 			#Column of datafile containing the x coordinate
YCOL = 0 			#y coordinate
ZCOL = 2 			#field data, or z coordinate

Nx = 0 				#Number of datapoint in the x direction
Nt = 0 				#Number of datapoint in the y direction

XSlice = []			#x values of data slice defined by two coordinates
YSlice = []			#y values
SSlice = []			#actual spacing between data points among line
DataSlice = []		#z values
AlignmentSlice = []	#calculated alignment along slice
ERotSlice = []		#calculates rotational potential energy along slice
ERotGradX = []		#Gradient of the rotational potential energy surface, x direction
ERotGradY = []		#Gradient of the rotational potential energy surface, y direction

AllAlignmentData = None	#Variables for 2D alignement calculations, if needed
AllAlignmentFxn = None
AllERotData = None
AllERotFxn = None
DataFit = None
ERotGradXFit = None
ERotGradYFit = None

out_directory = "output_data"	#directory where alignement data will be printed
if not os.path.exists(out_directory):
    os.makedirs(out_directory)	#create if doesn't exist, otherwise alignment program will crash

JSONInput = "/Users/joshuaszekely/Dropbox/Codes/MolecularAlignment2/inputs.json"
MOLECInput = "/Users/joshuaszekely/Dropbox/Codes/MolecularAlignment2/Molecules.json"
AlignmentProgPath = "/Users/joshuaszekely/Dropbox/Codes/MolecularAlignment2/MolecAlignment"
AlignmentOutPath = "/Users/joshuaszekely/Dropbox/Codes/MolecularAlignment2/output_data/"


##################################################################
#
# Prepare the json file data for use throughout the program
#
##################################################################

def removeComments(string): 	#Removes comments denoted by '//'
    string = re.sub(re.compile("//.*" ) ,"" ,string)
    return string

jsonin = open(JSONInput,"r")
jsontemp = open("inputs_temp.json","w")
for line in jsonin:
	jsontemp.write(removeComments(line))	#Parse through the json file, removes problematic comments
jsontemp.write("\n")
jsonin.close()
jsontemp.close()
json_data = open("inputs_temp.json","r")	#Open the copy of the file sans comments
alignment_params = json.load(json_data)		#Load the json file into a dictionary
#pprint.pprint(alignment_params)					#Prints json input if desired
json_data.close()
os.system("rm inputs_temp.json")

alignment_params['FDTDField']['File'] = os.path.abspath("temp_slice.txt") 	#Alter the file used to the slice file
alignment_params['MoleculeParameters']['ParamFile'] = MOLECInput			#Specify the molecular parameter file locations

######################################
#
# Prepare the molecule json file data
#
######################################

mjsonin = open(MOLECInput,"r")
mjsontemp = open("minputs_temp.json","w")
for line in mjsonin:
	mjsontemp.write(removeComments(line))
mjsontemp.write("\n")
mjsonin.close()
mjsontemp.close()
mjson_data = open("minputs_temp.json","r")
molecule_params = json.load(mjson_data)
mjson_data.close()
os.system("rm minputs_temp.json")

#print molecule_params[alignment_params['MoleculeParameters']['Molecule']]['Mass1']

###################################
#
# Other universally useful routines
#
###################################

def WriteJSON():	#Write the json parameter file as currently defined to a file
	global alignment_params
	with open('inputs_mod.json', 'w') as outfile:
		json.dump(alignment_params, outfile, sort_keys=True, indent=4, ensure_ascii=False)

WriteJSON()

def WriteRunFile(progpath, inppath):	#Write the bash script used to execute the alignment script
	cscript = open('temp_script.sh',"w")
	cscript.write("""
#!/bin/bash
source  /opt/intel/composer_xe_2013_sp1.1.103/mkl/bin/mklvars.sh intel64
%s %s """ % (progpath, inppath))
	cscript.close()


####################################################################################
#
# Class definition to split frames into separate panels which automatically resize
# This section is taken from the internets, not my own work
#
####################################################################################

class ProportionalSplitter(wx.SplitterWindow):
	def __init__(self,parent, id = -1, proportion=0.66, size = wx.DefaultSize, **kwargs):
		wx.SplitterWindow.__init__(self,parent,id,wx.Point(0, 0),size, wx.SP_NOBORDER,**kwargs)
		self.SetMinimumPaneSize(50) #the minimum size of a pane.
		self.proportion = proportion
		if not 0 < self.proportion < 1:
			raise ValueError, "proportion value for ProportionalSplitter must be between 0 and 1."
		self.ResetSash()
		self.Bind(wx.EVT_SIZE, self.OnReSize)
		self.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGED, self.OnSashChanged, id=id)
		##hack to set sizes on first paint event
		self.Bind(wx.EVT_PAINT, self.OnPaint)
		self.firstpaint = True

	def SplitHorizontally(self, win1, win2):
		if self.GetParent() is None: return False
		return wx.SplitterWindow.SplitHorizontally(self, win1, win2,
			int(round(self.GetParent().GetSize().GetHeight() * self.proportion)))

	def SplitVertically(self, win1, win2):
		if self.GetParent() is None: return False
		return wx.SplitterWindow.SplitVertically(self, win1, win2,
			int(round(self.GetParent().GetSize().GetWidth() * self.proportion)))

	def GetExpectedSashPosition(self):
		if self.GetSplitMode() == wx.SPLIT_HORIZONTAL:
			tot = max(self.GetMinimumPaneSize(),self.GetParent().GetClientSize().height)
		else:
			tot = max(self.GetMinimumPaneSize(),self.GetParent().GetClientSize().width)
		return int(round(tot * self.proportion))

	def ResetSash(self):
		self.SetSashPosition(self.GetExpectedSashPosition())

	def OnReSize(self, event):
		"Window has been resized, so we need to adjust the sash based on self.proportion."
		self.ResetSash()
		event.Skip()

	def OnSashChanged(self, event):
		"We'll change self.proportion now based on where user dragged the sash."
		pos = float(self.GetSashPosition())
		if self.GetSplitMode() == wx.SPLIT_HORIZONTAL:
			tot = max(self.GetMinimumPaneSize(),self.GetParent().GetClientSize().height)
		else:
			tot = max(self.GetMinimumPaneSize(),self.GetParent().GetClientSize().width)
		self.proportion = pos / tot
		event.Skip()

	def OnPaint(self,event):
		if self.firstpaint:
			if self.GetSashPosition() != self.GetExpectedSashPosition():
				self.ResetSash()
			self.firstpaint = False
		event.Skip()

##################################################################
#
# Main window containing the 2D plot, 1D slice, and control panel
#
##################################################################

class MainFrame(wx.Frame):
	def __init__(self, parent, title):
		wx.Frame.__init__(self, parent, title=title, size = (1000,600))
		global inputfile, dirname, Nx, Ny
		self.CreateStatusBar()

		#Setup the file menu commands
		filemenu = wx.Menu()
		menuAbout = filemenu.Append(10, "&About", " About this program")
		menuOpen = filemenu.Append(20, "&Open", " Open a file to plot")
		menuExit = filemenu.Append(30, "E&xit", " Quit the program")

		#Setup the view menu commands
		viewmenu = wx.Menu()
		menuViewPES = viewmenu.Append(40, "1D PES", " Calculate alignment data for 1D slices Data")
		menuViewTraj = viewmenu.Append(70, "Trajectory", " New window for calculating molecule trajectories")
		menuEditJSON = viewmenu.Append(80, "Edit JSON", " Editor for the calculation parameter json file")

		#Setup the data menu commands
		datamenu = wx.Menu()
		menuTransData = datamenu.Append(50,"Transpose", "Switch x and y columns in input file")
		menuFullCalc = datamenu.Append(60, "Full Alignment", "Calculate Rotational PES and Alignment for full 2D Dataset")

		#Add all submenus to the menu bar
		menuBar = wx.MenuBar()
		menuBar.Append(filemenu, "&File")
		menuBar.Append(datamenu, "&Data")
		menuBar.Append(viewmenu, "Window")
		self.SetMenuBar(menuBar)

		#split the window into three panels
		self.split1 = ProportionalSplitter(self,-1,0.35)
		self.split2 = ProportionalSplitter(self.split1,-1,0.45)

		self.panel1 = PlotPanel1(self.split1)	#Panel conatining the FDTD plot
		self.panel2 = PlotPanel2(self.split2)	#Other plotting window, 1D
		self.panel3 = ControlPanel(self.split2)	#Control panel

		#Assign the panels their appropriate locations
		self.split1.SplitVertically(self.panel1, self.split2)
		self.split2.SplitVertically(self.panel2, self.panel3)

		#Bind the menu items to their respective functions
		self.Bind(wx.EVT_MENU, self.OnOpen, menuOpen)
		self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
		self.Bind(wx.EVT_MENU, self.OnAbout, menuAbout)
		self.Bind(wx.EVT_MENU, self.ViewPES, menuViewPES)
		self.Bind(wx.EVT_MENU, self.GetAllAlignment, menuFullCalc)
		self.Bind(wx.EVT_MENU, self.ViewTraj, menuViewTraj)
		self.Bind(wx.EVT_MENU, self.EditJSON, menuEditJSON)

		self.Bind(wx.EVT_MENU, self.panel1.Transpose, menuTransData)

	def GetAllAlignment(self,e):	#Get Alignment data for the entire 2D data range, interpolate results
		global AllAlignmentFxn, AllERotFxn, AlignmentProgPath, AllERotData, AllAlignmentData, alignment_params, inputfile,Nx,Ny
		alignment_params['FDTDField']['File'] = os.path.abspath(inputfile) 		#Change json input file to full 2D file
		WriteJSON()																#Write json
		WriteRunFile(AlignmentProgPath,os.path.abspath('inputs_mod.json'))		#write temp bash script
		os.system("chmod +x temp_script.sh")									#make executable
		subprocess.call("./temp_script.sh", shell=True)							#run script
		os.system("rm temp_script.sh")											#remove temp

		#Import the calculated data and make global interpolations
		AllAlignmentData = np.genfromtxt("output_data/AlignmentData.txt", usecols=(3), dtype=float).reshape(Ny,Nx)
		AllAlignmentFxn = interpolate.interp2d(self.panel1.X,self.panel1.Y,AllAlignmentData,kind='cubic')
		AllERotData = np.genfromtxt("output_data/AlignmentData.txt", usecols=(4), dtype=float).reshape(Ny,Nx)
		AllERotFxn = interpolate.interp2d(self.panel1.X,self.panel1.Y,AllERotData,kind='cubic')

	def ViewPES(self,e):	#Open additional window to show calculated 1D plots
		self.fr = AlignmentFrame(self, title='1D Data')
		self.fr.Show()

	def ViewTraj(self,e):	#Open trajectory calculation window
		self.fr = TrajectoryFrame(self, title='Trajectory')
		self.fr.Show()

	def EditJSON(self,e):
		self.fr = JSONFrame(self, title="JSON Editor")
		self.fr.Show()

	def OnAbout(self, e):	#Display about program message
		dlg = wx.MessageDialog(self, "Molecule trapping in FDTD potentials\n v 1.0.0, Joshua E. Szekely")
		dlg.ShowModal()
		dlg.Destroy()

	def OnExit(self, e):	#Clean up temp files prior to exiting program
		try:
			os.system("rm inputs_mod.json")
			os.system("rm temp_slice.txt")
		except:
			pass
		self.Close(True)

	def OnOpen(self, e):	#Open new file and redraw where appropriate
		global inputfile, dirname
		dlg = wx.FileDialog(self, "Choose a file", dirname, "", "*.*", wx.OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			self.filename = dlg.GetFilename()					#Returns just the filename
			dirname = dlg.GetDirectory()						#Returns the rest of the filepath
			inputfile = os.path.join(dirname, self.filename)	#Join Paths
			self.panel1.LoadData(inputfile)				#Load the data to the 2D plotting panel
			self.panel1.axes.cla()								#Clear the subplot
			self.panel1.Draw()									#Replot with new data
		dlg.Destroy()

################################################
#
# Panel containing the 2D plot of the FDTD data
#
################################################

class PlotPanel1(wx.Panel):
	def __init__(self, parent, **kwargs):
		wx.Panel.__init__(self, parent, **kwargs)
		global inputfile
		self.activecoord = 1 		#Defines coordinate 1 to be "active", coordinate 2 inactive
		self.LoadData(inputfile)
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

	def LoadData(self, inputdata):	#Get the data needed for plotting
		global XCOL,YCOL,ZCOL, inputfile, Nx, Ny, alignment_params, DataFit,X,Y, Data

		self.X = np.genfromtxt(inputfile, usecols=(XCOL), dtype=int)
		self.Y = np.genfromtxt(inputfile, usecols=(YCOL), dtype=int)
		Nx = self.maxX = max(self.X)
		Ny = self.maxY = max(self.Y)
		Nx += 1
		Ny += 1

		#Get rid of multiples by redefining the X and Y axes
		self.X = np.linspace(0.0,self.maxX,Nx)
		self.Y = np.linspace(0.0,self.maxY,Ny)

		self.data = np.genfromtxt(inputfile, usecols=(ZCOL), dtype=float).reshape(Ny,Nx)

		try:
			self.resolution = alignment_params['FDTDField']['meep_res']
			self.aval = alignment_params['FDTDField']['meep_a']
		except:
			print "Error importing resolution and aval from json file"
			self.aval = 100.0
			self.resolution = 20.0

		self.X /= self.resolution
		self.Y /= self.resolution

		self.maxX = max(self.X)	#Need redefinition since array was altered
		self.maxY = max(self.Y)

		X,Y,Data = self.X, self.Y, self.data

		#Set plotting variables
		self.xmin = 0.0
		self.xmax = self.maxX
		self.ymin = 0.0
		self.ymax = self.maxY
		self.cmin = 0.0
		self.cmax = int(max(self.data.flatten()))
		if self.cmin == self.cmax:
			self.cmax += 0.5
		self.numxtics = 3
		self.numytics = 5
		self.numctics = 5
		self.color = 'jet'

		DataFit = interpolate.interp2d(self.X,self.Y,self.data,kind='cubic')

	def Draw(self):	#Plot the 2D data
		global fontsize
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

		self.axes.set_xticklabels(self.xticks,fontweight='bold',fontsize=fontsize)
		self.axes.set_yticklabels(self.yticks,fontweight='bold',fontsize=fontsize)

		self.axes.set_xlabel("X-plane / %i nm" % int(self.aval), fontweight='bold',fontsize=fontsize)
		self.axes.set_ylabel("Y-plane / %i nm" % int(self.aval), fontweight='bold',fontsize=fontsize)

		t = self.axes.set_title(r'$\mathbf{|\vec{E}|^2/ |\vec E_0|^2}$', fontweight='bold',fontsize=fontsize+1)
		t.set_y(1.04) 							#Offset the title so it doesn't overlap with the plot
		im = self.axes.imshow(self.data,cmap=self.color,origin='lower',alpha=0.9,extent=[self.xmin,self.xmax,self.ymin,self.ymax])
		im.set_clim(self.cmin,self.cmax) 		#Set the colorbar limits
		self.cbar = self.figure.colorbar(im, ticks = self.cticks)
		self.cbar.ax.set_yticklabels(self.cticks,fontweight='bold',fontsize=fontsize)

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
		self.line_segment = self.axes.plot([Pt1.x,Pt2.x],[Pt1.y,Pt2.y],color=linecolor)
		self.GetGrandParent().panel2.Draw()				#Draw and update the 1D Panel to show the slice
		self.GetGrandParent().panel2.canvas.draw()
		self.canvas.draw()								#Update 2D plot

	def Transpose(self,e):	#Flip x and y dimensions of the data set, reloads and replots data
		global XCOL, YCOL, inputfile, Nx, Ny
		XCOL,YCOL = YCOL,XCOL
		Nx,Ny = Ny,Nx
		self.LoadData(inputfile)
		self.Draw()

#######################################
#
# Panel containing the 1D data slice
#
######################################

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

#################################################
#
# Control panel for saving and manipulating plots
#
#################################################

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

###############################################################
#
# Window to plot the FDTD data slice, cos^2, and rotational PES
#
###############################################################

class AlignmentFrame(wx.Frame):
	def __init__(self, parent, title):
		wx.Frame.__init__(self, parent, title=title, size = (600,1000))
		global inputfile, dirname, XSlice, YSlice, DataSlice, SSlice
		self.CreateStatusBar()

		#Setup the file menu options, independent from the main window
		filemenu = wx.Menu()
		menuExit = filemenu.Append(100, "E&xit", " Quit the program")

		calcmenu = wx.Menu()
		menuCalc = calcmenu.Append(110, "&Calculate", " Calculate the Alignment")

		menuBar = wx.MenuBar()
		menuBar.Append(filemenu, "&File")
		menuBar.Append(calcmenu, "&Calculations")
		self.SetMenuBar(menuBar)

		self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
		self.Bind(wx.EVT_MENU, self.GetAlignment, menuCalc)

		self.split1 = ProportionalSplitter(self,-1,0.33)
		self.split2 = ProportionalSplitter(self.split1,-1,0.33)

		self.panel1 = PES1D(self.split1)	#Slice of data, same as in main window
		self.panel2 = PESCOS(self.split2)	#Rotational alignment of slice
		self.panel3 = PESEROT(self.split2)	#Rotational PES of slice

		self.split1.SplitHorizontally(self.panel1, self.split2)
		self.split2.SplitHorizontally(self.panel2, self.panel3)

		#Write a temp file that will be used to store the slice data for use by the alignment program
		self.tempfile = open("temp_slice.txt","w")
		for i in range(len(SSlice)):
			self.tempfile.write("0\t%i\t%f\n" % (i,DataSlice[i]))
		self.tempfile.close()

	def GetAlignment(self,e):
		global AlignmentSlice, ERotSlice, AlignmentProgPath, alignment_params
		alignment_params['FDTDField']['File'] = os.path.abspath("temp_slice.txt")
		WriteJSON()
		WriteRunFile(AlignmentProgPath, os.path.abspath('inputs_mod.json'))
		os.system("chmod +x temp_script.sh")
		subprocess.call("./temp_script.sh", shell=True)
		os.system("rm temp_script.sh")

		AlignmentSlice = np.genfromtxt("output_data/AlignmentData.txt",usecols=(3))
		ERotSlice = np.genfromtxt("output_data/AlignmentData.txt",usecols=(4))

		self.panel1.Draw()
		self.panel2.Draw()
		self.panel3.Draw()

		self.split1.ResetSash()
		self.split1.ResetSash()

	def OnExit(self, e):
		os.system('rm temp_slice.txt')
		self.Close(True)

###############################################################
#
# Top panel of alignment window, 1D slice of data
#
###############################################################

class PES1D(wx.Panel):
	def __init__(self, parent, **kwargs):
		wx.Panel.__init__(self, parent, **kwargs)
		global SSlice, DataSlice
		self.Draw()

	def Draw(self):
		global SSlice, DataSlice, alignment_params
		self.figure = plt.Figure()							#Figure needed to plot cbar
		self.axes = self.figure.add_subplot(111)			#For the heatmap
		self.canvas = FigureCanvas(self, -1, self.figure)	#No canvas, no plot
		self.axes.cla()
		try:	#The slice might not be definied initially so this catches the error
			self.axes.set_xlabel("Position / %i nm" % int(alignment_params['FDTDField']['meep_a']), fontweight='bold',fontsize=fontsize)
			self.axes.set_ylabel(r'$\mathbf{|\vec{E}|^2/ |\vec E_0|^2}$', fontweight='bold',fontsize=fontsize)
			t = self.axes.set_title("Field Enhancement", fontweight='bold',fontsize=fontsize+1)
			self.axes.plot(SSlice,DataSlice)
		except:
			pass
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		self.sizer.Add(self.canvas, 1, wx.ALL| wx.EXPAND)
		self.SetSizer(self.sizer)
		self.Fit()
		self.canvas.draw()

####################################################
#
# Middle panel of alignment window, <cos^2 theta>
#
####################################################

class PESCOS(wx.Panel):
	def __init__(self, parent, **kwargs):
		wx.Panel.__init__(self, parent, **kwargs)
		global SSlice, AlignmentSlice
		self.Draw()

	def Draw(self):
		global SSlice, AlignmentSlice, alignment_params
		self.figure = plt.Figure()							#Figure needed to plot cbar
		self.axes = self.figure.add_subplot(111)			#For the heatmap
		self.canvas = FigureCanvas(self, -1, self.figure)	#No canvas, no plot
		self.axes.cla()
		try:	#The slice might not be definied initially so this catches the error
			self.axes.set_xlabel("Position / %i nm" % int(alignment_params['FDTDField']['meep_a']), fontweight='bold',fontsize=fontsize)
			self.axes.set_ylabel(r'$\langle \cos^2 \theta \rangle$', fontweight='bold',fontsize=fontsize)
			t = self.axes.set_title('Alignment in Field', fontweight='bold',fontsize=fontsize+1)
			self.axes.plot(SSlice,AlignmentSlice)
		except:
			pass
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		self.sizer.Add(self.canvas, 1, wx.ALL| wx.EXPAND)
		self.SetSizer(self.sizer)
		self.Fit()
		self.canvas.draw()

####################################################
#
# Bottom panel of alignment window, rotational PES
#
####################################################

class PESEROT(wx.Panel):
	def __init__(self, parent, **kwargs):
		wx.Panel.__init__(self, parent, **kwargs)
		global SSlice, ERotSlice
		self.Draw()

	def Draw(self):
		global SSlice, ERotSlice, alignment_params
		self.figure = plt.Figure()							#Figure needed to plot cbar
		self.axes = self.figure.add_subplot(111)			#For the heatmap
		self.canvas = FigureCanvas(self, -1, self.figure)	#No canvas, no plot
		self.axes.cla()
		try:	#The slice might not be definied initially so this catches the error
			self.axes.set_xlabel("Position / %i nm" % int(alignment_params['FDTDField']['meep_a']), fontweight='bold',fontsize=fontsize)
			self.axes.set_ylabel("Energy / au", fontweight='bold',fontsize=fontsize)
			t = self.axes.set_title('Rotational Potential Energy', fontweight='bold',fontsize=fontsize+1)
			self.axes.plot(SSlice,ERotSlice)
		except:
			pass
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		self.sizer.Add(self.canvas, 1, wx.ALL| wx.EXPAND)
		self.SetSizer(self.sizer)
		self.Fit()
		self.canvas.draw()

####################################################
#
# New Window to to create trajectory images
#
####################################################

class TrajectoryFrame(wx.Frame):
	def __init__(self, parent, title):
		wx.Frame.__init__(self, parent, title=title, size = (600,500))
		global inputfile, dirname, DataFit, AllERotData, AllERotFxn, X, Y, alignment_params, ERotGradX, ERotGradY, ERotGradYFit, ERotGradXFit, molecule_params
		self.CreateStatusBar()

		#Setup the file menu options, independent from the main window
		filemenu = wx.Menu()
		menuExit = filemenu.Append(200, "E&xit", " Quit the program")

		calcmenu = wx.Menu()
		menuCalc = calcmenu.Append(210, "&Calculate", " Calculate the Trajectory")

		menuBar = wx.MenuBar()
		menuBar.Append(filemenu, "&File")
		menuBar.Append(calcmenu, "&Calculations")
		self.SetMenuBar(menuBar)

		self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
		self.Bind(wx.EVT_MENU, self.GetTrajectory, menuCalc)

		self.split1 = ProportionalSplitter(self,-1,0.65)

		self.panel1 = TrajectoryDispPanel(self.split1)	#2D plot of data, with trajectory
		self.panel2 = TrajectoryCtrlPanel(self.split1)	#Trajectory variables and controls

		self.split1.SplitVertically(self.panel1, self.panel2)

		#multiply dist by this num to go from plotted units to au
		self.distance2au = alignment_params['FDTDField']['meep_a']/0.05291772

		#Specify Trajectory initial conditions, provided in meep units
		self.x0 = 5.0
		self.y0 = 30.0
		self.px0 = 0.0
		self.py0 = 0.0

		self.maxX = max(X)
		self.maxY = max(Y)

		#Get the partial derivatives of the potential energy surface, interpolate
		ERotGradY, ERotGradX = np.gradient(AllERotData,alignment_params['FDTDField']['meep_a']/alignment_params['FDTDField']['meep_res']/0.05291772,alignment_params['FDTDField']['meep_a']/alignment_params['FDTDField']['meep_res']/0.05291772)
		ERotGradXFit = interpolate.interp2d(X,Y,ERotGradX,kind='cubic')
		ERotGradYFit = interpolate.interp2d(X,Y,ERotGradY,kind='cubic')

		#Get the particle mass, convert to au
		self.mass =  molecule_params[alignment_params['MoleculeParameters']['Molecule']]['Mass1'] + molecule_params[alignment_params['MoleculeParameters']['Molecule']]['Mass2']
		self.mass *= 1836.0

	def OnExit(self, e):
		self.Close(True)

	#Brief routine needed by the differential equation solver
	#Assumes that y array is in atomic units
	def ClassicalHamiltonian(self,y,t):
		#Factors in the periodicity of the potential
		if y[0]/self.distance2au > self.maxX:
			y[0] -= self.maxX*self.distance2au
			# y[4] += 1
		elif y[0]/self.distance2au < 0.0:
			y[0] += self.maxX*self.distance2au
			# y[4] -= 1

		return [y[2]/self.mass, y[3]/self.mass, -1.0*ERotGradXFit(y[0]/self.distance2au,y[1]/self.distance2au), -1.0*ERotGradYFit(y[0]/self.distance2au,y[1]/self.distance2au)]

	def GetTrajectory(self,e):
		#Convert all to au before propagation
		self.x0 *= self.distance2au
		self.y0 *= self.distance2au
		self.ICs = [self.x0, self.y0, self.px0, self.py0]
		self.t = np.linspace(0,8.0e8,801)/(2.41888432e-2) #go to 400 fs, convert to au as expected

		self.traj = odeint(self.ClassicalHamiltonian, self.ICs, self.t)

		self.X = [y[0]/self.distance2au for y in self.traj]
		self.Y = [y[1]/self.distance2au for y in self.traj]
		for x,y in zip(self.X,self.Y):
			print x,y
		self.panel1.AddLine(self.X,self.Y,2.41888432e-2*self.t)

####################################
#
# Frame for 2D Plot with Trajectory
#
####################################

class TrajectoryDispPanel(wx.Panel):
	def __init__(self, parent, **kwargs):
		wx.Panel.__init__(self, parent, **kwargs)
		global inputfile,XCOL,YCOL,ZCOL,X,Y,Data,AllERotData,ERotGradX,ERotGradY
		self.LoadParameters()
		self.Draw = self.DrawERot
		self.Draw()
		self.plotmax = 1.0
		self.plotmin = 0.0

	def OnPress(self,e):	#On mouse click, set the active coordinate
		try:
			print 'Press:   button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(e.button, e.x, e.y, e.xdata, e.ydata)
			self.GetGrandParent().x0, self.GetGrandParent().y0 = e.xdata,e.ydata
		except:
			pass

	def LoadParameters(self):	#Set up the plot like the main window
		try:
			self.resolution = alignment_params['FDTDField']['meep_res']
			self.aval = alignment_params['FDTDField']['meep_a']
		except:
			print "Error importing resolution and aval from json file"
			self.aval = 100.0
			self.resolution = 20.0

		self.maxX = max(X)	#Need redefinition since array was altered
		self.maxY = max(Y)

		#Set plotting variables
		self.xmin = 0.0
		self.xmax = self.maxX
		self.ymin = 0.0
		self.ymax = self.maxY
		self.cmin = 0.0
		self.cmax = int(max(AllERotData.flatten()))
		if self.cmin == self.cmax:
			self.cmax += 0.5
		self.numxtics = 3
		self.numytics = 5
		self.numctics = 5
		self.color = 'jet'

	def PlotERot(self,e):
		self.Draw = self.DrawERot
		self.plotmax = max(AllERotData.flatten())
		self.plotmin = min(AllERotData.flatten())
		self.Draw()

	def PlotERotGradX(self,e):
		self.Draw = self.DrawERotGradX
		self.plotmax = max(ERotGradX.flatten())
		self.plotmin = min(ERotGradX.flatten())
		self.Draw()

	def PlotERotGradY(self,e):
		self.Draw = self.DrawERotGradY
		self.plotmax = max(ERotGradY.flatten())
		self.plotmin = min(ERotGradY.flatten())
 		self.Draw()

	def DrawERot(self):	#Plot the 2D data
		global fontsize
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

		self.axes.set_xticklabels(self.xticks,fontweight='bold',fontsize=fontsize)
		self.axes.set_yticklabels(self.yticks,fontweight='bold',fontsize=fontsize)

		self.axes.set_xlabel("X-plane / %i nm" % int(self.aval), fontweight='bold',fontsize=fontsize)
		self.axes.set_ylabel("Y-plane / %i nm" % int(self.aval), fontweight='bold',fontsize=fontsize)

		t = self.axes.set_title(r'$\mathbf{|\vec{E}|^2/ |\vec E_0|^2}$', fontweight='bold',fontsize=fontsize+1)
		t.set_y(1.04) 							#Offset the title so it doesn't overlap with the plot
		im = self.axes.imshow(AllERotData,cmap=self.color,origin='lower',alpha=0.9,extent=[self.xmin,self.xmax,self.ymin,self.ymax])
		im.set_clim(self.cmin,self.cmax) 		#Set the colorbar limits
		self.cbar = self.figure.colorbar(im)#, ticks = self.cticks)
		#self.cbar.ax.set_yticklabels(self.cticks,fontweight='bold',fontsize=fontsize)

		self.canvas.draw() #Refreshed the canvas so that plot appears

		#Need to reset these to get the plot to diplay and resize properly
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		self.sizer.Add(self.canvas, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER)
		self.SetSizer(self.sizer)
		self.Fit()
		self.GetParent().ResetSash()
		self.cidpress = self.canvas.mpl_connect('button_press_event',self.OnPress)

	def DrawERotGradX(self):	#Plot the 2D data
		global fontsize
		try:	#If already plotted, clears to prevent unwanted artifacts
			self.figure.clf()
		except:
			pass
		self.figure = plt.Figure()							#Figure needed to plot cbar
		self.axes = self.figure.add_subplot(111)			#For the heatmap
		self.canvas = FigureCanvas(self, -1, self.figure)	#No canvas, no plot
		self.axes.cla()										#Clear if anything was left. I don't know why I need this, but it improves functionality
		self.xticks = np.linspace(self.xmin, self.xmax, self.numxtics)
		self.yticks = np.linspace(self.ymin, self.ymax, self.numytics)
		self.cticks = np.linspace(self.cmin,self.cmax,self.numctics)
		self.axes.set_xlim([self.xmin, self.xmax])
		self.axes.set_ylim([self.ymin, self.ymax])
		self.axes.xaxis.set_ticks(self.xticks)
		self.axes.yaxis.set_ticks(self.yticks)
		self.axes.tick_params(axis='x',which='both',bottom='on',top='off',direction='out',width=1.5)
		self.axes.tick_params(axis='y',which='both',left='on',right='off',direction='out',width=1.5)
		self.axes.set_xticklabels(self.xticks,fontweight='bold',fontsize=fontsize)
		self.axes.set_yticklabels(self.yticks,fontweight='bold',fontsize=fontsize)
		self.axes.set_xlabel("X-plane / %i nm" % int(self.aval), fontweight='bold',fontsize=fontsize)
		self.axes.set_ylabel("Y-plane / %i nm" % int(self.aval), fontweight='bold',fontsize=fontsize)
		t = self.axes.set_title(r'$\mathbf{|\vec{E}|^2/ |\vec E_0|^2}$', fontweight='bold',fontsize=fontsize+1)
		t.set_y(1.04) 							#Offset the title so it doesn't overlap with the plot
		im = self.axes.imshow(ERotGradX,cmap=self.color,origin='lower',alpha=0.9,extent=[self.xmin,self.xmax,self.ymin,self.ymax])
		im.set_clim(self.cmin,self.cmax) 		#Set the colorbar limits
		self.cbar = self.figure.colorbar(im)
		self.canvas.draw() #Refreshed the canvas so that plot appears
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		self.sizer.Add(self.canvas, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER)
		self.SetSizer(self.sizer)
		self.Fit()
		self.GetParent().ResetSash()
		self.cidpress = self.canvas.mpl_connect('button_press_event',self.OnPress)

	def DrawERotGradY(self):	#Plot the 2D data
		global fontsize
		try:	#If already plotted, clears to prevent unwanted artifacts
			self.figure.clf()
		except:
			pass
		self.figure = plt.Figure()							#Figure needed to plot cbar
		self.axes = self.figure.add_subplot(111)			#For the heatmap
		self.canvas = FigureCanvas(self, -1, self.figure)	#No canvas, no plot
		self.axes.cla()										#Clear if anything was left. I don't know why I need this, but it improves functionality
		self.xticks = np.linspace(self.xmin, self.xmax, self.numxtics)
		self.yticks = np.linspace(self.ymin, self.ymax, self.numytics)
		self.cticks = np.linspace(self.cmin,self.cmax,self.numctics)
		self.axes.set_xlim([self.xmin, self.xmax])
		self.axes.set_ylim([self.ymin, self.ymax])
		self.axes.xaxis.set_ticks(self.xticks)
		self.axes.yaxis.set_ticks(self.yticks)
		self.axes.tick_params(axis='x',which='both',bottom='on',top='off',direction='out',width=1.5)
		self.axes.tick_params(axis='y',which='both',left='on',right='off',direction='out',width=1.5)
		self.axes.set_xticklabels(self.xticks,fontweight='bold',fontsize=fontsize)
		self.axes.set_yticklabels(self.yticks,fontweight='bold',fontsize=fontsize)
		self.axes.set_xlabel("X-plane / %i nm" % int(self.aval), fontweight='bold',fontsize=fontsize)
		self.axes.set_ylabel("Y-plane / %i nm" % int(self.aval), fontweight='bold',fontsize=fontsize)
		t = self.axes.set_title(r'$\mathbf{|\vec{E}|^2/ |\vec E_0|^2}$', fontweight='bold',fontsize=fontsize+1)
		t.set_y(1.04) 							#Offset the title so it doesn't overlap with the plot
		im = self.axes.imshow(ERotGradY,cmap=self.color,origin='lower',alpha=0.9,extent=[self.xmin,self.xmax,self.ymin,self.ymax])
		im.set_clim(self.cmin,self.cmax) 		#Set the colorbar limits
		self.cbar = self.figure.colorbar(im)
		self.canvas.draw() #Refreshed the canvas so that plot appears
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		self.sizer.Add(self.canvas, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER)
		self.SetSizer(self.sizer)
		self.Fit()
		self.GetParent().ResetSash()
		self.cidpress = self.canvas.mpl_connect('button_press_event',self.OnPress)

	def AddLine(self, X, Y, Z):	#Add line to the 2D plot to show where the slice is taken
		# global linecolor
		# try:
		# 	self.line_segment.pop(0).remove()
		# except:
		# 	pass
		points = np.array([X,Y]).T.reshape(-1,1,2)
		segments = np.concatenate([points[:-1], points[1:]], axis=1)
		lc = LineCollection(segments, cmap=plt.get_cmap('Spectral'))
		lc.set_array(Z)
		lc.set_linewidth(2)
		self.axes.add_collection(lc)

		self.canvas.draw()								#Update 2D plot


#################################################
#
# Control panel for trajectory inputs
#
#################################################

class TrajectoryCtrlPanel(wx.Panel):
	def __init__(self, parent, **kwargs):
		wx.Panel.__init__(self, parent, **kwargs)
		self.grid = wx.GridBagSizer(hgap=5, vgap=5)	#Helps automatically place the controls

		#2D Save Plot Button
		self.save2Dbutton = wx.Button(self, label="Save PES")
		self.grid.Add(self.save2Dbutton, pos=(2,0), flag=wx.ALL)

		#Get Trajectory
		self.trajbutton = wx.Button(self, label="Trajectory")
		self.grid.Add(self.trajbutton, pos=(2,1), flag=wx.ALL)

		#Plot Choice Buttons
		self.RotPESbutton = wx.Button(self, label="Rot. PES")
		self.grid.Add(self.RotPESbutton, pos=(3,0), flag=wx.ALL)
		self.GradXbutton = wx.Button(self, label="X Grad.")
		self.grid.Add(self.GradXbutton, pos=(3,1), flag=wx.ALL)
		self.GradYbutton = wx.Button(self, label="Y Grad.")
		self.grid.Add(self.GradYbutton, pos=(3,2), flag=wx.ALL)

		# Sliders to set cbrange of 2D plot
		self.cminText = wx.StaticText(self, label="Color Bar Minimum:")
		self.grid.Add(self.cminText, pos=(4,0), flag=wx.ALL, span=(1,2))
		self.cmaxText = wx.StaticText(self, label="Color Bar Maximum:")
		self.grid.Add(self.cmaxText, pos=(6,0), flag=wx.ALL)
		#self.cmin = self.GetGrandParent().GetParent().panel1.cmin
		#self.cmax = self.GetGrandParent().GetParent().panel1.cmax
		self.cbminSlider = wx.Slider(self, -1, 0, 0, 100, size = (100,-1), style=wx.HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
		self.cbminSlider.SetTickFreq(5,1)
		self.cbmaxSlider = wx.Slider(self, -1, 100, 0, 100, size = (100,-1), style=wx.HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
		self.cbmaxSlider.SetTickFreq(5,1)
		self.SetSizer(self.grid)
		self.grid.Add(self.cbminSlider, pos=(5,0), span=(1,2), flag=wx.ALL)
		self.grid.Add(self.cbmaxSlider, pos=(7,0), span=(1,2), flag=wx.ALL)

		# Bind the previously defined controls to the proper functions
		self.Bind(wx.EVT_BUTTON, self.Save2DPlot, self.save2Dbutton)
		self.Bind(wx.EVT_SLIDER, self.CBRUpdate)
		self.Bind(wx.EVT_BUTTON, self.GetGrandParent().GetTrajectory, self.trajbutton)

		self.Bind(wx.EVT_BUTTON, self.GetGrandParent().panel1.PlotERot, self.RotPESbutton)
		self.Bind(wx.EVT_BUTTON, self.GetGrandParent().panel1.PlotERotGradX, self.GradXbutton)
		self.Bind(wx.EVT_BUTTON, self.GetGrandParent().panel1.PlotERotGradY, self.GradYbutton)

	def CBRUpdate(self,e):
		self.crange = float(self.GetGrandParent().panel1.plotmax) - float(self.GetGrandParent().panel1.plotmin)
		self.GetGrandParent().panel1.cmax = float(self.GetGrandParent().panel1.plotmin) + self.crange*(float(self.cbmaxSlider.GetValue())/100.0)
		self.GetGrandParent().panel1.cmin = float(self.GetGrandParent().panel1.plotmin) + self.crange*(float(self.cbminSlider.GetValue())/100.0)
		#self.GetGrandParent().panel1.cmin = float(self.cbminSlider.GetValue())
		#self.GetGrandParent().panel1.cmax = float(self.cbmaxSlider.GetValue())
		self.GetGrandParent().panel1.Draw()

	def Save2DPlot(self, e): #Save the 2D PES plot to a file
		dlg = wx.FileDialog(self, "Save 2D Plot", dirname, "", "*.*", wx.SAVE)
		if dlg.ShowModal() == wx.ID_OK:
			self.filename = dlg.GetFilename()
			self.dirname = dlg.GetDirectory()
			self.plot2Dfile = os.path.join(self.dirname, self.filename)
			self.GetGrandParent().panel1.figure.savefig(self.plot2Dfile,dpi=300)
		dlg.Destroy()

################################################################################################
#
# Next class calls on pprint function, but this function is needed to modify functionality
# Prints 'u' next to unicode objects by default when converting from json data structure to string
# Define new print funciton that supresses this output
#
################################################################################################

def my_safe_repr(object, context, maxlevels, level):
	typ = pprint._type(object)
	if typ is unicode:
		object = str(object)
	return pprint._safe_repr(object, context, maxlevels, level)

#################################################
#
# Window for modifying the c code calc parameters
#
#################################################

class JSONFrame(wx.Frame):
	def __init__(self, parent, title):
		wx.Frame.__init__(self, parent, title=title, size=(600,400))
		global JSONInput, alignment_params
		self.control = stc.StyledTextCtrl(self) #TE_MULTILINE)
		self.control.SetLexer(stc.STC_LEX_AUTOMATIC)
		self.CreateStatusBar()

		# Setting up the menu.
		filemenu= wx.Menu()

		menuExit = filemenu.Append(310,"E&xit", " Exit the window")

		# Creating the menubar.
		menuBar = wx.MenuBar()
		menuBar.Append(filemenu,"&File") 	# Adding the "filemenu" to the MenuBar
		self.SetMenuBar(menuBar)  			# Adding the MenuBar to the Frame content.

		# Set events.
		self.Bind(wx.EVT_MENU, self.OnExit, menuExit)

		#Load the JSON File into the editor
		printer = pprint.PrettyPrinter()
		printer.format = my_safe_repr
		self.control.SetValue(printer.pformat(alignment_params))

		self.Show(True)

 	def OnExit(self,e):
 		global alignment_params
 		try:
	 		alignment_params = json.loads(unicode(self.control.GetValue().replace('\'',"\"")))		#Load new parameters to json object
			WriteJSON()
		except:
			pass
		self.Close()  # Close the frame


##############################
#
# Finally, run the application
#
##############################

if __name__ == "__main__":
	app = wx.PySimpleApp()
	fr = MainFrame(None, title='Molecular Alignment in FDTD Fields')
	fr.Show()
	app.MainLoop()
