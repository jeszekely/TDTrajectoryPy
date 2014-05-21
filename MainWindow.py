import wx, os
import PropSplitter as ps
import Panel1Main as p1
import Panel2Main as p2
import ControlPanelMain as cp
import Rotational_PES as RPES
import Inputs as ip

class MainFrame(wx.Frame):
	def __init__(self, fdtdarray, parent, title):
		wx.Frame.__init__(self, parent, title=title, size = (1000,600))
		self.CreateStatusBar()
		self.fdtd = fdtdarray
		self.plotparams = ip.ImportJSON("PlottingParameters.json")
		#Setup the file menu commands
		filemenu  = wx.Menu()
		menuAbout = filemenu.Append(10, "&About", " About this program")
		menuOpen  = filemenu.Append(20, "&Open", " Open a file to plot")
		menuExit  = filemenu.Append(30, "E&xit", " Quit the program")

		#Setup the view menu commands
		viewmenu     = wx.Menu()
		menuViewPES  = viewmenu.Append(40, "1D PES", " Calculate alignment data for 1D slices Data")
		menuViewTraj = viewmenu.Append(70, "Trajectory", " New window for calculating molecule trajectories")
		menuEditJSON = viewmenu.Append(80, "Edit JSON", " Editor for the calculation parameter json file")

		#Setup the data menu commands
		datamenu      = wx.Menu()
		menuTransData = datamenu.Append(50,"Transpose", "Switch x and y columns in input file")
		menuFullCalc  = datamenu.Append(60, "Full Alignment", "Calculate Rotational PES and Alignment for full 2D Dataset")

		#Add all submenus to the menu bar
		menuBar = wx.MenuBar()
		menuBar.Append(filemenu, "&File")
		menuBar.Append(datamenu, "&Data")
		menuBar.Append(viewmenu, "Window")
		self.SetMenuBar(menuBar)

		#split the window into three panels
		self.split1 = ps.ProportionalSplitter(self,-1,0.35)
		self.split2 = ps.ProportionalSplitter(self.split1,-1,0.45)

		self.panel1 = p1.PlotPanel1(self.split1)	#Panel conatining the FDTD plot
		self.panel2 = p2.PlotPanel2(self.split2)	#Other plotting window, 1D
		self.panel3 = cp.ControlPanel(self.split2)	#Control panel

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
		self.alignFXN,self.enFXN = RPES.RotationalFxn(self.fdtd.params,self.fdtd.molParams)
		self.fdtd.getAlignment(self.alignFXN)
		self.fdtd.getRotationalEnergy(self.enFXN)

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
		dirname = os.getcwd()
		dlg = wx.FileDialog(self, "Choose a file", dirname, "", "*.*", wx.OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			self.filename = dlg.GetFilename()						#Returns just the filename
			dirname       = dlg.GetDirectory()						#Returns the rest of the filepath
			inputfile     = os.path.join(dirname, self.filename)	#Join Paths
			self.panel1.LoadData(inputfile)							#Load the data to the 2D plotting panel
			self.panel1.axes.cla()									#Clear the subplot
			self.panel1.Draw()										#Replot with new data
		dlg.Destroy()