#!/usr/bin/env/python

import numpy as np
import Inputs as IP
from scipy import interpolate
import constants as c

class Array2D:
	#Default values that will change later
	def __init__(self,filename=None):
		if filename is not None:
			self.setParameters(filename)
		else:
			self.aval = 100.0 #nm
			self.res  = 20.0 #resolution

	def loadTextFile(self,filename):
		self.X = np.genfromtxt(filename, usecols=(0), dtype=int)
		self.Y = np.genfromtxt(filename, usecols=(1), dtype=int)
		#Add one to account for indexing starting at zero
		self.nx = max(self.X)+1
		self.ny = max(self.Y)+1

		#Get rid of multiples by redefining the X and Y axes
		self.X = np.linspace(0.0,float(self.nx-1),self.nx)
		self.Y = np.linspace(0.0,float(self.ny-1),self.ny)

		self.EMF = np.genfromtxt(filename, usecols=(2), dtype=float).reshape(self.ny,self.nx)

		#Interpolate the field and get the gradient
		self.E = interpolate.interp2d(self.X,self.Y,self.EMF,kind='cubic')

	def setParameters(self,filename):
		self.params = IP.ImportJSON(filename)
		try:
			self.resolution = self.params['FDTDField']['meep_res']
			self.aval       = self.params['FDTDField']['meep_a']
		except:
			print "Error importing resolution and aval from json file"
			self.aval       = 100.0
			self.resolution = 20.0

	def getGradients(self):
		self.EMFy, self.EMFx = np.gradient(self.EMF,self.aval/self.resolution/c.LEN,self.aval/self.resolution/c.LEN)
		self.Ex = interpolate.interp2d(self.X,self.Y,self.EMFx,kind='cubic')
		self.Ey = interpolate.interp2d(self.X,self.Y,self.EMFy,kind='cubic')


class RotPES(Array2D):
	def __init__(self,filename=None):
		Array2D.__init__(self,filename)