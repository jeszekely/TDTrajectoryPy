#!/usr/bin/env/python

import numpy as np


class FDTDArray:
	#Default values that will change later
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

		self.EMF = np.genfromtxt(inputfile, usecols=(2), dtype=float).reshape(self.ny,self.nx)

	def setParameters(self,filename):
		try:
			self.resolution = alignment_params['FDTDField']['meep_res']
			self.aval = alignment_params['FDTDField']['meep_a']
		except:
			print "Error importing resolution and aval from json file"
			self.aval = 100.0
			self.resolution = 20.0