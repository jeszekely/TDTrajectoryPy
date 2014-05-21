#!/usr/bin/env python

import FDTD
import Inputs as IP
import Rotational_PES as RPES

import matplotlib.pyplot as plt
import matplotlib as mpl

from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

import wx
import MainWindow as mw

A = FDTD.Array2D("inputs.json")
# with PyCallGraph(output=GraphvizOutput()):
# 	A.loadTextFile("AgFilmX-3col.txt")
A.loadTextFile("AgFilmX-3col.txt")

#Calculate the rotational alignment and position dependent eigenenergy
# alignFXN,enFXN = RPES.RotationalFxn(A.params,molParams)
# A.getAlignment(alignFXN,1.0e10)
# A.getRotationalEnergy(enFXN,1.0e10)




if __name__ == "__main__":
	app = wx.PySimpleApp()
	fr = mw.MainFrame(A,None,title='Molecular Alignment in FDTD Fields')
	fr.Show()
	app.MainLoop()
