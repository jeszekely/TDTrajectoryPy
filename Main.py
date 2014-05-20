#!/usr/bin/env python

import FDTD
import Inputs as IP
import Rotational_PES as RPES

import matplotlib.pyplot as plt
import matplotlib as mpl

from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

A = FDTD.Array2D("inputs.json")
# with PyCallGraph(output=GraphvizOutput()):
# 	A.loadTextFile("AgFilmX-3col.txt")

A.loadTextFile("AgFilmX-3col.txt")

molParams = IP.ImportJSON(A.params['MoleculeParameters']['ParamFile'])

#Calculate the rotational alignment and position dependent eigenenergy
alignFXN,enFXN = RPES.RotationalFxn(A.params,molParams)
A.getAlignment(alignFXN,1.0e10)
A.getRotationalEnergy(enFXN,1.0e10)
