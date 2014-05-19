#!/usr/bin/env python

import FDTD
import Inputs as IP
import Rotational_PES as RPES

from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

A = FDTD.Array2D("inputs.json")
# with PyCallGraph(output=GraphvizOutput()):
# 	A.loadTextFile("AgFilmX-3col.txt")

A.loadTextFile("AgFilmX-3col.txt")

molParams = IP.ImportJSON(A.params['MoleculeParameters']['ParamFile'])
alignFXN,enFXN = RPES.RotationalFxn(A.params,molParams)