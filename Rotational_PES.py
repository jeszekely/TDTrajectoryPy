#!/usr/bin/env python

import numpy as np
import constants as c
import math
import sympy.physics.quantum.cg as WF
from scipy import interpolate

def kron_delta(A,B):
	if A == B:
		return 1.0
	else:
		return 0.0

def FMIME(J,K,M,Q,S,j,k,m):
	c1 = pow(-1.0,k+m)
	c2 = math.sqrt((2.0*J+1.0)*(2.0*j+1.0))
	J1 = WF.Wigner3j(J,M,2,Q,j,-1*m).doit()
	J2 = WF.Wigner3j(J,K,2,S,j,-1*k).doit()
	return c1*c2*J1*J2

def RotationalFxn(compparams,molparams):
	#Define important constants
	amass1 = molparams[compparams['MoleculeParameters']['Molecule']]['Mass1']
	amass2 = molparams[compparams['MoleculeParameters']['Molecule']]['Mass2']
	bond   = molparams[compparams['MoleculeParameters']['Molecule']]['BondLength']
	amass1 *= 1823.0;
	amass2 *= 1823.0;
	bond   /= 0.528;
	aperp  = molparams[compparams['MoleculeParameters']['Molecule']]['Polarizability']['XX']
	apara  = molparams[compparams['MoleculeParameters']['Molecule']]['Polarizability']['ZZ']
	Be     = pow(c.HBAR,2)*0.5*(1.0/pow(bond,2))*(amass1+amass2)/(amass2*amass1)

	intensities = []
	alignment   = []
	energies    = []

	power_factor = math.sqrt(math.sqrt(math.sqrt(10)))
	power_factor = 10.0
	inten = 1.0e3
	while (inten <= 1.0e16):
		intensities.append(inten/c.LASERINTEN)
		inten *= power_factor
	if compparams['CalculationParameters']['UseM'] == 0 and compparams['CalculationParameters']['UseOdd'] == 0:
		factor = 2.0
	else:
		factor = 1.0

	if compparams['CalculationParameters']['UseM'] == 0:
		NEq = compparams['CalculationParameters']['JStates']+1
	else:
		NEq = pow(compparams['CalculationParameters']['JStates']+1,2)

	#Rotational Basis
	RB = np.zeros((NEq,3))
	if compparams['CalculationParameters']['UseM'] == 0:
		kk = 0
		for ii in range(NEq):
			RB[kk][0] = factor*ii
			RB[kk][1] = 0
			RB[kk][2] = 0
			kk        += 1
	else:
		kk = 0
		for ii in range(0,compparams['CalculationParameters']['JStates']+1):
			for jj in range(-1*ii,ii+1):
				RB[kk][0] = ii
				RB[kk][1] = 0
				RB[kk][2] = jj
				kk        += 1

	#Define cosine squared matrix
	cossq = np.zeros((NEq,NEq))
	for ii in range(0,NEq):
		for jj in range(0,NEq):
			cossq[ii][jj] = (1.0/3.0)*kron_delta(RB[ii][0],RB[jj][0])*kron_delta(RB[ii][2],RB[jj][2]) + (2.0/3.0)*FMIME(RB[ii][0],RB[ii][1],RB[ii][2],0,0,RB[jj][0],RB[jj][1],RB[jj][2]);

	Htotal = np.zeros((NEq,NEq))
	for efield in intensities:
		for ii in range(0,NEq):
			for jj in range(0,NEq):
				E_JK           =  kron_delta(RB[ii][2],RB[jj][2])*kron_delta(RB[ii][0],RB[jj][0])*(Be*(RB[ii][0])*(RB[ii][0]+1) - 0.25*efield*aperp) #Be(J(J+1) - (1/4)|E|^2*a_perp
				coeff          = (-1.0/6.0)*abs(apara-aperp)*efield  #EField already squared here
				coupling       = FMIME(RB[ii][0],RB[ii][1],RB[ii][2],0,0,RB[jj][0],RB[jj][1],RB[jj][2]);
				Htotal[ii][jj] = kron_delta(RB[ii][1],RB[jj][1])*(E_JK + coeff*coupling);

		#Get eigenvalues and eigenvectors for the Hamiltonian
		EigVals,EigVecs = np.linalg.eigh(Htotal)
		idx = EigVals.argsort()
		eigvals = EigVals[idx]
		eigvecs = EigVecs[:,idx]
		energies.append(eigvals[0])
		gs = eigvecs[:,0]

		cossquared = np.dot(gs.T,np.dot(cossq,gs))
		alignment.append(cossquared)

	alignmentfxn = interpolate.interp1d(intensities,alignment)
	energyfxn    = interpolate.interp1d(intensities,energies)

	return alignmentfxn,energyfxn
