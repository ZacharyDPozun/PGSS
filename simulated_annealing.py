#!/usr/bin/env python

#####################################################
# The first line must be that line.                 # 
# It tells the computer this is a python script     #
#                                                   #
# The next thing we do is import all of the         #
# important python libraries that we will need      #
#####################################################
import ase
import tsase
import math
from qsc import QSC
import numpy
from numpy import random
from copy import deepcopy
from ase import units, io
from tsase import md 
from ase.optimize import FIRE
from ase.io.trajectory import PickleTrajectory
from commonfunctions import *
from InputVariables import *

######################################################
# Now we are going to define a few key functions     #
######################################################

#T = max temp, N = number of steps, xfactor1 = as a decimal; the fraction of the number of steps until the first point,
#yfactor 1 = as a decimal, xfactor2 = as a decimal;the fraction of the number of steps until the second point, yfactor2 = as a decimal
def generateTemperatures(T, N, xfactor1, yfactor1, xfactor2, yfactor2):
    temperature = [0]*N
    k = numpy.log(0.5/4.) / (N*(xfactor1 - xfactor2))
    for i in range(0, int(round(N*xfactor1))):
        temperature[i] = ((T*(1-yfactor1))/(-N*xfactor1))*i + T
    for i in range(int(round(N*xfactor1)),int(round(N*xfactor2))):
        temperature[i] = temperature[int(round(N*xfactor1)-1)] * numpy.exp(-k * (i - N * xfactor1))
    for i in range(int(round(N*xfactor2)),N):
        temperature[i] = (-(temperature[int(round(N*xfactor2) -1)]) / (N * (1 - xfactor2))) * (i - N * xfactor2) + temperature[int(round(N*xfactor2) -1)]
    return temperature

##
## This function will return the current temperature
## as a function of the total steps (NSteps), max
## temperature (TMax) and our current step (CurrentStep)
##

def temperature(CurrentStep,temperatures): 
	newTemp = temperatures[CurrentStep + 1]
	return newTemp

########################################################
# Now we can actually run our simulated annealing      #
# First we need to define our initial variables        #
########################################################

NSteps=35000
InitialTemp=5000
bestEnergy = 0.
totalMinimaFound = 0
atoms = makeBimetallic('POSCAR',NAtoms,Atom1,Atom2,CompAtom1)
filename = str(Atom1) + '_' + str(Atom2) + '_' + str(CompAtom1) + r'.traj' + str(NAtoms)

for i in range(NRuns): 
	xfactor1 = numpy.random.random()
	yfactor1 = numpy.random.random()
	xfactor2 = (numpy.random.random() * (1. - xfactor1)) + xfactor1
	yfactor2 = numpy.random.random() * yfactor1
	temperatures = generateTemperatures(InitialTemp,NSteps,xfactor1,yfactor1,xfactor2, yfactor2)
	# do our annealing according to the schedule set above
	atoms.center() # recenter the atoms every time, just in case
	for n in range(NSteps - 1):
		currentTemp = temperature(n, temperatures)
		dyn = tsase.md.nvtandersen(atoms, 5 * units.fs, units.kB * currentTemp)
		dyn.run(1)
	# do a last optimization of the structure
	dyn = FIRE(atoms)
	dyn.run()
	newEnergy = atoms.get_potential_energy()
	if (newEnergy < bestEnergy):
		bestEnergy = newEnergy
		line = str(totalMinimaFound) + "  " + str(atoms.get_potential_energy()) + "  " + str(i) +"\n"
		print line
		f = open('EnergyList.txt','a')
		f.write(line)
		f.close()
		minimaList = PickleTrajectory(filename,mode='a')
		minimaList.write(atoms)
		minimaList.close()
		totalMinimaFound += 1
                    
minimaList = PickleTrajectory(filename,mode='r')
atomslist = [atom for atom in minimaList]
ase.io.write('movie.xyz',atomslist,format='xyz') # write a movie file of our dynamics
