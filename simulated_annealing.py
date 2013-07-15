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
from qsc import QSC
import numpy
from numpy import random
from copy import deepcopy
from ase import units, io
from tsase import md 
from ase.optimize import FIRE
from ase.io.trajectory import PickleTrajectory
from commonfunctions import *

######################################################
# Now we are going to define a few key functions     #
######################################################

##
## This function will return the current temperature
## as a function of the total steps (NSteps), max
## temperature (TMax) and our current step (CurrentStep)
##
def temperature(NSteps, TMax, CurrentStep): 
	newTemp = TMax * (1 - (CurrentStep / NSteps))  # Linear Schedule
	return newTemp

	
########################################################
# Now we can actually run our simulated annealing      #
# First we need to define our initial variables        #
########################################################

NSteps = 30000 # Total steps PER anneal
InitialTemp = 1000 # Starting temperature
NRuns = 500 # Total anneals to do
bestEnergy = 0.
totalMinimaFound = 0


###########################################################
# Here is the main body of the code. We'll load our atoms #
# object and attach a calculator (QSC). Then we'll create #
# two nested loops to run NRuns anneals for NSteps each.  #
###########################################################
atoms = makeBimetallic('POSCAR',100,78,79,.50)
calc = QSC()
atoms.set_calculator(calc)
minimaList = PickleTrajectory('Pt75Au25.traj',mode='a')

for i in range(NRuns): 
	# do our annealing according to the schedule set above
	atoms.center() # recenter the atoms every time, just in case
	for n in range(NSteps):
		currentTemp = temperature(NSteps,InitialTemp, n)
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
		minimaList.write(atoms)
		totalMinimaFound += 1
		
minimaList.close()
minimaList = PickleTrajectory('Pt75Au25.traj',mode='r')

atomslist = [atom for atom in minimaList]
ase.io.write('movie.xyz',atomslist,format='xyz') # write a movie file of our dynamics

minimaList.close()

