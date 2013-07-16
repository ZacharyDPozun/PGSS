#!/usr/bin/env python

#####################################################
# The first line must be that line.                 # 
# It tells the computer this is a python script     #
#                                                   #
# The next thing we do is import all of the         #
# important python libraries that we will need      #
#####################################################

from basin_hopping import *
from commonfunctions import *
import ase
import tsase
from qsc import QSC
import numpy
from numpy import random
from copy import deepcopy
from ase import units, io, Atoms
from tsase import md 
from ase.optimize import FIRE
from ase.io.trajectory import PickleTrajectory
import math


######################################################
# Now we are going to define a few key functions     #
######################################################

#T = max temp, N = number of steps, xfactor1 = as a decimal; the fraction of the number of steps until the first point,
#yfactor 1 = as a decimal, xfactor2 = as a decimal;the fraction of the number of steps until the second point, yfactor2 = as a decimal
#the "factor" arguments represent ratios, so they should be less than 1. The factor arguments define the endpoints of the "pieces" of the piecewise function
def generateTemperatures(T, N, xfactor1, yfactor1, xfactor2, yfactor2):
    temperature = [0]*N
    k = (math.log(yfactor1*T) - math.log(yfactor2*T))/(N*(xfactor1 - xfactor2))
    A = (yfactor1*T)/(math.e**(k*(xfactor1*N)))
    for i in range(0, int(round(N*xfactor1))+1):
        temperature[i] = ((T*(1-yfactor1))/(-N*xfactor1))*i + T
    for i in range(int(round(N*xfactor1)+1),int(round(N*xfactor2))+1):
        temperature[i] = A * (math.e**(k*i))
    for i in range(int(round(N*xfactor2)+1),N):
        temperature[i] = ((T*yfactor2)/(N*(xfactor2-1)))*(i-N)
    return temperature

##
## This function will return the current temperature
## as a function of the total steps (NSteps), max
## temperature (TMax) and our current step (CurrentStep)
##

def temperature(CurrentStep,temperatures): 
	newTemp = temperatures[CurrentStep + 1]
	return newTemp

	
##
## This function takes our input atoms object and turns  
## it into a bimetallic alloy based on an input fraction 
## 
def makeBimetallic(atoms,element1,element2,fractionElement1): 
	numberAtoms = len(atoms)
	numberElementOne = numpy.int(len(atoms) * fractionElement1) # use integers not floats! 
	atomicNumberArray = numpy.ones(numberAtoms) * element2
	for i in range(numberElementOne):
		atomicNumberArray[i] = element1
	atoms.set_atomic_numbers(atomicNumberArray)
	return atoms

	
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

#atoms = makeBimetallic('POSCAR',100,78,79,.50)
#calc = QSC()
#atoms.set_calculator(calc)
#minimaList = PickleTrajectory('Pt75Au25.traj',mode='a')

#change the arguments in the line below to vary the temperature schedule
#temperatures = generateTemperatures(InitialTemp,NSteps,0.0,1.0,0.85, 0.30)

##for i in range(NRuns):
##	# do our annealing according to the schedule set above
##	atoms.center() # recenter the atoms every time, just in case
##	for n in range(NSteps - 1):
##		currentTemp = temperature(NSteps,InitialTemp, n)
##		dyn = tsase.md.nvtandersen(atoms, 5 * units.fs, units.kB * currentTemp)
##		dyn.run(1)
##	# do a last optimization of the structure
##	dyn = FIRE(atoms)
##	dyn.run()
##	newEnergy = atoms.get_potential_energy()
##	if (newEnergy < bestEnergy):
##		bestEnergy = newEnergy
##		line = str(totalMinimaFound) + "  " + str(atoms.get_potential_energy()) + "  " + str(i) +"\n"
##		print line
##		f = open('EnergyList.txt','a')
##		f.write(line)
##		f.close()
##		minimaList.write(atoms)
##		totalMinimaFound += 1
##
##minimaList.close()
##minimaList = PickleTrajectory('Pt75Au25.traj',mode='r')
##
##atomslist = [atom for atom in minimaList]
##ase.io.write('movie.xyz',atomslist,format='xyz') # write a movie file of our dynamics
##
##minimaList.close()

