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
from ase import io, optimize, md, units
from ase.optimize import FIRE
from ase.io.trajectory import PickleTrajectory
from commonfunctions import *

######################################################
# Now we are going to define a few key functions     #
######################################################

##
## This function makes our move. We give a number of atoms
## to move and it moves that many atoms each in a direction
## chosen from a Gaussian distribution in X Y and Z
##
def moveAtoms(numbertomove,atoms):
	totalAtoms = len(atoms)
	positions = atoms.get_positions()
	for i in range(numbertomove):
		atomtomove = numpy.random.randint(totalAtoms) # choose a random atom number
		displacementX = numpy.random.normal(scale = 1.0)
		displacementY = numpy.random.normal(scale = 1.0)
		displacementZ = numpy.random.normal(scale = 1.0)
		positions[atomtomove,0] += displacementX
		positions[atomtomove,1] += displacementY
		positions[atomtomove,2] += displacementZ
	atoms.set_positions(positions)
	print atoms.get_potential_energy()
	return atoms
	
def preventExplosions(atoms):
	positions = atoms.get_positions()
	for i in range(len(atoms)):
		for j in range(len(atoms)):
			distance = numpy.sqrt(numpy.dot((positions[j] - positions[i]),(positions[j] - positions[i])))
			if ((distance < 1.5) and (i != j)):
				print i, j, distance
				unitDirection = (positions[j] - positions[i]) / distance
				positions[j] += unitDirection
	atoms.set_positions(positions)
	return atoms
	
	
########################################################
# Now we can actually run our basin hopping            #
# First we need to define our initial variables        #
########################################################

# total moves made
NMoves = 10000
bestEnergy = 0.
totalMinimaFound = 0
sinceLastFind = 0

###########################################################
# Here is the main body of the code. We'll load our atoms #
# object and attach a calculator (QSC). 
###########################################################
atoms = makeBimetallic('POSCAR',100,78,79,0.5)
calc = QSC()
atoms.set_calculator(calc)
minimaList = PickleTrajectory('Pt75Au25.traj',mode='a')

for i in range(NMoves):
	numbertomove = numpy.random.randint(len(atoms))
	atoms = moveAtoms(numbertomove,atoms)
	atoms.center()
	atoms = preventExplosions(atoms)
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
		sinceLastFind = 0
	elif (sinceLastFind < 200): # if we haven't found a new minimum in 200 tries, start over
		atoms = ase.io.read('POSCAR')
		calc = QSC()
		atoms.set_calculator(calc)
		atoms = makeBimetallic(atoms,79,78,0.25)
		sinceLastFind = 0


minimaList.close()		
minimaList = PickleTrajectory('Pt75Au25.traj',mode='r')

atomslist = [atom for atom in minimaList]
ase.io.write('movie.xyz',atomslist,format='xyz') # write a movie file of our dynamics

minimaList.close()
