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

def rattleAtoms(atoms):
        seed = numpy.random.randint(100)
        atoms.rattle(stdev=0.5,seed=seed)
        return atoms

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
		
def shell_move(inAtom,atomIndex):  
  #  we're going to be changing the position of atomIndex inside inAtom
  #  make sure that you remove any crazy outliers before you do this
  #  or else it'll just make a bunch more outliers, which is a poor idea
  #  make sure atomIndex comes from range(len(inAtom.get_positions())) so we don't get out of bounds
  try:
    inCOM = inAtom.get_center_of_mass()
    inDistances = distanceCenter(inAtom)
    
    ninetyNinthRadius = stats.scoreatpercentile(inDistances,99)
    ninetyFifthRadius = stats.scoreatpercentile(inDistances,95) 
    outerFourRadius = ninetyNinthRadius - ninetyFifthRadius
    
    randomNewRadius = random.gauss( (ninetyNinthRadius+ninetyFifthRadius)/2 , (ninetyNinthRadius - ninetyFifthRadius)/2 )
    xFromCenter = random.uniform(0,randomNewRadius)
    randomNewRadius = ((randomNewRadius**2) - (xFromCenter**2))**0.5
    yFromCenter = random.uniform(0,randomNewRadius)
    zFromCenter = ((randomNewRadius**2) - (yFromCenter**2))**0.5
    
    newXPosition = inCOM[0] + random.choice([-1,1])*xFromCenter
    newYPosition = inCOM[1] + random.choice([-1,1])*yFromCenter
    newZPosition = inCOM[2] + random.choice([-1,1])*zFromCenter
    
    positionArray = inAtom.get_positions()
    positionArray[atomIndex] = (newXPosition,newYPosition,newZPosition)
    inAtom.set_positions(positionArray)
    
    return inAtom
    
  except IndexError:
    print "The index of the atom you wanted to move is too high or too low."
    print "Please check your function call of shell_move(a,b)"
    print "-Jeff"


minimaList.close()		
minimaList = PickleTrajectory('Pt75Au25.traj',mode='r')

atomslist = [atom for atom in minimaList]
ase.io.write('movie.xyz',atomslist,format='xyz') # write a movie file of our dynamics

minimaList.close()
