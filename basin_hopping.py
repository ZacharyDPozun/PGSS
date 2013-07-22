#!/usr/bin/env python

#####################################################
# The first line must be that line.                 #
# It tells the computer this is a python script     #
#                                                   #
# The next thing we do is import all of the         #
# important python libraries that we will need      #
#####################################################

from commonfunctions import nearlySphericalAtom, makeBimetallic, distanceCenter
import ase
import tsase
from qsc import QSC
import numpy
from numpy import *
from ase import io, optimize, md, units, Atoms
from ase.optimize import FIRE
from ase.io.trajectory import PickleTrajectory
from ase.md import VelocityVerlet
from scipy import stats
import operator
from matplotlib import pyplot
import pylab
from mpl_toolkits.mplot3d import Axes3D
from InputVariables import *


######################################################
# Now we are going to define a few key functions     #
######################################################

##
## This function makes our move. We give a number of atoms
## to move and it moves that many atoms each in a direction
## chosen from a Gaussian distribution in X Y and Z
##

def mainBasinLoop(symbol1, symbol2, elementN1, elementN2, numberOfType1, numberOfType2, radius, percentType1):
  """symbol1 and symbol2 should be STRINGS.
  elementN1 and elementN2 should be ATOMIC NUMBERS of those elements.
  numberOfType1 and numberOfType2 are the actual # of those atoms in the molecule
  radius can be a decimal
  percent should be expressed as a number between 0 and 1.

  After that, it will create five new sets of 400 molecules,
  each set being the result of that first set having one of our small perturbations
  performed on it, which are also then minimized, perturbed, etc., until 50,000
  total original orientations have been minimized."""

  bigKickResults, round1PE = [], []

  creationString = symbol1 + str(elementN1) + symbol2 + str(elementN2)
  print creationString

  baseAtom = nearlySphericalAtom(str(creationString),radius,elementN1+elementN2)
  baseAtom.set_pbc((1,1,1))
  baseAtom.set_cell((100,100,100))

  baseAtom = makeBimetallic(baseAtom,numberOfType1+numberOfType2,elementN1,elementN2,percentType1)
  baseAtom = preventExplosions(baseAtom)
  baseAtom = preventExplosions(baseAtom)

  for x in range(200):
    bigKickResults.append(shake(baseAtom))
    bigKickResults.append(switchAtoms(baseAtom))

  for x in range(len(bigKickResults)):
    print bigKickResults[x]
    print x
    bigKickResults[x] = optimizeMolecule(bigKickResults[x],40)
    round1PE.append(bigKickResults[x].get_potential_energy())

  minimumPE = min(round1PE)
  minimumIndex = round1PE.index(minimumPE)
  bestAtom = bickKickResults[minimumIndex]

  creationString2 = creationString + '.traj'
  print creationString2
  minimaList = PickleTrajectory(str(creationString2),atoms=bestAtom,mode = 'a')
  minimaList.close()

  #smallKicks(bigKickResults,0)


def smallKicks(moleculeList, inTreeLevel):
  inTreeLevel += 1

  if inTreeLevel == 4:
    return None

  else:
    list1, list2, list3, list4, list5 = [],[],[],[],[]
    p1, p2, p3, p4, p5 = [],[],[],[],[]
    nAtoms = moleculeList[0].get_number_of_atoms()
    for molecule in moleculeList:

      #This is not pythonic
      list1.append(shell_move(molecule,random.randint(0,nAtoms)))
      list2.append(HighEnergyMove(molecule))
      list3.append(ball_move(molecule,random.randint(0,nAtoms)))
      list4.append(smallSwitchAtoms(molecule))
      list5.append(moveAtoms(2,molecule))

    for x in range(len(moleculeList)):
      pass
      #RUN THE OPTIMIZER THINGY AGAIN
      #list1[x] = resultOfOptimization(list1[x])
      #             .
      #             .
      #             .
      #listn[x] = resultOfOptimization(listn[x])

    #ZOMG SAVE YOUR MEMORYYYY
    del moleculeList

    #This is also decidedly unpythonic
    for x in range(len(list1)):
      p1.append(list1[x].get_potential_energy())
      p2.append(list2[x].get_potential_energy())
      p3.append(list3[x].get_potential_energy())
      p4.append(list4[x].get_potential_energy())
      p5.append(list5[x].get_potential_energy())

    #Pickle this part and make it not suck instead of doing this
    print min(p1),min(p2),min(p3),min(p4),min(p5)

    smallKicks(list1,inTreeLevel)
    smallKicks(list2,inTreeLevel)
    smallKicks(list3,inTreeLevel)
    smallKicks(list4,inTreeLevel)
    smallKicks(list5,inTreeLevel)

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

    randomNewRadius = random.normal( (ninetyNinthRadius+ninetyFifthRadius)/2 , (ninetyNinthRadius - ninetyFifthRadius)/2 )
    xFromCenter = random.uniform(0,randomNewRadius)
    randomNewRadius = ((randomNewRadius**2) - (xFromCenter**2))**0.5
    yFromCenter = random.uniform(0,randomNewRadius)
    zFromCenter = ((randomNewRadius**2) - (yFromCenter**2))**0.5

    newXPosition = inCOM[0] + plusOrMinus()*xFromCenter
    newYPosition = inCOM[1] + plusOrMinus()*yFromCenter
    newZPosition = inCOM[2] + plusOrMinus()*zFromCenter

    positionArray = inAtom.get_positions()
    positionArray[atomIndex] = (newXPosition,newYPosition,newZPosition)
    inAtom.set_positions(positionArray)

    return inAtom

  except IndexError:
    print "The index of the atom you wanted to move is too high or too low."
    print "Please check your function call of shell_move(a,b)"
    print "-Jeff"

def jolt(inAtom,magnitude):
  """Moves all the atoms in inAtom by a Gaussian distribution
  with mean of 0 and StdDev of magnitude.
  A small magnitude is a small kick,
  a larger magnitude is a big kick."""
  posList = inAtom.get_positions()
  for x in range(len(posList)):
    posList[x][0] = posList[x][0] + random.normal(0,magnitude)
    posList[x][1] = posList[x][1] + random.normal(0,magnitude)
    posList[x][2] = posList[x][2] + random.normal(0,magnitude)
  inAtom.set_positions(posList)
  return inAtom

def HighEnergyMove2(atoms):
    for cntr in range (len(atoms)):
        if (atoms[cntr].getNeighboringAtoms(.5) > 3):
            decideSingleMove()
        elif (atoms[cntr].getNeighboringAtoms(.5) <= 3):
            shake(atoms)
            

def ball_move(inAtom,atomIndex):
  """takes an atom defined by atomIndex inside of inAtom
  and moves it somewhere within the core of the atom randomly.
  Atoms will almost always end up inside the sphere which
  contains 85% of the atoms, centered at the center of mass."""
  #  we're going to be changing the position of atomIndex inside inAtom
  #  we'll take atom of index atomIndex and throw it somewhere inside the core
  #  make sure that you remove any crazy outliers before you do this
  #  or else it'll just make a bunch more outliers, which is a poor idea
  try:
    #get all the distances from the center of mass
    inCOM = inAtom.get_center_of_mass()
    inDistances = distanceCenter(inAtom)
    #figure out the distance from the core to the 85th percentile
    #we'll consider "the core" to be the sphere which contains 85% of the atoms
    eightyFifthRadius = stats.scoreatpercentile(inDistances,85)
    #pick a new distance from center somewhere inside that 85th percentile limit
    randomNewRadius = random.normal(eightyFifthRadius/2, eightyFifthRadius/3 )
    xFromCenter = random.uniform(0,randomNewRadius)
    randomNewRadius = ((randomNewRadius**2) - (xFromCenter**2))**0.5
    yFromCenter = random.uniform(0,randomNewRadius)
    zFromCenter = ((randomNewRadius**2) - (yFromCenter**2))**0.5
    newXPosition = inCOM[0] + plusOrMinus()*xFromCenter
    newYPosition = inCOM[1] + plusOrMinus()*yFromCenter
    newZPosition = inCOM[2] + plusOrMinus()*zFromCenter
    positionArray = inAtom.get_positions()
    positionArray[atomIndex] = (newXPosition,newYPosition,newZPosition)
    inAtom.set_positions(positionArray)
    return inAtom
  except IndexError:
    print "The index of the atom you wanted to move is too high or too low."
    print "Please check your function call of ball_move(a,b)"
    print "-Jeff"

def shake(atoms):
        seed = random.randint(100)
        atoms.rattle(stdev=0.5,seed=seed)
        return atoms

def switchAtoms(atoms):
        numbers = atoms.get_atomic_numbers()
        random.shuffle(numbers)
        atoms.set_atomic_numbers(numbers)
        return atoms

def smallSwitchAtoms(atoms):
  numbers = atoms.get_atomic_numbers()
  firstAtom = random.randint(0,len(numbers))
  secondAtom = random.randint(0,len(numbers))
  while numbers[firstAtom] == numbers[secondAtom]:
    firstAtom = random.randint(0,len(numbers))
    secondAtom = random.randint(0,len(numbers))
  numbers[firstAtom],numbers[secondAtom] = numbers[secondAtom],numbers[firstAtom]
  atoms.set_atomic_numbers(numbers)
  return atoms

def moveAtoms(numbertomove,atoms):
        totalAtoms = len(atoms)
        positions = atoms.get_positions()
        for i in range(numbertomove):
                atomtomove = random.randint(totalAtoms) # choose a random atom number
                displacementX = random.normal(scale = 1.0)
                displacementY = random.normal(scale = 1.0)
                displacementZ = random.normal(scale = 1.0)
                positions[atomtomove,0] += displacementX
                positions[atomtomove,1] += displacementY
                positions[atomtomove,2] += displacementZ
        atoms.set_positions(positions)
        print atoms.get_potential_energy()
        return atoms
        
def HighEnergyMove(molecule):
  posList = molecule.get_positions()
  for atom1 in posList:
    lessThanThreeA = 0
    for atom2 in posList:
      distanceBetween = numpy.linalg.norm(atom1-atom2)
      if distanceBetween < 3.0:
        lessThanThreeA += 1
    if lessThanThreeA < 3:
      pass
    else:
      chosenOne = posList.index(atom1)
      ball_move(molecule,chosenOne)

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

#bestEnergy = 0.
#totalMinimaFound = 0
#sinceLastFind = 0

###########################################################
# Here is the main body of the code. We'll load our atoms #
# object and attach a calculator (QSC).
###########################################################

##atoms = makeBimetallic('POSCAR',NAtoms,Atom1,Atom2,CompAtom1)
##calc = QSC()
##atoms.set_calculator(calc)
##minimaList = PickleTrajectory('Pt75Au25.traj',mode='a')

##for i in range(NMoves):
##      numbertomove = random.randint(len(atoms))
##      atoms = moveAtoms(numbertomove,atoms)
##      atoms.center()
##
##      atoms = preventExplosions(atoms)
##      # do a last optimization of the structure
##      dyn = FIRE(atoms)
##      dyn.run()
##      newEnergy = atoms.get_potential_energy()
##      if (newEnergy < bestEnergy):
##              bestEnergy = newEnergy
##              line = str(totalMinimaFound) + "  " + str(atoms.get_potential_energy()) + "  " + str(i) +"\n"
##              print line
##              f = open('EnergyList.txt','a')
##              f.write(line)
##              f.close()
##              minimaList.write(atoms)
##              totalMinimaFound += 1
##              sinceLastFind = 0
##      elif (sinceLastFind < 200): # if we haven't found a new minimum in 200 tries, start over
##              atoms = ase.io.read('POSCAR')
##              calc = QSC()
##              atoms.set_calculator(calc)
##              atoms = makeBimetallic(atoms,79,78,0.25)
##              sinceLastFind = 0




##minimaList.close()
##minimaList = PickleTrajectory('Pt75Au25.traj',mode='r')

##atomslist = [atom for atom in minimaList]
##ase.io.write('movie.xyz',atomslist,format='xyz') # write a movie file of our dynamics

##minimaList.close()


def optimizeMolecule(molecule,NMoves):

  bestEnergy = 0.
  totalMinimaFound = 0
  sinceLastFind = 0

  calc = QSC()
  molecule.set_calculator(calc)
  minimaList = PickleTrajectory('Au75Pt25.traj',mode='a')

  for i in range(NMoves):
        molecule = newMove(molecule)
        molecule = preventExplosions(molecule)
        molecule.center()
        sinceLastFind += 1
        # do a last optimization of the structure
        dyn = FIRE(molecule)
        dyn.run()
        newEnergy = molecule.get_potential_energy()

        if (newEnergy < bestEnergy):
                bestEnergy = newEnergy
                optimizedMolecule = molecule
                line = str(totalMinimaFound) + "  " + str(molecule.get_potential_energy()) + "  " + str(i) +"\n"
                print line
                f = open('EnergyList.txt','a')
                f.write(line)
                f.close()
		print str(optimizedMolecule)
		print str(molecule)
                #minimaList.write(molecule)
                totalMinimaFound += 1
                sinceLastFind = 0
        elif (sinceLastFind < 200):
          pass

  minimaList.close()
  minimaList = PickleTrajectory('Pt75Au25.traj',mode='r')

  atomslist = [atom for atom in minimaList]
  ase.io.write('movie.xyz',atomslist,format='xyz') # write a movie file of our dynamics

  minimaList.close()
  
  return optimizedMolecule

def newMove(molecule):
  decision = random.randint(1,6)
  nAtoms = molecule.get_number_of_atoms()
  if decision == 1:
    molecule = shell_move(molecule,random.randint(0,nAtoms))
  elif decision == 2:
    molecule = HighEnergyMove(molecule)
  elif decision == 3:
    molecule = ball_move(molecule,random.randint(0,nAtoms))
  elif decision == 4:
    molecule = smallSwitchAtoms(molecule)
  else:
    molecule = moveAtoms(2,molecule)

  return molecule

