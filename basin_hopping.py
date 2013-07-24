#!/usr/bin/env python

#####################################################
# The first line must be that line.                 #
# It tells the computer this is a python script     #
#                                                   #
# The next thing we do is import all of the         #
# important python libraries that we will need      #
#####################################################



from commonfunctions import nearlySphericalAtom, makeBimetallic, distanceCenter, plusOrMinus
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

#some global variables for convenience:
listLength = 6
treeDepth = 2
FIREMoves = 3


######################################################
# Now we are going to define a few key functions     #
######################################################


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
  creationString2 = creationString + '.traj'

  baseAtom = nearlySphericalAtom(str(creationString),radius,elementN1+elementN2)
  baseAtom.set_pbc((1,1,1))
  baseAtom.set_cell((100,100,100))
  calc = QSC()
  baseAtom.set_calculator(calc)

  baseAtom = makeBimetallic(baseAtom,numberOfType1+numberOfType2,elementN1,elementN2,percentType1)
  baseAtom = preventExplosions(baseAtom)
  baseAtom = preventExplosions(baseAtom)

  for x in range(listLength):
    bigKickResults.append(baseAtom.copy())

  for x in range(listLength/2):
    bigKickResults[x] = shake(bigKickResults[x])

  for x in range(listLength/2, listLength):
    bigKickResults[x] = switchAtoms(bigKickResults[x])

  for x in range(len(bigKickResults)):
    bigKickResults[x] = optimizeMolecule(bigKickResults[x],FIREMoves,creationString2)
    round1PE.append(bigKickResults[x].get_potential_energy())

  minimumPE = min(round1PE)
  minimumIndex = round1PE.index(minimumPE)
  bestAtom = bigKickResults[minimumIndex]

  minimaList = PickleTrajectory(str(creationString2),mode='a')
  minimaList.write(bestAtom)
  minimaList.close()

  smallKicks(bigKickResults,0,creationString2)


def smallKicks(moleculeList, inTreeLevel, creationString):
  inTreeLevel += 1

  if inTreeLevel == treeDepth:
    print " #YOLODUBSTEP  #WUBWUBWUB "
    return None

  else:
    list1, list2, list3, list4, list5 = [],[],[],[],[]
    p1, p2, p3, p4, p5 = [],[],[],[],[]
    nAtoms = moleculeList[0].get_number_of_atoms() #they should all have the same number

    for x in range(len(moleculeList)):
      list1.append(moleculeList[x].copy())
      list2.append(moleculeList[x].copy())
      list3.append(moleculeList[x].copy())
      list4.append(moleculeList[x].copy())
      list5.append(moleculeList[x].copy())

    for x in range(len(list1)): #but they should all be the same length
      list1[x] = shell_move(list1[x],random.randint(0,nAtoms))
      list2[x] = HighEnergyMove(list2[x])
      list3[x] = ball_move(list3[x],random.randint(0,nAtoms))
      list4[x] = smallSwitchAtoms(list4[x])
      list5[x] = moveAtoms(2,list5[x])


    for x in range(len(list1)): #again, should all be the same length
      list1[x] = optimizeMolecule(list1[x],FIREMoves,creationString)
      list2[x] = optimizeMolecule(list2[x],FIREMoves,creationString)
      list3[x] = optimizeMolecule(list3[x],FIREMoves,creationString)
      list4[x] = optimizeMolecule(list4[x],FIREMoves,creationString)
      list5[x] = optimizeMolecule(list5[x],FIREMoves,creationString)

    del moleculeList[:]
    del moleculeList

    #This is also decidedly unpythonic
    for x in range(len(list1)):
      p1.append(list1[x].get_potential_energy())
      p2.append(list2[x].get_potential_energy())
      p3.append(list3[x].get_potential_energy())
      p4.append(list4[x].get_potential_energy())
      p5.append(list5[x].get_potential_energy())

    #p1min = (MIN)imum PE.  p1mi = (M)inimum (I)ndex.  p1b = (B)est.
    p1min, p2min, p3min, p4min, p5min = min(p1),min(p2),min(p3),min(p4),min(p5)
    p1mi, p2mi, p3mi, p4mi, p5mi = p1.index(p1min), p2.index(p2min), p3.index(p3min), p4.index(p4min), p5.index(p5min)
    p1b, p2b, p3b, p4b, p5b = p1[p1mi],p2[p2mi],p3[p3mi],p4[p4mi],p5[p5mi]

    minimaList = PickleTrajectory(str(creationString),mode='a')
    minimaList.write(p1b)
    minimaList.write(p2b)
    minimaList.write(p3b)
    minimaList.write(p4b)
    minimaList.write(p5b)
    minimaList.close()

    smallKicks(list1,inTreeLevel,creationString)
    smallKicks(list2,inTreeLevel,creationString)
    smallKicks(list3,inTreeLevel,creationString)
    smallKicks(list4,inTreeLevel,creationString)
    smallKicks(list5,inTreeLevel,creationString)

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
      return molecule
    else:
      chosenOne = numpy.where(posList==atom1)
      molecule = ball_move(molecule,chosenOne)
      return molecule

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


def optimizeMolecule(molecule,NMoves,creationString):

  bestEnergy = 0.
  totalMinimaFound = 0
  sinceLastFind = 0

  calc = QSC()
  molecule.set_calculator(calc)
  minimaList = PickleTrajectory(str(creationString),mode='a', molecule)

  for i in range(NMoves):
        molecule = newMove(molecule)
        molecule = preventExplosions(molecule)
        molecule.center()
        sinceLastFind += 1
        # do a last optimization of the structure
        dyn = FIRE(molecule)
        dyn.run(steps = 2000)
        newEnergy = molecule.get_potential_energy()

        if (newEnergy < bestEnergy):
                bestEnergy = newEnergy
                optimizedMolecule = molecule
                line = str(totalMinimaFound) + "  " + str(molecule.get_potential_energy()) + "  " + str(i) +"\n"
                print line
                f = open('EnergyList.txt','a')
                f.write(line)
                f.close()
                minimaList.write(optimizedMolecule)
                totalMinimaFound += 1
                sinceLastFind = 0
        elif (sinceLastFind < 200):
          break


  minimaList.close()

  minimaList = PickleTrajectory(creationString,mode='r')

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
