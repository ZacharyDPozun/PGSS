import ase
import tsase
from qsc import QSC
from numpy import *
import numpy
from ase import io, optimize, md, units, Atoms
from ase.optimize import FIRE
from ase.io.trajectory import PickleTrajectory
from ase.md import VelocityVerlet


def create_sample_atom(inNumberOfAtoms,inAtomType):
  """Creates a mostly spherical atom with inNumberOfAtoms atoms
  and all of element inAtomType.  inAtomType MUST BE A STRING.
  For example, 'Pt79' for Platinum and 79 atoms."""
  positionList = []
  halfOfRadius = 3
  standardDeviation = 1.65
  for x in range(0,inNumberOfAtoms):
    atomX = random.normal(halfOfRadius,standardDeviation)*plusOrMinus()
    atomY = random.normal(halfOfRadius,standardDeviation)*plusOrMinus()
    atomZ = random.normal(halfOfRadius,standardDeviation)*plusOrMinus()
    positionList.append((atomX,atomY,atomZ))
  sample_atom = Atoms(inAtomType,numpy.asarray(positionList))
  return sample_atom
  
def nearlySphericalAtom(definingString,inRadius,number):
  """definingString should be something like " 'Pt80' ", for example.
  inRadius should be the radius that you wish to have for the atom.
  number should be the number of atoms that will be in the molecule."""
  positionList = []
  for x in range(number):
    xDistance = random.uniform(0,inRadius)*plusOrMinus()
    remainingX = ((inRadius**2) - (xDistance**2))**0.5
    yDistance = random.uniform(0,remainingX)*plusOrMinus()
    zDistance = (((remainingX**2) - (yDistance**2))**0.5)*plusOrMinus() + random.normal(0,0.1)
    coordinates = (xDistance,yDistance,zDistance)
    positionList.append(coordinates)
  newAtom = Atoms(definingString,positionList)
  return newAtom

def makeBimetallic(filename,numberAtoms,element1,element2,fractionElement1):
        atoms = ase.io.read(filename,format='vasp')
	for i in range(numberAtoms, len(atoms)):
		atoms.pop()
        numberElementOne = numpy.int(len(atoms) * fractionElement1) # use integers not floats!
        atomicNumberArray = numpy.ones(numberAtoms) * element2
        for i in range(numberElementOne):
                atomicNumberArray[i] = element1
        atoms.set_atomic_numbers(atomicNumberArray)
	atoms.center()
        return atoms
        
def distanceCenter(atoms):
	distanceArray = numpy.zeros(len(atoms))
	cx = atoms.get_center_of_mass()[0]
	cy = atoms.get_center_of_mass()[0]
	cz = atoms.get_center_of_mass()[0]
	for i in range(len(atoms)):
                px = atoms.get_positions()[i,0]
		py = atoms.get_positions()[i,1]
		pz = atoms.get_positions()[i,2]
		distance = ((px-cx)**2+(py-cy)**2+(pz-cz)**2)**0.5
		distanceArray[i] = distance
	return distanceArray

def visualize_atom(inAtom):
  positionList = inAtom.get_positions()
  xList, yList, zList = [],[],[]
  for x in range(len(positionList)):
    xList.append(positionList[x][0])
    yList.append(positionList[x][1])
    zList.append(positionList[x][2])
  #
  #
  fig = pylab.figure()
  ax = Axes3D(fig)
  ax.scatter(xList,yList,zList)
  pyplot.show()

def plusOrMinus():
    a = random.uniform(0,1)
    if (a>0.5):
        return 1
    else:
        return -1
        
def totalCost(molecule):
  """Pass in the molecule and we return the price per gram"""
  totalAtoms = molecule.get_number_of_atoms()
  #price per gram of common metals
  Al = 0.00220462
  Ag = 0.62308
  Au = 41.25
  Cu = 0.00661
  Ni = 0.01543
  Pd = 24.25085
  Pt = 45.3008
  #counters for the common metals
  alCounter = 0
  agCounter = 0
  auCounter = 0
  cuCounter = 0
  niCounter = 0
  pdCounter = 0
  ptCounter = 0
  #get the chemical symbols of our molecule and count them up
  symbolsList = molecule.get_chemical_symbols()
  alCounter = symbolsList.count("Al")
  agCounter = symbolsList.count("Ag")
  auCounter = symbolsList.count("Au")
  cuCounter = symbolsList.count("Cu")
  niCounter = symbolsList.count("Ni")
  pdCounter = symbolsList.count("Pd")
  ptCounter = symbolsList.count("Pt")
  #find the percentage for each metal
  x = alCounter/totalAtoms
  y = agCounter/totalAtoms
  z = auCounter/totalAtoms
  w = cuCounter/totalAtoms
  r = niCounter/totalAtoms
  n = pdCounter/totalAtoms
  m = ptCounter/totalAtoms
  #finally find the cost per gram for the entire molecule
  cost = (Al*x + Ag*y + Au*z + Cu*w + Ni*r + Pd*n + Pt*m)
  return cost
