import ase
import tsase
from qsc import QSC
import numpy
import random
from numpy import random
from ase import io, optimize, md, units, Atoms
from ase.optimize import FIRE
from ase.io.trajectory import PickleTrajectory
from ase.md import VelocityVerlet

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

def create_sample_atom(inNumberOfAtoms,inAtomType):
  """Creates a mostly spherical atom with inNumberOfAtoms atoms
  and all of element inAtomType.  inAtomType MUST BE A STRING.
  For example, 'Pt79' for Platinum and 79 atoms."""
  positionList = []
  halfOfRadius = 3
  standardDeviation = 1.65
  for x in range(0,inNumberOfAtoms):
    atomX = random.gauss(halfOfRadius,standardDeviation)*random.choice([-1,1])
    atomY = random.gauss(halfOfRadius,standardDeviation)*random.choice([-1,1])
    atomZ = random.gauss(halfOfRadius,standardDeviation)*random.choice([-1,1])
    positionList.append((atomX,atomY,atomZ))
  sample_atom = Atoms(inAtomType,numpy.asarray(positionList))
  return sample_atom
