import ase
import tsase
from qsc import QSC
import numpy
from ase import io

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

def nearestNeighbor(molecule):
    MinDistanceToX = []
    positionList = molecule.get_positions()
    for i in range(0,len(positionList)):
        distanceList = []
        for j in range (0,len(positionList)):
            distance = numpy.linalg.norm(positionList[i]-positionList[j])
            pair = (distance,i,j)
            distanceList.append(pair)
        distanceList.sort()
        MinDistanceToX.append(distanceList[1])

    return MinDistanceToX
