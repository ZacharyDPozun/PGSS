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

def returnFinalBestAtom(trajectoryFile,imageFile):
  candidates = PickleTrajectory(str(trajectoryFile),mode='r')
  molecules = [molecule for molecule in candidates]
  potentialEnergyList = [molecule.get_potential_energy() for molecule in molecules]
  minimumPotentialEnergy = min(potentialEnergyList)
  indexOfMinimumPotentialEnergy = potentialEnergyList.index(minimumPotentialEnergy)
  bestAtom = molecules[indexOfMinimumPotentialEnergy].copy()
  
  bestAtom.set_calculator(QSC())
  
  bookFile = open('BestEnergies.txt', mode = 'a')
  line = str(trajectoryFile) + " " + str(bestAtom.get_potential_energy()) + " " + str(indexOfMinimumPotentialEnergy) + "\n"
  bookFile.write(line)
  bookFile.close()
  
  imageFileName = str(imageFile) + '.png'
  
  ase.io.write(str(imageFileName),bestAtom,format = 'png')

  return bestAtom

def writeBestEnergyIntoFile(trajectoryFile,molecule):
  bookkeepingFile = open('BestEnergies.txt', mode = 'a')
  candidates = PickleTrajectory(str(trajectoryFile),mode='r')
  molecules = [molecule for molecule in candidates]
  potentialEnergyList = [molecule.get_potential_energy() for molecule in molecules]
  minimumPotenialEnergy = min(potentialEnergyList)
  indexOfMinimumPotentialEnergy = potentialEnergyList.index(minimumPotenialEnergy)
  line = str(trajectoryFile) + str(molecule.get_potential_energy()) + str(indexOfMinimumPotentialEnergy)
  bookkeepingFile.write(line)
  bookkeepingFile.close()
