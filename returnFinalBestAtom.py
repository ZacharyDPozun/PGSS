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

def returnFinalBestAtom(trajectoryFile):
  candidates = PickleTrajectory(str(trajectoryFile),mode='r')
  molecules = [molecule for molecule in candidates]
  potentialEnergyList = [molecule.get_potential_energy() for molecule in molecules]
  minimumPotenialEnergy = min(potentialEnergyList)
  indexOfMinimumPotentialEnergy = potentialEnergyList.index(minimumPotenialEnergy)
  bestAtom = molecules[indexOfMinimumPotentialEnergy].copy()
  print "Yes We Can"
  return bestAtom
