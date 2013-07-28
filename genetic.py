#!/usr/bin/env python

import ase
from qsc import QSC
import numpy
from ase import io, optimize, md, units, Atoms
from ase.optimize import FIRE
from ase.io.trajectory import PickleTrajectory
#from inputvariables import *
from commonfunctions import *
import random
from copy import deepcopy

### Input Variables



### Set up functions

def createparticle(type1, type2, numbertype1):
	calc = QSC()
	atoms = ase.io.read('POSCAR_genetic',format='vasp')	
	shell = deepcopy(atoms)
	core = deepcopy(atoms)
	for i in range(85,225):
		core.pop(85)
	for i in range(85):
		shell.pop(0)
	dummyarray = core.get_chemical_symbols()
	for i in range(numbertype1):
		dummyarray[i] = type1
	for i in range(numbertype1,85):
		dummyarray[i] = type2
	random.shuffle(dummyarray)
	core.set_chemical_symbols(dummyarray)
	for i in range(len(shell)):
		core.append(shell.pop(0))
	core.set_calculator(calc)
	opt = FIRE(core)
	opt.run(steps=10000)
	return core

atoms= createparticle('Pt','Pt',3)

Ptcohesive=coreshellCohesive(atoms)
Ptcost=totalCost(atoms)

print Ptcohesive,Ptcost