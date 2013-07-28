#!/usr/bin/env python

import ase
from qsc import QSC
import numpy
from numpy import random
from ase import io, optimize, md, units, Atoms
from ase.optimize import FIRE
from ase.optimize.bfgslinesearch import BFGSLineSearch
from ase.io.trajectory import PickleTrajectory
#from inputvariables import *
from commonfunctions import *
import random
from copy import deepcopy

### Input Variables

generations = 1
NPsPerGen = 3
FractionCost = 0.5
f = open('test.txt','w')

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
	opt = BFGSLineSearch(core)
	opt.run(steps=7)
	opt =  FIRE(core)
	opt.run(steps=2000)
	return core


#### set it all up

ElementsDictionary = {1:'Rh',2:'Ir',3:'Ni',4:'Pd',5:'Pt',6:'Cu',7:'Ag',8:'Au'}
atoms= createparticle('Pt','Pt',3)
Ptcohesive=coreshellCohesive(atoms)
atoms =  createparticle('Cu','Cu',3)
Mincost=totalCost(atoms)

# set up generation 1
ParticleList = []
CostList = []
CohesiveList = []
f.write("Generation 0 " + "\n")
f.write('\n')
for j in range(NPsPerGen):
	type1 = numpy.random.randint(1,9) 
	type2 = numpy.random.randint(1,9)
	numbertype1 = numpy.random.randint(86)
	ParticleList.append(createparticle(ElementsDictionary[type1],ElementsDictionary[type2],numbertype1))
	ParticleList[j].gene = (ElementsDictionary[type1],ElementsDictionary[type2],numbertype1)
	CostList.append(totalCost(ParticleList[j]))
	CohesiveList.append(coreshellCohesive(ParticleList[j]))
	
for i in range(generations):
	ObjectFunction = []
	for j in range(NPsPerGen):
		costpart = (CostList[j] - Mincost) * 0.1
		cohesivepart = CohesiveList[j] - Ptcohesive
		ObjectFunction.append(FractionCost * costpart + (1 - FractionCost) * cohesivepart)
		f.write(ParticleList[j].gene[0]+ " " + ParticleList[j].gene[1] + " " + str(ParticleList[j].gene[2]) +  " " + str(CostList[j]) + " " + str(CohesiveList[j]) + " " + str(ObjectFunction[j]) + "\n")
	rankings = argsort(ObjectFunction)
	f.write("Average cost: " + str(sum(CostList) / NPsPerGen) + ", Average Core-Shell Binding: " + str(sum(CohesiveList) / NPsPerGen) + "\n")
	f.write('\n')
	




f.close()
