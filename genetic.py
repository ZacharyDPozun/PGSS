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
import os

### Input Variables

generations = 5
NPsPerGen = 35
FractionCost = 0.5
f = open('Output.txt','w')

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
	opt.run(steps=8)
	opt =  FIRE(core)
	opt.run(steps=2000)
	return core

def breed(particles,objects,gen):
	sortedIndices = numpy.argsort(numpy.asarray(objects))
	particlestobreed = []
	explambda = numpy.log(10) / (NPsPerGen - 1)
	for i in range(len(particles)):
		trialvar = numpy.random.random()
		if (trialvar <= numpy.exp(-1 * explambda * i)):
			particlestobreed.append(particles[sortedIndices[i]])
	NewGen = []
	for i in range(NPsPerGen):
		first = numpy.random.randint(len(particlestobreed))
		second = numpy.random.randint(len(particlestobreed))
		gene1 = numpy.random.randint(2)
		gene2 = numpy.random.randint(2)
		if (gene1 == 0):
			type1 = particlestobreed[first].gene[0]
		else:	
			type1 = particlestobreed[second].gene[0]
		if (gene2 == 0):
			type2 = particlestobreed[first].gene[1]
		else:
			type2 = particlestobreed[second].gene[1]
		if (particlestobreed[first].gene[2] < particlestobreed[second].gene[2]):
			minval = particlestobreed[first].gene[2]
			maxval = particlestobreed[second].gene[2]
		else:
			minval = particlestobreed[second].gene[2]
			maxval = particlestobreed[first].gene[2]
		numbertype1 = numpy.random.randint(minval, maxval + 1)
		NewGen.append(createparticle(type1,type2,numbertype1))
		traj = PickleTrajectory('generation'+str(gen)+ '_' + str(i) + '.traj','w')
		NewGen[i].gene = (type1,type2,numbertype1)
        	traj.write(ParticleList[i])
		traj.close()
	return NewGen

def mutate(atoms):
	type = numpy.random.randint(4)
	if (type == 0):
		newtype = numpy.random.randint(8) + 1
		type2 = atoms.gene[1]
		numbertype1 = atoms.gene[2]
		atoms = createparticle(ElementsDictionary[newtype], type2,numbertype1)
		atoms.gene = (ElementsDictionary[newtype],type2,numbertype1)
	elif (type == 1):
                newtype	= numpy.random.randint(8) + 1
                type1 =	atoms.gene[0]
                numbertype1 = atoms.gene[2]
                atoms =	createparticle(type1,ElementsDictionary[newtype],numbertype1)                 
                atoms.gene = (type1,ElementsDictionary[newtype],numbertype1)
	elif (type == 2):
		type1 = atoms.gene[0]
		type2 = atoms.gene[1]
		numbertype1 = atoms.gene[2]
		atoms = createparticle(type2,type1,numbertype1) 
                atoms.gene = (type2,type1,numbertype1)
	elif (type == 3):
		newtype1 = numpy.random.randint(8) + 1
		newtype2 = numpy.random.randint(8) + 1
		numbertype1 = numpy.random.randint(86)
		atoms = createparticle(ElementsDictionary[newtype1],ElementsDictionary[newtype2],numbertype1)
		atoms.gene = (ElementsDictionary[newtype1],ElementsDictionary[newtype2],numbertype1)
	return atoms
		
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
for j in range(NPsPerGen):
	type1 = numpy.random.randint(1,9) 
	type2 = numpy.random.randint(1,9)
	numbertype1 = numpy.random.randint(86)
	ParticleList.append(createparticle(ElementsDictionary[type1],ElementsDictionary[type2],numbertype1))
	ParticleList[j].gene = (ElementsDictionary[type1],ElementsDictionary[type2],numbertype1)
	filename = 'generation0_' + str(j) + ".traj"
	traj = PickleTrajectory(filename,'w')
	traj.write(ParticleList[j])
	traj.close()
	CostList.append(totalCost(ParticleList[j]))
	CohesiveList.append(coreshellCohesive(ParticleList[j]))
traj.close()	

for i in range(generations):
	ObjectFunction = []
	CostList = []
	CohesiveList = []
	for j in range(NPsPerGen):
		CostList.append(totalCost(ParticleList[j]))
        	CohesiveList.append(coreshellCohesive(ParticleList[j]))
		costpart = (CostList[j] - Mincost) * 0.1
		cohesivepart = numpy.abs(CohesiveList[j] - Ptcohesive)
		ObjectFunction.append(FractionCost * costpart + (1 - FractionCost) * cohesivepart)
		f.write(ParticleList[j].gene[0]+ " " + ParticleList[j].gene[1] + " " + str(ParticleList[j].gene[2]) +  " " + str(CostList[j]) + " " + str(CohesiveList[j]) + " " + str(ObjectFunction[j]) + "\n")
	f.write("Average cost: " + str(sum(CostList) / NPsPerGen) + ", Average Core-Shell Binding: " + str(sum(CohesiveList) / NPsPerGen) + "\n")
	f.write('\n')
	## New Generation starts here
	NewGen = breed(ParticleList,ObjectFunction, (i+1))
	f.write("Generation " + str(i + 1) + '\n')
	for j in range(NPsPerGen):
		tryval = numpy.random.random()
		if (tryval <= 0.1):
			NewGen[j] = mutate(NewGen[j])
	ParticleList = NewGen

# analyze the final generation
ObjectFunction = []
CostList = []
CohesiveList = []
for j in range(NPsPerGen):
	CostList.append(totalCost(ParticleList[j]))
	CohesiveList.append(coreshellCohesive(ParticleList[j]))
	costpart = (CostList[j] - Mincost) * 0.1
	cohesivepart = CohesiveList[j] - Ptcohesive
	ObjectFunction.append(FractionCost * costpart + (1 - FractionCost) * cohesivepart)
	f.write(ParticleList[j].gene[0]+ " " + ParticleList[j].gene[1] + " " + str(ParticleList[j].gene[2]) +  " " + str(CostList[j]) + " " + str(CohesiveList[j]) + " " + str(ObjectFunction[j]) + "\n")
f.write("Average cost: " + str(sum(CostList) / NPsPerGen) + ", Average Core-Shell Binding: " + str(sum(CohesiveList) / NPsPerGen) + "\n")
f.write('\n')
f.close()
