#!/usr/bin/env python

import os

dirs = [i for i in range(10,201,2)]

for dir in dirs:
	os.mkdir(str(dir))
	os.chdir(str(dir)) 
	f = os.popen('cp ~/PGSS/POSCAR ~/PGSS/simulated_annealing.py ~/PGSS/commonfunctions.py .')
	f.close()
	f = open('InputVariables.py','w')
	f.write('#!/usr/bin/env python\n')
	f.write(' \n')
	f.write('NAtoms = ' + str(dir) + '\n')
	f.write('Atom1=78\n')
	f.write('Atom2=79\n')
	f.write('CompAtom1=0.5\n')
	f.write('NRuns = 200\n')
	f.write(' \n')
	f.close()
	f = open('anneal.sub','w')
	f.write(r'#!/bin/bash\n')
	f.write(r'#' + '\n')
	f.write(r'#$ -cwd'+ '\n')
	f.write(r'#$ -j y'+ '\n')
	f.write(r'#$ -V'+ '\n')
	f.write(r'#$ -m es'+ '\n')
	f.write(r'#$ -o ll_out'+ '\n')
	f.write(r'#$ -S /bin/bash'+ '\n')
	f.write(r'#$ -N simulated_annealing'+str(dir)+' \n')
	f.write(' \n')
	f.write(r'./simulated_annealing.py' + '\n')
	f.close()
	f = os.popen('qsub anneal.sub')
	line = f.readline()
	f.close()
	print line
	os.chdir('..')
