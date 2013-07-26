#!/usr/bin/env python

import os

dirs = [i for i in range(10,201,2)]

for dir in dirs:
	os.mkdir(str(dir))
	os.chdir(str(dir)) 
	f = os.popen('cp ~/PGSS/returnFinalBestAtom.py ~/PGSS/basin_hopping.py ~/PGSS/commonfunctions.py .')
	f.close()
	f = open('run_hopping.py','w')
	f.write('#!/usr/bin/env python\n')
	atoms = dir / 2
	f.write(' \n')
	f.write('from basin_hopping import * \n')
	f.write(' \n')
	f.write("mainBasinLoop('Pt', 'Au', 78, 79, " + str(atoms) + ',' + str(atoms) + ", 6., 0.5)")
	f.write(' \n')
	f.close()
	f = open('hop.sub','w')
	f.write(r'#!/bin/bash\n')
	f.write(r'#' + '\n')
	f.write(r'#$ -cwd'+ '\n')
	f.write(r'#$ -j y'+ '\n')
	f.write(r'#$ -V'+ '\n')
	f.write(r'#$ -m es'+ '\n')
	f.write(r'#$ -o ll_out'+ '\n')
	f.write(r'#$ -S /bin/bash'+ '\n')
	f.write(r'#$ -N basin_hopping_Au'+str(dir)+' \n')
	f.write(' \n')
	f.write(r'./run_hopping.py' + '\n')
	f.close()
	f = os.popen('chmod 774 run_hopping.py')
	f.close() 
	f = os.popen('qsub hop.sub')
	line = f.readline()
	f.close()
	print line
	os.chdir('..')
