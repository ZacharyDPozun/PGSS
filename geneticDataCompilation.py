#!/usr/bin/env python

import os

#Notes on file paths:
#first CD into ~zach/Basin_Hopping/Pt
#From there, the choices are [Ag, Au, Cu, Ir, Ni, Pd, Pt, Rh]
#For each here, the choices are [25, 50, 75]
#In here, the numbers change per folder (so get a list of the stuff in that folder)
#Then, get the line from best energies and draw out the second element, divide by number of atoms

os.chdir('/home/zach/Genetic')

compiledEnergies = open('fullOutputGenetic.txt','a')

subdirectories = os.listdir(os.getcwd())
subdirectories.remove('submit_hopping.py*') #remove *????
for alphaValue in subdirectories:
  os.chdir(str(alphaValue))
  readFile = open('Output.txt', mode = 'r')
  for line in readFile:
    compiledEnergies.writelines(str(line))
  os.chdir('..')

  compiledEnergies.close()
  readFile.close()
