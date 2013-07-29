#!/usr/bin/env python

import os

#Notes on file paths:
#first CD into ~zach/Basin_Hopping/Pt
#From there, the choices are [Ag, Au, Cu, Ir, Ni, Pd, Pt, Rh]
#For each here, the choices are [25, 50, 75]
#In here, the numbers change per folder (so get a list of the stuff in that folder)
#Then, get the line from best energies and draw out the second element, divide by number of atoms

os.chdir('/home/zach/Simulated_Annealing/Pt')

compiledEnergies = open('bestEnergiesPerAlloyAnnealing.txt','a')

elementChoices = ['Au', 'Ni', 'Pd','Pt']
percentPtChoies = [50]

for element in elementChoices:
  os.chdir(str(element))
  for percentage in percentPtChoies:
    if (element == 'Pt' and percentage != 50):
      break
    os.chdir(str(percentage))
    subdirectories = os.listdir(os.getcwd())
    try:
      subdirectories.remove('submit_annealing.py')
    except:
      1 == 1
    try:
      subdirectories.remove('POSCAR')
    except:
      1 == 1
    for alloySize in subdirectories:
      #We want to return... Pt, Au, 160, 50, -353
      print str(alloySize), percentage, element
      os.chdir(str(alloySize))
      writeString = 'Pt ' + str(element) + ' ' + str(percentage)
      writeString += ' ' + str(alloySize)
      readFile = open('EnergyList.txt', mode = 'r')
      line  = 'empty'
      while (line != ''):
        line = readFile.readline()
        if (line != ''): 
          newLine = line.split()
      print newLine
      bestEnergy = float(newLine[1])
      bestEnergy /= float(alloySize)
      writeString += ' ' + str(bestEnergy) + '\n'
      compiledEnergies.write(writeString)
      os.chdir('..')
    os.chdir('..')
  os.chdir('..')
