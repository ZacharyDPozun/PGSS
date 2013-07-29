#!/usr/bin/env python

import os

#Notes on file paths:
#first CD into ~zach/Basin_Hopping/Pt
#From there, the choices are [Ag, Au, Cu, Ir, Ni, Pd, Pt, Rh]
#For each here, the choices are [25, 50, 75]
#In here, the numbers change per folder (so get a list of the stuff in that folder)
#Then, get the line from best energies and draw out the second element, divide by number of atoms

os.chdir('/home/zach/Basin_Hopping/Pt')

compiledEnergies = open('bestEnergiesPerAlloy.txt','a')

elementChoices = ['Ag', 'Au', 'Cu', 'Ir', 'Ni', 'Pd', 'Rh']
percentPtChoies = [25, 50, 75]

for element in elementChoices:
  os.chdir(str(element))
  for percentage in percentPtChoies:
    if (element == 'Pt' and percentage != 50):
      break
    os.chdir(str(percentage))
    subdirectories = os.listdir(os.getcwd())
    subdirectories.remove('submit_hopping.py')
    for alloySize in subdirectories:
      #We want to return... Pt, Au, 160, 50, -353
      print str(alloySize), percentage, element
      os.chdir(str(alloySize))
      writeString = 'Pt ' + str(element) + ' ' + str(percentage)
      writeString += ' ' + str(alloySize)
      readFile = open('BestEnergies.txt', mode = 'r')
      for line in readFile:
        newLine = line.split()
      bestEnergy = float(newLine[1])
      bestEnergy /= float(alloySize)
      writeString += ' ' + str(bestEnergy) + '\n'
      compiledEnergies.write(writeString)
      os.chdir('..')
    os.chdir('..')
  os.chdir('..')
