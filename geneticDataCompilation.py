#!/usr/bin/env python

import os

compiledEnergies = open('fullOutputGenetic.txt',mode = 'a')

os.chdir('/home/zach/Genetic')

subdirectories = os.listdir(os.getcwd())
subdirectories.remove('submit_hopping.py') #remove *????
for alphaValue in subdirectories[:94]:
  os.chdir(str(alphaValue))
  readFile = open('Output.txt', mode = 'r')
  dataLines = [line for line in readFile if len(line)==6]
  averageLines = [line for line in readFile if len(line)==7]
  allLines = dataLines+averageLines
  for line in allLines:
    line.append(alphaValue)
  for line in allLines:
    compiledEnergies.writelines(str(line))
  os.chdir('..')

compiledEnergies.close()
readFile.close()
