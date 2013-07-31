#!/usr/bin/env python

import os

os.chdir('/home/zach/Genetic')

compiledEnergies = open('fullOutputGenetic.txt',mode = 'a')

subdirectories = os.listdir(os.getcwd())
subdirectories.remove('submit_hopping.py') #remove *????
for alphaValue in subdirectories[:94]:
  os.chdir(str(alphaValue))
  readFile = open('Output.txt', mode = 'r')
  for line in readFile:
    compiledEnergies.writelines(str(line))
  os.chdir('..')

compiledEnergies.close()
readFile.close()
