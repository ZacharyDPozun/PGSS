#!/usr/bin/env python

import os

compiledEnergies = open('fullOutputGenetic.txt',mode = 'a')

os.chdir('/home/zach/Genetic')

subdirectories = os.listdir(os.getcwd())
subdirectories.remove('submit_hopping.py') #remove *????
for alphaValue in subdirectories:
  os.chdir(str(alphaValue))
  readFile = open('Output.txt', mode = 'r')
  for line in readFile:
    compiledEnergies.write(str(line))
  os.chdir('..')

compiledEnergies.close()
readFile.close()
