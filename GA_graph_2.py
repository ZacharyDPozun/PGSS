#!/usr/bin/env python

#
# GENETIC ALGORITHM
# ALPHA (X) vs %COMPOSITION OF ELEMENT
# FOUR PLOTS:  PT, AG, NI, CU

LENGTH = 93

from matplotlib import pyplot
import pylab
from mpl_toolkits.mplot3d import Axes3D
import numpy
import matplotlib.pyplot as plt

def makeSimpleGraph():
  baseFile = open('betterFilev3.txt',mode = 'r')
  lines = []

  for line in baseFile:
    group = line.split()
    lines.append(group)


  lines = [line for line in lines if len(line)==8]

  platinumPairs = [((line[0]=='Pt')*int(line[2])+(line[1]=='Pt')*(85-int(line[2])),int(line[7])) for line in lines if line[6]=='20']
  silverPairs = [((line[0]=='Ag')*int(line[2])+(line[1]=='Ag')*(85-int(line[2])),int(line[7])) for line in lines if line[6]=='20']
  nickelPairs = [((line[0]=='Ni')*int(line[2])+(line[1]=='Ni')*(85-int(line[2])),int(line[7])) for line in lines if line[6]=='20']
  copperPairs = [((line[0]=='Cu')*int(line[2])+(line[1]=='Cu')*(85-int(line[2])),int(line[7])) for line in lines if line[6]=='20']

  platY = []
  silverY = []
  nickelY = []
  copperY = []


  for k in range(LENGTH):
    platY.append(sum( piece[0] for piece in platinumPairs if piece[1]==k )/78.75 )
    silverY.append(sum( piece[0] for piece in silverPairs if piece[1]==k )/78.75 )
    nickelY.append(sum( piece[0] for piece in nickelPairs if piece[1]==k )/78.75 )
    copperY.append(sum( piece[0] for piece in copperPairs if piece[1]==k )/78.75 )

  xRange = range(LENGTH)

  plt.title("Some stuff")
  plt.plot(xRange,platY,color='b',label = "Pt")
  plt.plot(xRange,silverY,color='r',label= "Ag")
  plt.plot(xRange,nickelY,color='g',label="Ni")
  plt.plot(xRange,copperY,color='y',label="Cu")
  plt.legend()
  plt.show()

makeSimpleGraph()
