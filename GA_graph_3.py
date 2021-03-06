#!/usr/bin/env python

#
# GENETIC ALGORITHM
# ALPHA (X) vs %COMPOSITION OF ELEMENT
# FOUR PLOTS:  PT, AG, NI, CU

LENGTH = 101

from matplotlib import pyplot
import pylab
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import spline

def makeSimpleGraph():
  baseFile = open('betterFile3.txt',mode = 'r')
  lines = []

  for line in baseFile:
    group = line.split()
    lines.append(group)


  lines = [line for line in lines if len(line)==8]

  gen = str(20)

  platinumPairs = [((line[0]=='Pt')*int(line[2])+(line[1]=='Pt')*(85-int(line[2])),int(line[7])) for line in lines if line[6]==gen]
  silverPairs = [((line[0]=='Ag')*int(line[2])+(line[1]=='Ag')*(85-int(line[2])),int(line[7])) for line in lines if line[6]==gen]
  nickelPairs = [((line[0]=='Ni')*int(line[2])+(line[1]=='Ni')*(85-int(line[2])),int(line[7])) for line in lines if line[6]==gen]
  copperPairs = [((line[0]=='Cu')*int(line[2])+(line[1]=='Cu')*(85-int(line[2])),int(line[7])) for line in lines if line[6]==gen]

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

  xNew = np.linspace(min(xRange),max(xRange),20)
  power_smooth1 = spline(xRange,platY,xNew)
  power_smooth2 = spline(xRange,silverY,xNew)
  power_smooth3 = spline(xRange,nickelY,xNew)
  power_smooth4 = spline(xRange,copperY,xNew)


  plt.title("% Composition per Alpha Value")

  #plt.plot(xRange,platY,color='b',label = "Pt")
  #plt.plot(xRange,silverY,color='r',label= "Ag")
  #plt.plot(xRange,nickelY,color='g',label="Ni")
  #plt.plot(xRange,copperY,color='y',label="Cu")
  plt.plot(xNew,power_smooth1,color='b',label="Pt")
  plt.plot(xNew,power_smooth2,color='r',label="Ag")
  plt.plot(xNew,power_smooth3,color='g',label="Ni")
  plt.plot(xNew,power_smooth4,color='y',label="Cu")
  plt.xlabel("Alpha Value")
  plt.ylabel("Percent Composition")
  plt.legend(loc="center left",bbox_to_anchor=(1,0.5))
  plt.show()

makeSimpleGraph()
