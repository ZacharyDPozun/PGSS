#!/usr/bin/env python

#
# GENETIC ALGORITHM
# SHELL BINDING (X) VS COST (Y)
#

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

  xValues = [float(element[4]) for element in lines if element[6]=='20']
  yValues = [float(element[3]) for element in lines if element[6]=='20']

  print "oh hell no"
  print xValues[0],yValues[0]
  print xValues[1],yValues[1]
  print len(yValues)

  plt.title("Strength vs. Cost")
  plt.scatter(xValues,yValues)
  plt.xlabel("Core Binding Strength")
  plt.ylabel("Cost ($/Gram)")
  plt.show()



makeSimpleGraph()
