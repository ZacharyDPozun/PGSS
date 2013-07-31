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
  baseFile = open('betterFile3.txt',mode = 'r')
  lines = []

  for line in baseFile:
    group = line.split()
    lines.append(group)


  lines = [line for line in lines if len(line)==8]

  xValues = [abs(float(element[4])-0.562537) for element in lines if element[6]=='20']
  yValues = [abs(float(element[3])-28.18966) for element in lines if element[6]=='20']
  xValues2 = [abs(float(element[4])-0.562537) for element in lines if element[6]!='20']
  yValues2 = [abs(float(element[3])-28.18966) for element in lines if element[6]!='20']


  print xValues[0],yValues[0]
  print xValues[1],yValues[1]
  print len(yValues)

  plt.title("Strength vs. Cost")
  plt.scatter(xValues2,yValues2,color='r')
  plt.scatter(xValues,yValues,color='b')
  plt.xlabel("Core Binding Strength")
  plt.ylabel("Cost ($/Gram)")
  plt.show()



makeSimpleGraph()
