#!/usr/bin/env python

#################################################################
## Notes:  Format will be as follows on the line below:
## Pt SecondElement PercentPt AlloySize PEPerAtom
#################################################################

from matplotlib import pyplot
import pylab
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def makeGraphs(makePictures=False):

  baseFile = open('BHBestEnergiesWithPAPECost.txt', mode = 'r')
  lines = []

  for line in baseFile:
    group = line.split()
    lines.append(group)

  ptPercent25 = [element for element in lines if element[2]=='75']

  Pt25Ag75 = [element for element in ptPercent25 if element[1]=='Ag']
  Pt25Au75 = [element for element in ptPercent25 if element[1]=='Au']
  Pt25Cu75 = [element for element in ptPercent25 if element[1]=='Cu']
  Pt25Ir75 = [element for element in ptPercent25 if element[1]=='Ir']
  Pt25Ni75 = [element for element in ptPercent25 if element[1]=='Ni']
  Pt25Pd75 = [element for element in ptPercent25 if element[1]=='Pd']
  Pt25Rh75 = [element for element in ptPercent25 if element[1]=='Rh']

  Pt25Ag75.sort(key=lambda x: int(x[3]))
  Pt25Au75.sort(key=lambda x: int(x[3]))
  Pt25Cu75.sort(key=lambda x: int(x[3]))
  Pt25Ir75.sort(key=lambda x: int(x[3]))
  Pt25Ni75.sort(key=lambda x: int(x[3]))
  Pt25Pd75.sort(key=lambda x: int(x[3]))
  Pt25Rh75.sort(key=lambda x: int(x[3]))

  Ag75xyz, Au75xyz, Cu75xyz, Ir75xyz, Ni75xyz, Pd75xyz, Rh75xyz = [],[],[],[],[],[],[]

  for element in range(len(Pt25Ag75)):
    Ag75xyz.append([int(Pt25Ag75[element][3]),Pt25Ag75[element][1],float(Pt25Ag75[element][5])])
    Au75xyz.append([int(Pt25Au75[element][3]),Pt25Au75[element][1],float(Pt25Au75[element][5])])
    Cu75xyz.append([int(Pt25Cu75[element][3]),Pt25Cu75[element][1],float(Pt25Cu75[element][5])])
    Ir75xyz.append([int(Pt25Ir75[element][3]),Pt25Ir75[element][1],float(Pt25Ir75[element][5])])
    Ni75xyz.append([int(Pt25Ni75[element][3]),Pt25Ni75[element][1],float(Pt25Ni75[element][5])])
    Pd75xyz.append([int(Pt25Pd75[element][3]),Pt25Pd75[element][1],float(Pt25Pd75[element][5])])
    Rh75xyz.append([int(Pt25Rh75[element][3]),Pt25Rh75[element][1],float(Pt25Rh75[element][5])])

  for element in range(len(Ag75xyz)):
    Ag75xyz[element][1] = 0
    Au75xyz[element][1] = 1
    Cu75xyz[element][1] = 2
    Ir75xyz[element][1] = 3
    Ni75xyz[element][1] = 4
    Pd75xyz[element][1] = 5
    Rh75xyz[element][1] = 6

  x1,x2,x3,x4,x5,x6,x7 = [],[],[],[],[],[],[]
  y1,y2,y3,y4,y5,y6,y7 = [],[],[],[],[],[],[]
  z1,z2,z3,z4,z5,z6,z7 = [],[],[],[],[],[],[]

  for element in range(len(Ag75xyz)):
    x1.append(Ag75xyz[element][0])
    y1.append(Ag75xyz[element][1])
    z1.append(Ag75xyz[element][2])
    x2.append(Au75xyz[element][0])
    y2.append(Au75xyz[element][1])
    z2.append(Au75xyz[element][2])
    x3.append(Cu75xyz[element][0])
    y3.append(Cu75xyz[element][1])
    z3.append(Cu75xyz[element][2])
    x4.append(Ir75xyz[element][0])
    y4.append(Ir75xyz[element][1])
    z4.append(Ir75xyz[element][2])
    x5.append(Ni75xyz[element][0])
    y5.append(Ni75xyz[element][1])
    z5.append(Ni75xyz[element][2])
    x6.append(Pd75xyz[element][0])
    y6.append(Pd75xyz[element][1])
    z6.append(Pd75xyz[element][2])
    x7.append(Rh75xyz[element][0])
    y7.append(Rh75xyz[element][1])
    z7.append(Rh75xyz[element][2])

  plt.title("Basin Hopping with Cost: 75% Platinum")
  plt.plot(x1,z1,label="Silver")
  plt.plot(x2,z2,label="Gold")
  plt.plot(x3,z3,label="Copper")
  plt.plot(x4,z4,label="Iridium")
  plt.plot(x5,z5,label="Nickel")
  plt.plot(x6,z6,label="Palladium")
  plt.plot(x7,z7,label="Rhodium")
  plt.xlabel('# of Atoms')
  plt.ylabel('Cost/(PE per Atom)')
  #plt.legend()
  plt.show()
  plt.autoscale_view()

makeGraphs(False)
