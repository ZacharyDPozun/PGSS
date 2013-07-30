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

def basinVsAnnealGraphs(inElement="Au",make3DGraph=False):

  BHFile = open('bestEnergiesPerAlloy.txt',mode='r')
  BHlines = []


  SAFile = open('bestEnergiesPerAlloyAnnealing.txt',mode='r')
  SAlines = []

  for line in BHFile:
    group = line.split()
    BHlines.append(group)

  for line in SAFile:
    group = line.split()
    SAlines.append(group)

  ptPercent50BH = [element for element in BHlines if element[2]=='50']
  ptPercent50SA = [element for element in SAlines if element[2]=='50']

  Au50BH = [element for element in ptPercent50BH if element[1]=='Au']
  Ni50BH = [element for element in ptPercent50BH if element[1]=='Ni']
  Pd50BH = [element for element in ptPercent50BH if element[1]=='Pd']
  Au50SA = [element for element in ptPercent50SA if element[1]=='Au']
  Ni50SA = [element for element in ptPercent50SA if element[1]=='Ni']
  Pd50SA = [element for element in ptPercent50SA if element[1]=='Pd']

  Au50BH.sort(key=lambda x: int(x[3]))
  Ni50BH.sort(key=lambda x: int(x[3]))
  Pd50BH.sort(key=lambda x: int(x[3]))
  Au50SA.sort(key=lambda x: int(x[3]))
  Ni50SA.sort(key=lambda x: int(x[3]))
  Pd50SA.sort(key=lambda x: int(x[3]))

  x1BH,x2BH,x3BH,x1SA,x2SA,x3SA = [],[],[],[],[],[]
  y1BH,y2BH,y3BH,y1SA,y2SA,y3SA = [],[],[],[],[],[]
  z1BH,z2BH,z3BH,z1SA,z2SA,z3SA = [],[],[],[],[],[]

  for k in range(len(Au50BH)):
    x1BH.append(int(Au50BH[k][3]))
    x2BH.append(int(Ni50BH[k][3]))
    x3BH.append(int(Pd50BH[k][3]))
    x1SA.append(int(Au50SA[k][3]))
    x2SA.append(int(Ni50SA[k][3]))
    x3SA.append(int(Pd50SA[k][3]))
    y1BH.append(0)
    y2BH.append(0)
    y3BH.append(0)
    y1SA.append(1)
    y2SA.append(1)
    y3SA.append(1)
    z1BH.append(float(Au50BH[k][4]))
    z2BH.append(float(Ni50BH[k][4]))
    z3BH.append(float(Pd50BH[k][4]))
    z1SA.append(float(Au50SA[k][4]))
    z2SA.append(float(Ni50SA[k][4]))
    z3SA.append(float(Pd50SA[k][4]))

  if make3DGraph:
    fig = pylab.figure()
    ax = Axes3D(fig)
    if inElement == "Au":
      ax.plot(x1BH,y1BH,z1BH,label="BH")
      ax.plot(x1SA,y1SA,z1SA,label="SA")
      pyplot.title("Gold alloys", fontsize=12)
    elif inElement == "Ni":
      ax.plot(x2BH,y2BH,z2BH,label="BH")
      ax.plot(x2SA,y2SA,z2SA,label="SA")
      pyplot.title("Nickel alloys", fontsize=12)
    elif inElement == "Pd":
      ax.plot(x3BH,y3BH,z3BH,label="BH")
      ax.plot(x3SA,y3SA,z3SA,label="SA")
      pyplot.title("Palladium alloys", fontsize=12)
    else:
      pass
    ax.legend(prop={'size':9})
    pylab.yticks([0,1],["BH","SA"])
    ax.set_xlabel("# of Atoms")
    ax.set_ylabel("BH or SA")
    ax.set_zlabel('PE per Atom')
    pyplot.show()
  else:
  #make everything an array
  #keep track of what's what
  #plot and give options
    x1BHa = np.array(x1BH)
    y1BHa = np.array(z1BH)
    x2BHa = np.array(x2BH)
    y2BHa = np.array(z2BH)
    x3BHa = np.array(x3BH)
    y3BHa = np.array(z3BH)
    x1SAa = np.array(x1SA)
    y1SAa = np.array(z1SA)
    x2SAa = np.array(x2SA)
    y2SAa = np.array(z2SA)
    x3SAa = np.array(x3SA)
    y3SAa = np.array(z3SA)

    if inElement=="Au":
      plt.title("Gold alloys")
      plt.plot(x1BHa,y1BHa,label = "Basin")
      plt.plot(x1SAa,y1SAa,label = "Anneal")
    if inElement=="Ni":
      plt.title("Nickel alloys")
      plt.plot(x1BHa,y1BHa,label = "Basin")
      plt.plot(x1SAa,y1SAa,label = "Anneal")
    if inElement=="Pd":
      plt.title("Palladium alloys")
      plt.plot(x1BHa,y1BHa,label = "Basin")
      plt.plot(x1SAa,y1SAa,label = "Anneal")
    plt.legend()
    plt.show()

basinVsAnnealGraphs("Pd")