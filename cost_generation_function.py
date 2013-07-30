#!/usr/bin/env python

def stupidCostFunction(intakeList):
  """Takes in a list from our text files
  and returns the price per unit of electronegativity
  per atom.  List should be in the form demonstrated
  in the code below
  [Pt, element2, percentPt, #Atoms, PE]"""
  element1 = "Pt"
  element2 = str(intakeList[2])
  percentPt = int(intakeList[3])
  nAtoms = int(intakeList[4])
  percent2 = 100 - percentPt
  potentialPerAtom = float(intakeList[5])

  #element prices per gram
  Rh = 76.52
  Ir = 23.00
  Ag = 0.62308
  Au = 41.25
  Cu = 0.00661
  Ni = 0.01543
  Pd = 24.25085
  Pt = 45.3008


  if str(element2) == "Rh":
    price = percentPt*Pt + percent2*Rh
  elif str(element2) == "Ir":
    price = percentPt*Pt + percent2*Ir
  elif str(element2) == "Ag":
    price = percentPt*Pt + percent2*Ag
  elif str(element2) == "Au":
    price = percentPt*Pt + percent2*Au
  elif str(element2) == "Cu":
    price = percentPt*Pt + percent2*Cu
  elif str(element2) == "Ni":
    price = percentPt*Pt + percent2*Ni
  else:
    price = percentPt*Pt + percent2*Pd

  return -1*price/potentialPerAtom


fullFile = open('bestEnergiesPerAlloyAnnealing.txt', mode = 'r')
newFile = open('SABestEnergiesWithPAPECost.txt', mode = 'a')

for line in fullFile:
  newList = line.split()
  newList.append(str(stupidCostFunction(newList)))
  newFile.write(str(newList))
  newFile.write('\n')

fullFile.close()
newFile.close()

