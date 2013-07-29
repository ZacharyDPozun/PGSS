#!/usr/bin/env python

#################################################################
## Notes:  Format will be as follows on the line below:
## Pt SecondElement PercentPt AlloySize PEPerAtom
#################################################################

baseFile = open('bestEnergiesPerAlloy.txt', mode = 'r')
lines = []

for line in baseFile:
  group = line.split()
  lines.append(group)

ptPercent25 = [element for element in lines if element[2]=='25']

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