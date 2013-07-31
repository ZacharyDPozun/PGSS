
import os

def readGeneticOutputText():

  f = open("fullOutputGenetic.txt", mode = "r")
  parent = os.path.split(os.getcwd())[-1]

  lines = f.readlines()

  for line in lines:
    lines[lines.index(line)] = line.split()

  for line in lines:
    if len(line)==2:
      generation = int(line[1])
    line.append(generation)
    #line.append(str(parent))

  newLines = []

  newLines = [line for line in lines if (line[0]!='Average') and len(line)==7]
  generationAverages = []

  currentGeneration = 0
  currentAlpha = -1

  for k in range(len(newLines)):
    if k%38==36:
      currentGeneration+=1
    if k%735==0:
      currentAlpha+=1
    newLines[k][2] = int(newLines[k][2])
    newLines[k][3] = float(newLines[k][3])
    newLines[k][4] = float(newLines[k][4])
    newLines[k][5] = float(newLines[k][5])
    #newLines.append(currentGeneration)
    newLines[k].append(currentAlpha)

  currentGeneration = 0
  currentAlpha = 0

  for k in range(len(lines)):
    if k%38==36:
      line = [float(lines[k][2][:-1]),float(lines[k][6][:-1])]
      line.append(currentGeneration%21)
      currentGeneration += 1
      line.append(k/834)
      #line.append(int(parent))
      generationAverages += [line]

  generationAverages.remove(generationAverages[-1])

  finalList = newLines + generationAverages

  return finalList

x = readGeneticOutputText()

betterFile = open('betterFile.txt', mode = 'a')

for line in x:
  betterFile.writelines(str(line))
  betterFile.write('\n')

betterFile.close()

print "done"