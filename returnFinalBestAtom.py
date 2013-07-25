
def returnFinalBestAtom(trajectoryFile):
  candidates = PickleTrajectory(str(trajectoryFile),mode='r')
  molecules = [molecule for molecule in candidates]
  potentialEnergyList = [molecule.get_potential_energy() for molecule in molecules]
  minimumPotenialEnergy = min(potentialEnergyList)
  indexOfMinimumPotentialEnergy = potentialEnergyList.index(minimumPotenialEnergy)
  bestAtom = molecules[indexOfMinimumPotentialEnergy].copy()
  print "Yes We Can"
  return bestAtom
