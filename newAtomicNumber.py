def newAtomicNumber(atoms):
	numbers=atoms.get_atomic_numbers()
	random.shuffle(numbers)
	atoms.set_atomic_numbers(numbers)
	return atoms.get_atomic_numbers()

newAtomicNumber(atoms)
