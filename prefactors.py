from math import sqrt, radians, sin, cos
from scipy import constants
from numpy import linspace, column_stack

energy = linspace(0.01, 12.00, 1200)

Ab = 0.53e-8 # Angstroms
theta = radians(65)
phi = radians(30)

#cpf = -ci*(Ab/a0)**2/(2.*Area)
#cpf = cpf * (2.*Ry)**2

def electrostatic_units(energy): 
	complex_esu = 1j * ((2 * constants.value("Rydberg constant times hc in eV")) ** 5) * ((Ab / (constants.value("lattice parameter of silicon") * 100)) ** 5) / ((2 * sqrt(3)) / ((2 * sqrt(2)) ** 2))
	factor = (complex_esu * 2.08e-15 * (((constants.value("lattice parameter of silicon") * 100) / 1e-8) ** 3)) / (energy ** 3)
	return factor

#pfesu = (32 * (constants.pi ** 3)) / (1e-7 * 1e-21 * ((constants.c * 100) ** 3) * (cos(theta) ** 2) * (constants.value("Planck constant over 2 pi in eV s") ** 2))

print electrostatic_units(energy)