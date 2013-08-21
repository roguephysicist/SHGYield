from math import sqrt, radians, sin, cos
from scipy import constants

hbar = 1.0545e-27 # erg s
sc = 2.997925e10 # cm/s
Ab = 0.53e-8 # Angstroms
a0 = 5.34e-8 # A for Si
Ry = 13.606 # eV Rydberg
Area = (1 / (2 * sqrt(2)) ** 2) * 2 * sqrt(3)
theta = radians(65)
phi = radians(30)
#gamma = (1e-40) * (1e7) * (1.6021e-12) / ((hbar ** 2) * (sc ** 3))
#gamma = gamma * 32 * (constants.pi ** 3)  # 32pi^3 \approx 1000
#gamma = gamma * (a0 ** 2) * (Ab ** 3) * ((2 * Ry) ** 3)
#cpf = -1j * (Ab / a0) ** 2 / (2 * Area)
#cpf = cpf * (2 * Ry) ** 2

fv1 = (2 * Ry) ** 9 / (a0 / Ab) ** 4
fv1 = fv1 * 1e7 * 1.6021e-12
fv2 = Ab ** 5 / ((hbar ** 2) * (sc ** 3))
gamma = fv1 * fv2
gamma = (gamma * 32 * (constants.pi ** 3) / Area ** 2) # cm^2/Watt

hbeV = 6.582e-16 # hbar in eV.s
gammaesu = (32 * (constants.pi ** 3) / ((hbeV ** 2) * (sc ** 3)))
cpf = 1j * (2 * Ry) ** 5 * (Ab / a0) ** 5 / Area # line 68
qz = cos(theta)
Q  = sin(theta)
sect2 = 1 / (qz ** 2)
pf = gamma * sect2
pfesu = gammaesu * sect2 * 1.e7
pfesu = pfesu * 1e21 # R in 10^{-21} cm^2/W

#line 112

print gamma, gammaesu
print cpf
print qz, Q
print pf, pfesu