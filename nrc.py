
import math
from scipy import constants # physical constants (Scipy will need to be installed on medusa)

theta = 65 # angle of incidence in degrees, math.radians can convert (doc 9.2.4)
phi = 30 # azimuthal angle
omega = 2 # for testing purposes
epsilon = 1
n_0 = constants.N_A # electronic density - for testing purposes

radthet = math.radians(theta)
radphi = math.radians(phi)

def rif_constants(theta, omega):
	const = (32 * (constants.pi ** 3) * (omega ** 2)) / (((n_0 * constants.e) ** 2) * (constants.c ** 3) * (math.cos(radthet)) ** 2)
	return const

def wave_vector(omega, epsilon):
	kzj = (omega / constants.c) * math.sqrt(epsilon - (math.sin(radthet) ** 2))
	return kzj

# print rif_(theta, omega)