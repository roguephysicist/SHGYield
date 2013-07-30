import math
from scipy import constants # physical constants (Scipy will need to be installed on medusa)

## User Input
# Files
file_path = 0 # full path to 'res' folder for desired structure
kpts = 0 # kpoints of the case you want to use
ecut = 0 # ecut of the case you want to use
# Polarization
entry = "p" # 's' or 'p'
exit = "p" # 's' or 'p'
# Angles
theta_deg = 65 # angle of incidence in degrees
phi_deg = 30 # azimuthal angle

## Constants
pi = constants.pi
e = constants.e
c = constants.c
n_0 = constants.N_A # electronic density - for testing purposes

## Responses
chi_3perp = 0 # zzz
chi_perp2par = 0 # zxx
chi_2parperp = 0 # xxz
chi_3par = 0 # xxx

## Housekeeping
theta = math.radians(theta_deg)
phi = math.radians(phi_deg)
omega = 2 # for testing purposes
epsilon = 1 # for testing purposes

## Functions
def rif_constants():
	const = (32 * pi ** 3 * omega ** 2) / (n_0 * e ** 2 * c ** 3 * math.cos(theta) ** 2)
	return const

def wave_vector(layer):
	if layer == "s":
		k = (omega / c) * math.sqrt(epsilon - math.sin(theta) ** 2)
	elif layer == "b":
		k = (omega / c) * math.sqrt(epsilon - math.sin(theta) ** 2)
	return k

def fresnel(polarization, material):
	if polarization == "s" and material == "vs":
		t = (2 * math.cos(theta)) / (math.cos(theta) + wave_vector("s"))
	elif polarization == "s" and material == "sb":
		t = (2 * wave_vector("s")) / (wave_vector("s") + wave_vector("b"))
	elif polarization == "p" and material == "vs":
		t = (2 * math.cos(theta)) / (epsilon * math.cos(theta) + wave_vector("s"))
	elif polarization == "p" and material == "sb":
		t = (2 * wave_vector("s")) / (epsilon * wave_vector("s") + epsilon * wave_vector("b"))
	return t

def r_factors():
	if entry == "p" and exit == "p":
		r = math.sin(theta) * epsilon * (((math.sin(theta) ** 2) * (epsilon ** 2) * chi_3perp) + (wave_vector("b") ** 2) * (epsilon ** 2) * chi_perp2par) + epsilon * epsilon * wave_vector("b") * wave_vector("b") * (-2 * math.sin(theta) * epsilon * chi_2parperp + wave_vector("b") * epsilon * chi_3par * math.cos(3 * phi))
	elif entry == "s" and exit == "p":
		r = math.sin(theta) * epsilon * chi_perp2par - wave_vector("b") * epsilon * chi_3par * math.cos(3 * phi)
	elif entry == "p" and exit == "s":
		r = -(wave_vector("b") ** 2) * (epsilon ** 2) * chi_3par * math.sin(3 * phi)
	elif entry == "s" and exit == "s":
		r = chi_3par * math.sin(3 * phi)
	return r

def rif_all():
 	Rif = rif_constants() * math.fabs(fresnel(exit, "vs") * fresnel(exit, "sb") * ((fresnel(entry, "vs") * fresnel(entry, "sb")) ** 2) * r_factors()) ** 2
 	return Rif

print rif_all()
