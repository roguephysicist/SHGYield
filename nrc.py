import math
import csv
from itertools import izip
import numpy as np
from scipy import constants

## User Input
# Input
file_path = "./sample_data/" # full path to 'res' folder for desired structure
kpts = 514 # kpoints of the case you want to use
ecut = 10 # ecut of the case you want to use
# Polarization
entry = "s" # 's' or 'p'
exit = "s" # 's' or 'p'
# Angles
theta_deg = 65 # angle of incidence in degrees
phi_deg = 30 # azimuthal angle

## Housekeeping
n_0 = constants.N_A # electronic density - for testing purposes
theta = math.radians(theta_deg)
phi = math.radians(phi_deg)

## Functions
def nonlinear_reflection_coefficient(polar_in, polar_out):
	energy = [row[0] for row in values("epsilon")]
	epsilon = [float(row[1]) for row in values("epsilon")]
	zzz = [float(row[0]) for row in values("zzz")]
	zxx = [float(row[0]) for row in values("zxx")]
	xxz = [float(row[0]) for row in values("xxz")]
	xxx = [float(row[0]) for row in values("xxx")]
	nrc = [[] for l in range(len(energy))]
	i = 0
	for energy, epsilon, zzz, zxx, xxz, xxx in izip(energy, epsilon, zzz, zxx, xxz, xxx):
		R = rif_constants(float(energy)) * math.fabs(fresnel(polar_out, "vs", float(energy), epsilon) * fresnel(polar_out, "sb", float(energy), epsilon) * ((fresnel(polar_in, "vs", float(energy), epsilon) * fresnel(polar_in, "sb", float(energy), epsilon)) ** 2) * r_factors(polar_in, polar_out, zzz, zxx, xxz, xxx, float(energy), epsilon)) ** 2
		nrc[i].append(energy)
		nrc[i].append(R)
		i = i + 1
	return nrc

def matrix_tests():
	matrix = np.array(values("epsilon"), float)
	return matrix

def values(component):
	if component == "epsilon":
		path = file_path + "chi1.sm_xx_" + str(kpts) + "_" + str(ecut) + "-nospin_scissor_0_Nc_29"
		with open(path, 'rb') as csvfile:
			data = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
			value = [[row[0], row[1]] for row in data]
	elif component == "zzz" or component == "zxx" or component == "xxz" or component == "xxx":
		path = file_path + "shgC.sm_" + component + "_" + str(kpts) + "_half-slab_" + str(ecut) + "-nospin_scissor_0_Nc_29"
		with open(path, 'rb') as csvfile:
			data = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
			value = [[row[0], row[3]] for row in data]
	return value

def rif_constants(energy):
	const = (32 * constants.pi ** 3 * energy ** 2) / (n_0 * constants.e ** 2 * constants.c ** 3 * math.cos(theta) ** 2)
	return const

def fresnel(polarization, material, energy, epsilon): # REMOVE PLUS ONE FROM LINE 69 ASAP
	if polarization == "s" and material == "vs":
		t = (2 * math.cos(theta)) / (math.cos(theta) + wave_vector(energy, epsilon))
	elif polarization == "s" and material == "sb":
		t = (2 * wave_vector(energy, epsilon)) / (1 + wave_vector(energy, epsilon) + wave_vector(energy, epsilon))
	elif polarization == "p" and material == "vs":
		t = (2 * math.cos(theta)) / (epsilon * math.cos(theta) + wave_vector(energy, epsilon))
	elif polarization == "p" and material == "sb":
		t = (2 * wave_vector(energy, epsilon)) / (epsilon * wave_vector(energy, epsilon) + epsilon * wave_vector(energy, epsilon))
	return t

def r_factors(polar_in, polar_out, triperp, perpbipar, biparperp, tripar, energy, epsilon):
	if polar_in == "p" and polar_out == "p":
		r = math.sin(theta) * epsilon * (((math.sin(theta) ** 2) * (epsilon ** 2) * triperp) + (wave_vector(energy, epsilon) ** 2) * (epsilon ** 2) * perpbipar) + epsilon * epsilon * wave_vector(energy, epsilon) * wave_vector(energy, epsilon) * (-2 * math.sin(theta) * epsilon * biparperp + wave_vector(energy, epsilon) * epsilon * tripar * math.cos(3 * phi))
	elif polar_in == "s" and polar_out == "p":
		r = math.sin(theta) * epsilon * perpbipar - wave_vector(energy, epsilon) * epsilon * tripar * math.cos(3 * phi)
	elif polar_in == "p" and polar_out == "s":
		r = -(wave_vector(energy, epsilon) ** 2) * (epsilon ** 2) * tripar * math.sin(3 * phi)
	elif polar_in == "s" and polar_out == "s":
		r = tripar * math.sin(3 * phi)
	return r

def wave_vector(energy, epsilon): # the absolute value thing is just temporary
	k = (energy / constants.c) * math.sqrt(math.fabs(epsilon - math.sin(theta) ** 2))
	return k

def output_stage():
	out_file = file_path + "R" + entry + exit + "_" + str(kpts) + "_" + str(ecut)
	with open(out_file, 'wb') as csvfile:
 		output = csv.writer(csvfile, delimiter=' ')
 		output.writerows(nonlinear_reflection_coefficient(entry, exit))

output_stage()
#print type(matrix_tests())
#print type(values("zzz")[0][0])


