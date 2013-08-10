import math
import random
import csv
from numpy import *
from scipy import constants
from itertools import izip

########### user input ###########
#### XYZ file disorder
disorder_path = "./sample_data/disorder/si_h_14.xyz"
l = 2.351 # bond length in angstrom for Si
atoms = [2, 3, 4] # atoms you want to disorder
disorder_amount = [0.0, 0.0, 0.0] # amount you want to disorder each atom
repeat = 3 # number of times you want to repeat the disordering

#### Nonlinear reflection coefficient
nrc_path = "./sample_data/res" # full path to 'res' folder for desired structure
kpts = 1219 # kpoints of the case you want to use
ecut = 15 # ecut of the case you want to use
# Polarization
entry = "p" # 's' or 'p'
exit = "p" # 's' or 'p'
# Angles
theta_deg = 65 # angle of incidence in degrees
phi_deg = 30 # azimuthal angle
# Misc
n_0 = 50 # electronic density
########### end user input ###########

## Housekeeping
theta = math.radians(theta_deg)
phi = math.radians(phi_deg)

## Functions
def load_matrix(f):
	data = loadtxt(f)
	return data

def save_matrix(f, data):
	save = column_stack(data)
	savetxt(f, save)

def disorder(file):
	atoms[:] = [x - 1 for x in atoms]
	bl_bohr = l * constants.angstrom / constants.value("Bohr radius")
	xyz = load_matrix(file)
	data = xyz[atoms]
	final = zip(data, disorder_amount)
	for atom in range(0, len(final)):
		polar = random.random() * constants.pi
		azimuthal = random.random() * 2 * constants.pi
		cx = final[atom][0][0] + (final[atom][1] * (bl_bohr / 2) * math.sin(polar) * math.cos(azimuthal))
		cy = final[atom][0][1] + (final[atom][1] * (bl_bohr / 2) * math.sin(polar) * math.sin(azimuthal))
		cz = final[atom][0][2] + (final[atom][1] * (bl_bohr / 2) * math.cos(polar))
		print cx, cy, cz

def nonlinear_reflection_coefficient(polar_in, polar_out):
	energy = [row[0] for row in values("chi1")]
	chi1 = [float(row[1]) for row in values("chi1")]
	zzz = [float(row[1]) for row in values("zzz")]
	zxx = [float(row[1]) for row in values("zxx")]
	xxz = [float(row[1]) for row in values("xxz")]
	xxx = [float(row[1]) for row in values("xxx")]
	nrc = [[] for l in range(len(energy))]
	i = 0
	for energy, chi1, zzz, zxx, xxz, xxx in izip(energy, chi1, zzz, zxx, xxz, xxx):
		R = rif_constants(float(energy)) * math.fabs(fresnel(polar_out, "vs", float(energy), chi1) * fresnel(polar_out, "sb", float(energy), chi1) * ((fresnel(polar_in, "vs", float(energy), chi1) * fresnel(polar_in, "sb", float(energy), chi1)) ** 2) * r_factors(polar_in, polar_out, zzz, zxx, xxz, xxx, float(energy), chi1)) ** 2
		nrc[i].append(energy)
		nrc[i].append(R)
		i = i + 1
	return nrc

def values(component):
	if component == "chi1":
		path = file_path + "chi1.sm_xx_yy_zz_" + str(kpts) + "_" + str(ecut) + "-nospin_scissor_0_Nc_29"
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

def epsilon(chi1):
	epsilon = 1 + (4 * constants.pi * chi1)
	return epsilon

def fresnel(polarization, material, energy, chi1): # Remove plus ones!!!
	if polarization == "s" and material == "vs":
		t = (2 * math.cos(theta)) / (math.cos(theta) + wave_vector(energy, chi1))
	elif polarization == "s" and material == "sb":
		t = (2 * wave_vector(energy, chi1)) / (1 + wave_vector(energy, chi1) + wave_vector(energy, chi1))
	elif polarization == "p" and material == "vs":
		t = (2 * math.cos(theta)) / (epsilon(chi1) * math.cos(theta) + wave_vector(energy, chi1))
	elif polarization == "p" and material == "sb":
		t = (2 * wave_vector(energy, chi1)) / (1 + epsilon(chi1) * wave_vector(energy, chi1) + epsilon(chi1) * wave_vector(energy, chi1))
	return t

def r_factors(polar_in, polar_out, triperp, perpbipar, biparperp, tripar, energy, chi1):
	if polar_in == "p" and polar_out == "p":
		r = math.sin(theta) * epsilon(chi1) * (((math.sin(theta) ** 2) * (epsilon(chi1) ** 2) * triperp) + (wave_vector(energy, chi1) ** 2) * (epsilon(chi1) ** 2) * perpbipar) + epsilon(chi1) * epsilon(chi1) * wave_vector(energy, chi1) * wave_vector(energy, chi1) * (-2 * math.sin(theta) * epsilon(chi1) * biparperp + wave_vector(energy, chi1) * epsilon(chi1) * tripar * math.cos(3 * phi))
	elif polar_in == "s" and polar_out == "p":
		r = math.sin(theta) * epsilon(chi1) * perpbipar - wave_vector(energy, chi1) * epsilon(chi1) * tripar * math.cos(3 * phi)
	elif polar_in == "p" and polar_out == "s":
		r = -(wave_vector(energy, chi1) ** 2) * (epsilon(chi1) ** 2) * tripar * math.sin(3 * phi)
	elif polar_in == "s" and polar_out == "s":
		r = tripar * math.sin(3 * phi)
	return r

def wave_vector(energy, chi1): # the absolute value thing is just temporary
	k = (energy / constants.c) * math.sqrt(math.fabs(epsilon(chi1) - math.sin(theta) ** 2))
	return k

def output_stage():
	out_file = file_path + "R" + entry + exit + "_" + str(kpts) + "_" + str(ecut)
	with open(out_file, 'wb') as csvfile:
 		output = csv.writer(csvfile, delimiter=' ')
 		output.writerows(nonlinear_reflection_coefficient(entry, exit))

disorder(disorder_path)
#output_stage()
#print type(matrix_tests())
#print type(values("zzz")[0][0])


