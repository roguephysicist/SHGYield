import math
import csv
from scipy import constants # physical constants (Scipy will need to be installed on medusa)

## User Input
# Files
file_path = "/Users/sma/Downloads/res/" # full path to 'res' folder for desired structure
kpts = 514 # kpoints of the case you want to use
ecut = 25 # ecut of the case you want to use
# Polarization
entry = "p" # 's' or 'p'
exit = "p" # 's' or 'p'
# Angles
theta_deg = 65 # angle of incidence in degrees
phi_deg = 30 # azimuthal angle

## Housekeeping
n_0 = constants.N_A # electronic density - for testing purposes
theta = math.radians(theta_deg)
phi = math.radians(phi_deg)
omega = 1
epsilon = 1


## Functions
def csv_import(file):
	with open(file, 'rb') as csvfile:
			reader = csv.reader(csvfile, delimiter=' ', quoting=csv.QUOTE_NONE, skipinitialspace=True)
			for row in reader:
				comp_list = row[3]
			return comp_list

def chi(component):
		path = file_path + "shgC.sm_" + component + "_" + str(kpts) + "_half-slab_" + str(ecut) + "-nospin_scissor_0_Nc_29"
		chi = csv_import(path)
		#chi = 1
		return chi

def rif_constants():
	const = (32 * constants.pi ** 3 * omega ** 2) / (n_0 * constants.e ** 2 * constants.c ** 3 * math.cos(theta) ** 2)
	return const

def wave_vector(layer):
	if layer == "s":
		k = (omega / constants.c) * math.sqrt(epsilon - math.sin(theta) ** 2)
	elif layer == "b":
		k = (omega / constants.c) * math.sqrt(epsilon - math.sin(theta) ** 2)
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

def r_factors(p_entry, p_exit):
	if p_entry == "p" and p_exit == "p":
		r = math.sin(theta) * epsilon * (((math.sin(theta) ** 2) * (epsilon ** 2) * chi("zzz")) + (wave_vector("b") ** 2) * (epsilon ** 2) * chi("zxx")) + epsilon * epsilon * wave_vector("b") * wave_vector("b") * (-2 * math.sin(theta) * epsilon * chi("xxz") + wave_vector("b") * epsilon * chi("xxx") * math.cos(3 * phi))
	elif p_entry == "s" and p_exit == "p":
		r = math.sin(theta) * epsilon * chi("zxx") - wave_vector("b") * epsilon * chi("xxx") * math.cos(3 * phi)
	elif p_entry == "p" and p_exit == "s":
		r = -(wave_vector("b") ** 2) * (epsilon ** 2) * chi("xxx") * math.sin(3 * phi)
	elif p_entry == "s" and p_exit == "s":
		r = chi("xxx") * math.sin(3 * phi)
	return r

def nonlinear_reflection_coefficient(polar_in, polar_out):
	Rif = rif_constants() * math.fabs(fresnel(polar_out, "vs") * fresnel(polar_out, "sb") * ((fresnel(polar_in, "vs") * fresnel(polar_in, "sb")) ** 2) * r_factors(polar_in, polar_out)) ** 2
	return Rif

# print nonlinear_reflection_coefficient(entry, exit)
# print chi("zzz")