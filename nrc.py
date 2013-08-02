import math
import csv
from itertools import izip

## User Input
# Input
file_path = "/Users/sma/Downloads/res/" # full path to 'res' folder for desired structure
kpts = 514 # kpoints of the case you want to use
ecut = 25 # ecut of the case you want to use
# Polarization
entry = "s" # 's' or 'p'
exit = "s" # 's' or 'p'
# Angles
theta_deg = 65 # angle of incidence in degrees
phi_deg = 30 # azimuthal angle

## Housekeeping
n_0 = 6.02214129e+23 # electronic density - for testing purposes
theta = math.radians(theta_deg)
phi = math.radians(phi_deg)

## Functions
def nonlinear_reflection_coefficient(polar_in, polar_out):
	nrc = []
	energy = [row[0] for row in values("epsilon")]
	epsilon = [row[1] for row in values("epsilon")]
	zzz = [row[0] for row in values("zzz")]
	zxx = [row[0] for row in values("zxx")]
	xxz = [row[0] for row in values("xxz")]
	xxx = [row[0] for row in values("xxx")]
	for energy, epsilon, zzz, zxx, xxz, xxx in izip(energy, epsilon, zzz, zxx, xxz, xxx):
		R = rif_constants(float(energy)) * math.fabs(fresnel(polar_out, "vs", float(energy), float(epsilon)) * fresnel(polar_out, "sb", float(energy), float(epsilon)) * ((fresnel(polar_in, "vs", float(energy), float(epsilon)) * fresnel(polar_in, "sb", float(energy), float(epsilon))) ** 2) * r_factors(polar_in, polar_out, float(zzz), float(zxx), float(xxz), float(xxx), float(energy), float(epsilon))) ** 2
		nrc.append(str(R))
	return nrc

def values(component):
	if component == "epsilon":
		path = file_path + "chi1.sm_xx_" + str(kpts) + "_" + str(ecut) + "-nospin_scissor_0_Nc_29"
		value = input_stage(path, "epsilon")
	else:
		path = file_path + "shgC.sm_" + component + "_" + str(kpts) + "_half-slab_" + str(ecut) + "-nospin_scissor_0_Nc_29"
		value = input_stage(path, "chi")
	return value

def input_stage(file, type):
	with open(file, 'rb') as csvfile:
		data = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
		if type == "epsilon":
			comp_list=[[row[0], row[2]] for row in data]
		elif type == "chi":
			comp_list=[[row[3]] for row in data]
	return comp_list

def rif_constants(energy):
	const = (32 * 3.14159265359 ** 3 * energy ** 2) / (n_0 * 1.602176565e-19 ** 2 * 299792458.0 ** 3 * math.cos(theta) ** 2)
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
	k = (energy / 299792458.0) * math.sqrt(math.fabs(epsilon - math.sin(theta) ** 2))
	return k

def output_stage():
	energy = [row[0] for row in values("epsilon")]
	#print omega
	print nonlinear_reflection_coefficient(entry, exit)
	#out_file = file_path + "R" + entry + exit + "_" + str(kpts) + "_" + str(ecut)
	#omega = [row[0] for row in values("epsilon")]
	#results = nonlinear_reflection_coefficient(entry, exit)
	#newdata = []
	#with open(out_file, 'wb') as csvfile:
 	#	output = csv.writer(csvfile)
	#	for c1, c2 in zip(omega, results):
	#		c = "%s  %s" % (c1, c2)
	#		newdata.append(c)
	#	print newdata
 	#	#for row in new_data:
 	#	#	output.writerow(row)

#print type(nonlinear_reflection_coefficient(entry, exit))
output_stage()

