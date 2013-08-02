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
# Output
out_file_path = file_path + "R" + entry + exit + "_" + str(kpts) + "_" + str(ecut)

## Housekeeping
n_0 = 6.02214129e+23 # electronic density - for testing purposes
theta = math.radians(theta_deg)
phi = math.radians(phi_deg)

## Functions
def nonlinear_reflection_coefficient(polar_in, polar_out):
	nrc = []
	omeg = [row[0] for row in values("linear")]
	epsi = [row[1] for row in values("linear")]
	chi_zzz = [row[1] for row in values("zzz")]
	chi_zxx = [row[1] for row in values("zxx")]
	chi_xxz = [row[1] for row in values("xxz")]
	chi_xxx = [row[1] for row in values("xxx")]
	for omega, epsilon, zzz, zxx, xxz, xxx in izip(omeg, epsi, chi_zzz, chi_zxx, chi_xxz, chi_xxx):
		R = rif_constants(float(omega)) * math.fabs(fresnel(polar_out, "vs", float(omega), float(epsilon)) * fresnel(polar_out, "sb", float(omega), float(epsilon)) * ((fresnel(polar_in, "vs", float(omega), float(epsilon)) * fresnel(polar_in, "sb", float(omega), float(epsilon))) ** 2) * r_factors(polar_in, polar_out, float(zzz), float(zxx), float(xxz), float(xxx), float(omega), float(epsilon))) ** 2
		nrc.append(R)
	omeg.extend(nrc)
	return omeg

def values(component):
		if component == "linear":
			path = file_path + "chi1.sm_xx_" + str(kpts) + "_" + str(ecut) + "-nospin_scissor_0_Nc_29"
			epsilon = input_stage(path, "epsilon")
			return epsilon
		else:
			path = file_path + "shgC.sm_" + component + "_" + str(kpts) + "_half-slab_" + str(ecut) + "-nospin_scissor_0_Nc_29"
			chi = input_stage(path, "chi")
			return chi

def input_stage(file, type):
	if type == "epsilon":
		with open(file, 'rb') as csvfile:
			data = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
			comp_list=[[row[0], row[2]] for row in data]
			return comp_list
	elif type == "chi":
		with open(file, 'rb') as csvfile:
			data = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
			comp_list=[[row[0], row[3]] for row in data]
			return comp_list

def rif_constants(frequency):
	const = (32 * 3.14159265359 ** 3 * frequency ** 2) / (n_0 * 1.602176565e-19 ** 2 * 299792458.0 ** 3 * math.cos(theta) ** 2)
	return const

def fresnel(polarization, material, frequency, linear): # REMOVE PLUS ONE FROM LINE 69 ASAP
	if polarization == "s" and material == "vs":
		t = (2 * math.cos(theta)) / (math.cos(theta) + wave_vector("s", frequency, linear))
	elif polarization == "s" and material == "sb":
		t = (2 * wave_vector("s", frequency, linear)) / (1 + wave_vector("s", frequency, linear) + wave_vector("b", frequency, linear))
	elif polarization == "p" and material == "vs":
		t = (2 * math.cos(theta)) / (linear * math.cos(theta) + wave_vector("s", frequency, linear))
	elif polarization == "p" and material == "sb":
		t = (2 * wave_vector("s", frequency, linear)) / (linear * wave_vector("s", frequency, linear) + linear * wave_vector("b", frequency, linear))
	return t

def r_factors(p_entry, p_exit, triperp, perpbipar, biparperp, tripar, frequency, linear):
	if p_entry == "p" and p_exit == "p":
		r = math.sin(theta) * linear * (((math.sin(theta) ** 2) * (linear ** 2) * triperp) + (wave_vector("b", frequency, linear) ** 2) * (linear ** 2) * perpbipar) + linear * linear * wave_vector("b", frequency, linear) * wave_vector("b", frequency, linear) * (-2 * math.sin(theta) * linear * biparperp + wave_vector("b", frequency, linear) * linear * tripar * math.cos(3 * phi))
	elif p_entry == "s" and p_exit == "p":
		r = math.sin(theta) * linear * perpbipar - wave_vector("b", frequency, linear) * linear * tripar * math.cos(3 * phi)
	elif p_entry == "p" and p_exit == "s":
		r = -(wave_vector("b", frequency, linear) ** 2) * (linear ** 2) * tripar * math.sin(3 * phi)
	elif p_entry == "s" and p_exit == "s":
		r = tripar * math.sin(3 * phi)
	return r

def wave_vector(layer, frequency, linear): # the absolute value thing is just temporary
	if layer == "s":
		k = (frequency / 299792458.0) * math.sqrt(math.fabs(linear - math.sin(theta) ** 2))
	elif layer == "b":
		k = (frequency / 299792458.0) * math.sqrt(math.fabs(linear - math.sin(theta) ** 2))
	return k

def output_stage(energies, outfile):
	results = nonlinear_reflection_coefficient(entry, exit)
 	with open(outfile, 'wb') as csvfile:
 		output = csv.writer(csvfile)
 		for row in results:
 			output.writerow(row)

print nonlinear_reflection_coefficient(entry, exit)
