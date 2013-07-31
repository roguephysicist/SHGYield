import math
import csv
from itertools import izip

## User Input
# Input
file_path = "/Users/sma/Downloads/res/" # full path to 'res' folder for desired structure
kpts = 514 # kpoints of the case you want to use
ecut = 25 # ecut of the case you want to use
# Polarization
entry = "p" # 's' or 'p'
exit = "p" # 's' or 'p'
# Angles
theta_deg = 65 # angle of incidence in degrees
phi_deg = 30 # azimuthal angle
# Output
out_file_path = file_path + "R" + entry + exit + "_" + str(kpts) + "_" + str(ecut)

## Housekeeping
n_0 = 6.02214129e+23 # electronic density - for testing purposes
theta = math.radians(theta_deg)
phi = math.radians(phi_deg)
omega = 1
epsilon = 1

## Functions
def nonlinear_reflection_coefficient(polar_in, polar_out):
	nrc = []
	ener = [row[0] for row in values("energy")]
	chi_zzz = [row[1] for row in values("zzz")]
	chi_zxx = [row[1] for row in values("zxx")]
	chi_xxz = [row[1] for row in values("xxz")]
	chi_xxx = [row[1] for row in values("xxx")]
	for energy, zzz, zxx, xxz, xxx in izip(ener, chi_zzz, chi_zxx, chi_xxz, chi_xxx):
		R = rif_constants() * math.fabs(fresnel(polar_out, "vs") * fresnel(polar_out, "sb") * ((fresnel(polar_in, "vs") * fresnel(polar_in, "sb")) ** 2) * r_factors(polar_in, polar_out, float(zzz), float(zxx), float(xxz), float(xxx))) ** 2
		nrc.append(R)
	return nrc

def values(component):
		if component == "energy":
			path = file_path + "shgC.sm_zzz_" + str(kpts) + "_half-slab_" + str(ecut) + "-nospin_scissor_0_Nc_29"
			energy = input_stage(path)
			return energy
		elif component == "linear":
			return
		else:
			path = file_path + "shgC.sm_" + component + "_" + str(kpts) + "_half-slab_" + str(ecut) + "-nospin_scissor_0_Nc_29"
			chi = input_stage(path)
			return chi

def input_stage(file):
	with open(file, 'rb') as csvfile:
		data = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
		comp_list=[[row[0], row[3]] for row in data]
		return comp_list

def rif_constants():
	const = (32 * 3.14159265359 ** 3 * omega ** 2) / (n_0 * 1.602176565e-19 ** 2 * 299792458.0 ** 3 * math.cos(theta) ** 2)
	return const

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

def r_factors(p_entry, p_exit, triperp, perpbipar, biparperp, tripar):
	if p_entry == "p" and p_exit == "p":
		r = math.sin(theta) * epsilon * (((math.sin(theta) ** 2) * (epsilon ** 2) * triperp) + (wave_vector("b") ** 2) * (epsilon ** 2) * perpbipar) + epsilon * epsilon * wave_vector("b") * wave_vector("b") * (-2 * math.sin(theta) * epsilon * biparperp + wave_vector("b") * epsilon * tripar * math.cos(3 * phi))
	elif p_entry == "s" and p_exit == "p":
		r = math.sin(theta) * epsilon * perpbipar - wave_vector("b") * epsilon * tripar * math.cos(3 * phi)
	elif p_entry == "p" and p_exit == "s":
		r = -(wave_vector("b") ** 2) * (epsilon ** 2) * tripar * math.sin(3 * phi)
	elif p_entry == "s" and p_exit == "s":
		r = tripar * math.sin(3 * phi)
	return r

def wave_vector(layer):
	if layer == "s":
		k = (omega / 299792458.0) * math.sqrt(epsilon - math.sin(theta) ** 2)
	elif layer == "b":
		k = (omega / 299792458.0) * math.sqrt(epsilon - math.sin(theta) ** 2)
	return k

def output_stage(energies, results, outfile):
 	with open(outfile, 'wb') as csvfile:
 		output = csv.writer(csvfile)
 		for row in results:
 			output.writerow(row)

print type(nonlinear_reflection_coefficient(entry, exit))
