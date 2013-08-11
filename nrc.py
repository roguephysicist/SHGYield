import math
import random
from numpy import *
from scipy import constants

########### user input ###########
#### XYZ file disorder
disorder_in = "./sample_data/disorder/si_h_14.xyz" # full path to xyz file you want to disorder
l = 2.351 # bond length in angstrom for Si
atoms = [2, 3, 4] # atoms you want to disorder
disorder_amount = [0.1, 0.1, 0.1] # amount you want to disorder each atom
repeat = 1 # number of times you want to repeat the disordering

#### Nonlinear reflection coefficient
nrc_path = "./sample_data/res/" # full path to 'res' folder for desired structure
kpts = 145 # kpoints of the case you want to use
ecut = 10 # ecut of the case you want to use
# Polarization
entry = "s" # 's' or 'p'
exit = "s" # 's' or 'p'
# Angles
theta_deg = 65 # angle of incidence in degrees
phi_deg = 30 # azimuthal angle
# Misc
n_0 = 50 # electronic density
########### end user input ###########

## Housekeeping
theta = math.radians(theta_deg)
phi = math.radians(phi_deg)

########### functions ###########
def disorder():
	count = 1
	while count <= repeat:
		selected = [x - 1 for x in atoms]
		bl_bohr = l * constants.angstrom / constants.value("Bohr radius")
		xyz = load_matrix("disorder")
		data = xyz[selected]
		new_xyz = xyz.copy()
		final = zip(data, disorder_amount)
		for atom in range(0, len(final)):
			polar = random.random() * constants.pi
			azimuthal = random.random() * 2 * constants.pi
			new_xyz[selected[atom], 0] = final[atom][0][0] + (final[atom][1] * (bl_bohr / 2) * math.sin(polar) * math.cos(azimuthal))
			new_xyz[selected[atom], 1] = final[atom][0][1] + (final[atom][1] * (bl_bohr / 2) * math.sin(polar) * math.sin(azimuthal))
			new_xyz[selected[atom], 2] = final[atom][0][2] + (final[atom][1] * (bl_bohr / 2) * math.cos(polar))
		out = "./sample_data/disorder/si_h_14_mod_" + str(count).zfill(3) + ".xyz"
		save_matrix(out, new_xyz)
		count += 1

def nonlinear_reflection_coefficient():
	e, chi1 = load_matrix("chi1")
	twoe, chi1twoe = load_matrix("twoe")
	zzz = load_matrix("zzz")
	zxx = load_matrix("zxx")
	xxz = load_matrix("xxz")
	xxx = load_matrix("xxx")
	R = rif_constants(e) * fabs(fresnel(exit, "vs", twoe, chi1twoe) * fresnel(exit, "sb", twoe, chi1twoe) * ((fresnel(entry, "vs", e, chi1) * fresnel(entry, "sb", e, chi1)) ** 2) * r_factors(entry, exit, zzz, zxx, xxz, xxx, e, twoe, chi1, chi1twoe)) ** 2
	R = column_stack((e, R))
	out = nrc_path + "R" + entry + exit + "_" + str(kpts) + "_" + str(ecut)
	save_matrix(out, R)

def rif_constants(energy):
	const = (32 * (constants.pi ** 3) * (energy ** 2)) / ((n_0 * constants.e ** 2) * (constants.c ** 3) * (math.cos(theta) ** 2))
	return const

def fresnel(polarization, material, energy, chi1):
	if polarization == "s" and material == "vs":
		t = (2 * math.cos(theta)) / (math.cos(theta) + wave_vector(energy, chi1))
	elif polarization == "s" and material == "sb":
		t = (2 * wave_vector(energy, chi1)) / (wave_vector(energy, chi1) + wave_vector(energy, chi1))
	elif polarization == "p" and material == "vs":
		t = (2 * math.cos(theta)) / (epsilon(chi1) * math.cos(theta) + wave_vector(energy, chi1))
	elif polarization == "p" and material == "sb":
		t = (2 * wave_vector(energy, chi1)) / (epsilon(chi1) * wave_vector(energy, chi1) + epsilon(chi1) * wave_vector(energy, chi1))
	return t

def r_factors(polar_in, polar_out, triperp, perpbipar, biparperp, tripar, energy, twoe, chi1, chi1twoe):
	if polar_in == "p" and polar_out == "p":
		r = math.sin(theta) * epsilon(chi1twoe) * (((math.sin(theta) ** 2) * (epsilon(chi1) ** 2) * triperp) + (wave_vector(energy, chi1) ** 2) * (epsilon(chi1) ** 2) * perpbipar) + epsilon(chi1) * epsilon(chi1twoe) * wave_vector(energy, chi1) * wave_vector(twoe, chi1twoe) * (-2 * math.sin(theta) * epsilon(chi1) * biparperp + wave_vector(energy, chi1) * epsilon(chi1) * tripar * math.cos(3 * phi))
	elif polar_in == "s" and polar_out == "p":
		r = math.sin(theta) * epsilon(chi1twoe) * perpbipar - wave_vector(twoe, chi1twoe) * epsilon(chi1twoe) * tripar * math.cos(3 * phi)
	elif polar_in == "p" and polar_out == "s":
		r = -(wave_vector(energy, chi1) ** 2) * (epsilon(chi1) ** 2) * tripar * math.sin(3 * phi)
	elif polar_in == "s" and polar_out == "s":
		r = tripar * math.sin(3 * phi)
	return r

def epsilon(chi1):
	epsilon = 1 + (4 * constants.pi * chi1)
	return epsilon

def wave_vector(energy, chi1):
	k = (energy / constants.c) * sqrt(epsilon(chi1) - (math.sin(theta) ** 2))
	return k

def load_matrix(comp):
	if comp == "disorder":
		data = loadtxt(disorder_in)
	elif comp == "chi1":
		f = nrc_path + "chi1.kk_xx_yy_zz_" + str(kpts) + "_" + str(ecut) + "-nospin_scissor_0_Nc_29"
		data = genfromtxt(f, skip_footer=1000, unpack=True, usecols=[0, 6])
	elif comp == "twoe":
		f = nrc_path + "chi1.kk_xx_yy_zz_" + str(kpts) + "_" + str(ecut) + "-nospin_scissor_0_Nc_29"
		new = column_stack(loadtxt(f, unpack=True, usecols=[0, 6]))
		wanted = array(2 * new[:,0]).reshape(-1).tolist()
		data = column_stack(new[logical_or.reduce([new[:,0] == x for x in wanted])])
	elif comp == "zzz" or comp == "zxx" or comp == "xxz" or comp == "xxx":
		f = nrc_path + "shgC.kk_" + comp + "_" + str(kpts) + "_half-slab_" + str(ecut) + "-nospin_scissor_0_Nc_29"
		data = genfromtxt(f, skip_footer=1000, unpack=True, usecols=[4])
	return data

def save_matrix(f, data):
	savetxt(f, data, fmt=('%5.14e'), delimiter='\t')

#disorder()
nonlinear_reflection_coefficient()