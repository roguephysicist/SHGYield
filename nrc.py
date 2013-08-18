"""
nrc.py is a python program designed to calculate the Nonlinear reflection
coefficient for silicon surfaces. It works in conjunction with the matrix
elements calculated using ABINIT, and open source ab initio software,
and TINIBA, our in-house optical calculation software.

The work codified in this software can be found in Phys.Rev.B66,195329(2002).
"""

import math
import random
from scipy import constants, interpolate
from numpy import loadtxt, savetxt, column_stack, absolute, sqrt, linspace
#import matplotlib.pyplot as plt

########### user input ###########
#### XYZ file disorder
DISORDER_PATH = "./sample_data/disorder/si_h_14.xyz"
BOND_LENGTH = 2.351
ATOMS_DISORDERED = [2, 3, 4]
DISORDER_AMOUNT = [0.1, 0.1, 0.1]
REPEAT = 1

#### Nonlinear reflection coefficient
OUT_PATH = "./sample_data/Rpp_145_10"
CHI1_PATH = "./sample_data/chi1.kk_xx_yy_zz_3107_25-nospin_scissor_0.5_Nc_8"
ZZZ_PATH = "./sample_data/shgC.kk_zzz_145_half-slab_10-nospin_scissor_0.5_Nc_29"
ZXX_PATH = "./sample_data/shgC.kk_zxx_145_half-slab_10-nospin_scissor_0.5_Nc_29"
XXZ_PATH = "./sample_data/shgC.kk_xxz_145_half-slab_10-nospin_scissor_0.5_Nc_29"
XXX_PATH = "./sample_data/shgC.kk_xxx_145_half-slab_10-nospin_scissor_0.5_Nc_29"
# Polarization
ENTRY = "s" # 's' or 'p'
EXIT = "s" # 's' or 'p'
# Angles
THETA_DEG = 65 # angle of incidence in degrees
PHI_DEG = 30 # azimuthal angle
# Misc
ELEC_DENS = 1.5e10 # electronic density
########### end user input ###########

## Housekeeping
THETA_RAD = math.radians(THETA_DEG)
PHI_RAD = math.radians(PHI_DEG)

########### functions ###########
def disorder():
    """ extracts selected rows from input, converts to matrix,
    disorders, and returns matrix which is written to file """
    count = 1
    while count <= REPEAT:
        selected = [x - 1 for x in ATOMS_DISORDERED]
        bl_bohr = BOND_LENGTH * constants.angstrom /\
                     constants.value("Bohr radius")
        xyz = loadtxt(DISORDER_PATH)
        data = xyz[selected]
        new_xyz = xyz.copy()
        final = zip(data, DISORDER_AMOUNT)
        for atom in range(0, len(final)):
            polar = random.random() * constants.pi
            azimuthal = random.random() * 2 * constants.pi
            new_xyz[selected[atom], 0] = final[atom][0][0] + (final[atom][1] *
                    (bl_bohr / 2) * math.sin(polar) * math.cos(azimuthal))
            new_xyz[selected[atom], 1] = final[atom][0][1] + (final[atom][1] *
                    (bl_bohr / 2) * math.sin(polar) * math.sin(azimuthal))
            new_xyz[selected[atom], 2] = final[atom][0][2] + (final[atom][1] *
                    (bl_bohr / 2) * math.cos(polar))
        out = "si_h_14_mod_" + str(count).zfill(3) + ".xyz"
        save_matrix(out, new_xyz)
        count += 1

def nonlinear_reflection():
    """ calls the different math functions and returns matrix,
    which is written to file """
    oneenergy = linspace(0, 20, 2001)
    twoenergy = 2 * oneenergy
    nrc = rif_constants(oneenergy) * absolute(fresnel_vs(EXIT, twoenergy) *
          fresnel_sb(EXIT, twoenergy) * ((fresnel_vs(ENTRY, oneenergy) *
          fresnel_sb(ENTRY, oneenergy)) ** 2) *
          reflection_components(ENTRY, EXIT, oneenergy, twoenergy)) ** 2
    print len(nrc)
    print nrc
    nrc = column_stack((oneenergy, nrc))
    save_matrix(OUT_PATH, nrc)
    #plt.plot(oneenergy,nrc)
    #plt.show()

def epsilon(energy):
    """ math to convert from chi1 to epsilon """
    energies = linspace(0, 20, 2001)
    chi1 = load_complex_matrix(CHI1_PATH)
    linear = 1 + (4 * constants.pi * chi1)
    interpolated = interpolate.InterpolatedUnivariateSpline(energies, linear)
    return interpolated(energy)

def wave_vector(energy):
    """ math for wave vectors """
    k = (energy / constants.c) * sqrt(epsilon(energy) -
        (math.sin(THETA_RAD) ** 2))
    return k

def rif_constants(energy):
    """ math for constant term """
    const = (32 * (constants.pi ** 3) * (energy ** 2)) / ((ELEC_DENS *
        constants.e ** 2) * (constants.c ** 3) * (math.cos(THETA_RAD) ** 2))
    return const

def fresnel_vs(polarization, energy):
    """ math for fresnel factors from vacuum to surface """
    if polarization == "s":
        fresnel = (2 * math.cos(THETA_RAD)) / (math.cos(THETA_RAD) +
                   wave_vector(energy))
    elif polarization == "p":
        fresnel = (2 * math.cos(THETA_RAD)) / (epsilon(energy) *
                   math.cos(THETA_RAD) + wave_vector(energy))
    return fresnel

def fresnel_sb(polarization, energy):
    """ math for fresnel factors from surface to bulk """
    if polarization == "s":
        fresnel = (2 * wave_vector(energy)) / (wave_vector(energy)
                   + wave_vector(energy))
    elif polarization == "p":
        fresnel = (2 * wave_vector(energy)) / (epsilon(energy) *
        wave_vector(energy) + epsilon(energy) * wave_vector(energy))
    return fresnel

def reflection_components(polar_in, polar_out, energy, twoenergy):
    """ math for different r factors. loads in different component matrices """
    xxx = load_complex_matrix(XXX_PATH)
    if polar_in == "p" and polar_out == "p":
        zzz = load_complex_matrix(ZZZ_PATH)
        zxx = load_complex_matrix(ZXX_PATH)
        xxz = load_complex_matrix(XXZ_PATH)
        r_factor = math.sin(THETA_RAD) * epsilon(twoenergy) * \
                (((math.sin(THETA_RAD) ** 2) * (epsilon(energy) ** 2) * zzz) +
                (wave_vector(energy) ** 2) * (epsilon(energy) ** 2) * zxx) \
                 + epsilon(energy) * epsilon(twoenergy) * \
                 wave_vector(energy) * wave_vector(twoenergy) * \
                 (-2 * math.sin(THETA_RAD) * epsilon(energy) * xxz +
                wave_vector(energy) * epsilon(energy) * xxx *
                math.cos(3 * PHI_RAD))
    elif polar_in == "s" and polar_out == "p":
        zxx = load_complex_matrix(ZXX_PATH)
        r_factor = math.sin(THETA_RAD) * epsilon(twoenergy) * zxx - \
               wave_vector(twoenergy) * epsilon(twoenergy) * \
               xxx * math.cos(3 * PHI_RAD)
    elif polar_in == "p" and polar_out == "s":
        r_factor = -(wave_vector(energy) ** 2) * (epsilon(energy) ** 2) * \
                     xxx * math.sin(3 * PHI_RAD)
    elif polar_in == "s" and polar_out == "s":
        r_factor = xxx * math.sin(3 * PHI_RAD)
    return r_factor

def load_complex_matrix(in_file):
    """ loads files into matrices and extracts columns """
    real, imaginary = loadtxt(in_file, unpack=True, usecols=[3, 4])
    data = real + 1j * imaginary
    return data

def save_matrix(out_file, data):
    """ saves matrix to file """
    savetxt(out_file, data, fmt=('%5.14e'), delimiter='\t')

#disorder()
#nonlinear_reflection()
ENERGY = linspace(0, 20, 2001)
print len(epsilon(ENERGY))
print type(epsilon(ENERGY))
