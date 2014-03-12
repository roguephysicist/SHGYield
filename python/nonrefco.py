"""
nrc.py is a python program designed to calculate the Nonlinear reflection
coefficient for silicon surfaces. It works in conjunction with the matrix
elements calculated using ABINIT, and open source ab initio software,
and TINIBA, our in-house optical calculation software.

The work codified in this software can be found in Phys.Rev.B66, 195329(2002).
"""

import sys
from math import sin, cos, radians
from scipy import constants, interpolate
from numpy import loadtxt, savetxt, column_stack, absolute, \
                  sqrt, linspace, ones, complex128

########### user input ###########
LAYERS = sys.argv[1] # reads layers from command line
#KPOINTS = sys.argv[1] # reads kpoints from command line
#ECUT = sys.argv[1] # reads ecut from command line
COND = sys.argv[2] # reads N_c from command line
#LAYERS = 12
KPOINTS = 514
ECUT = 15
#COND = 75
SCIS = 0.83808
OUT = "../calculated/nrc/"
# Angles
THETA_RAD = radians(65)
PHI_RAD = radians(30)
# Misc
ELEC_DENS = 1e-28 # electronic density and scaling factor (1e-7 * 1e-21)
ENERGIES = linspace(0.01, 20, 2000)

#if LAYERS == "12":
#    SCIS = 1.8679
#elif LAYERS == "18":
#    SCIS = 1.91908
#elif LAYERS == "24":
#    SCIS = 1.94138
#elif LAYERS == "30":
#    SCIS = 1.95309
#elif LAYERS == "36":
#    SCIS = 1.96003

""" these are the paths to the appropriate response files """
CHI1 = "../calculated/bulk_chi1/chi1.sm_xx_yy_zz_3107_25-nospin_scissor_0.83808_Nc_26"
#CHI1 = "../calculated/" + str(LAYERS) + "/chi1.sm_xx_yy_zz_" + str(KPOINTS)\
#        + "_" + str(ECUT) + "-nospin_scissor_" + str(SCIS) + "_Nc_" + str(COND)
ZZZ = "../calculated/shg/" + str(LAYERS) + "/shgC.sm_zzz_" + str(KPOINTS) + \
"_half-slab_" + str(ECUT) + "-nospin_scissor_" + str(SCIS) + "_Nc_" + str(COND)
ZXX = "../calculated/shg/" + str(LAYERS) + "/shgC.sm_zxx_" + str(KPOINTS) + \
"_half-slab_" + str(ECUT) + "-nospin_scissor_" + str(SCIS) + "_Nc_" + str(COND)
XXZ = "../calculated/shg/" + str(LAYERS) + "/shgC.sm_xxz_" + str(KPOINTS) + \
"_half-slab_" + str(ECUT) + "-nospin_scissor_" + str(SCIS) + "_Nc_" + str(COND)
XXX = "../calculated/shg/" + str(LAYERS) + "/shgC.sm_xxx_" + str(KPOINTS) + \
"_half-slab_" + str(ECUT) + "-nospin_scissor_" + str(SCIS) + "_Nc_" + str(COND)

########### functions ###########
def nonlinear_reflection():
    """ calls the different math functions and returns matrix,
    which is written to file """
    onee = linspace(0.01, 20, 2000)
    twoe = 2 * onee
    polarization = [["p", "p"], ["p", "s"], ["s", "p"], ["s", "s"]]
    for state in polarization:
        nrc = rif_constants(onee) * absolute(((fresnel_vs(state[1], twoe) *
              fresnel_sb(state[1], twoe) * ((fresnel_vs(state[0], onee) *
              fresnel_sb(state[0], onee)) ** 2)) / 2j) *
              reflection_components(state[0], state[1], onee, twoe)) ** 2
        nrc = column_stack((onee, nrc))
        out = OUT + str(LAYERS) + "_" + "R" + state[0] + state[1] + "_k" +\
              str(KPOINTS).zfill(4) + "_e" + str(ECUT) + "_" + str(SCIS) +\
              "_" + str(COND)
        save_matrix(out, nrc)

def chi_one(part, energy):
    """ creates spline from real part of chi1 matrix"""
    chi1 = load_matrix_linear(CHI1)
    interpolated = \
    interpolate.InterpolatedUnivariateSpline(ENERGIES, getattr(chi1, part))
    return interpolated(energy)

def epsilon(energy):
    """ combines splines for real and imaginary parts of chi1 """
    chi1 = chi_one("real", energy) + 1j * chi_one("imag", energy)
    linear = 4 * constants.pi * chi1
    return linear

def wave_vector(energy):
    """ math for wave vectors """
    k = sqrt(epsilon(energy) - (sin(THETA_RAD) ** 2))
    return k

def rif_constants(energy):
    """ math for constant term """
    const = (32 * (constants.pi ** 3) * (energy ** 2)) / (ELEC_DENS *
              ((constants.c * 100) ** 3) * (cos(THETA_RAD) ** 2) *
              (constants.value("Planck constant over 2 pi in eV s") ** 2))
    return const

def electrostatic_units(energy):
    """ coefficient to convert to appropriate electrostatic units """
    area = (1 / ((2 * sqrt(2)) ** 2)) * 2 * sqrt(3)
    factor = (1j * ((2 *
              constants.value("Rydberg constant times hc in eV")) ** 5) *
              1e-5 * 2.08e-15 *
            ((constants.value("lattice parameter of silicon") / 1e-10) ** 3))\
              / (area * ((energy) ** 3))
    return factor

def fresnel_vs(polarization, energy):
    """ math for fresnel factors from vacuum to surface """
    if polarization == "s":
        fresnel = (2 * cos(THETA_RAD)) / (cos(THETA_RAD) +
                   wave_vector(energy))
    elif polarization == "p":
        fresnel = (2 * cos(THETA_RAD)) / (epsilon(energy) *
                   cos(THETA_RAD) + wave_vector(energy))
    return fresnel

def fresnel_sb(polarization, energy):
    """ math for fresnel factors from surface to bulk. Fresnel model """
    if polarization == "s":
        fresnel = ones(2000, dtype=complex128)
        #fresnel = (2 * wave_vector(energy)) / (wave_vector(energy)
        #           + wave_vector(energy))
    elif polarization == "p":
        fresnel = 1 / epsilon(energy)
        #fresnel = (2 * wave_vector(energy)) / (epsilon(energy) *
        #wave_vector(energy) + epsilon(energy) * wave_vector(energy))
    return fresnel

def reflection_components(polar_in, polar_out, energy, twoenergy):
    """ math for different r factors. loads in different component matrices.
        the bms program had electrostatic units multiplying each - caution"""
    zzz = load_matrix(ZZZ)
    zxx = load_matrix(ZXX)
    xxz = load_matrix(XXZ)
    xxx = load_matrix(XXX)
    if polar_in == "p" and polar_out == "p":
        r_factor = sin(THETA_RAD) * epsilon(twoenergy) * \
                (((sin(THETA_RAD) ** 2) * (epsilon(energy) ** 2) * zzz) +
                (wave_vector(energy) ** 2) * (epsilon(energy) ** 2) * zxx) \
                 + epsilon(energy) * epsilon(twoenergy) * \
                 wave_vector(energy) * wave_vector(twoenergy) * \
                 (-2 * sin(THETA_RAD) * epsilon(energy) * xxz +
                wave_vector(energy) * epsilon(energy) * xxx *
                cos(3 * PHI_RAD))
    elif polar_in == "s" and polar_out == "p":
        r_factor = sin(THETA_RAD) * epsilon(twoenergy) * zxx - \
               wave_vector(twoenergy) * epsilon(twoenergy) * \
               xxx * cos(3 * PHI_RAD)
    elif polar_in == "p" and polar_out == "s":
        r_factor = -(wave_vector(energy) ** 2) * (epsilon(energy) ** 2) * \
                     xxx * sin(3 * PHI_RAD)
    elif polar_in == "s" and polar_out == "s":
        r_factor = xxx * sin(3 * PHI_RAD)
    return r_factor

def load_matrix(in_file):
    """ loads files into matrices and extracts columns """
    real1w, imaginary1w, real2w, imaginary2w = loadtxt(in_file, unpack=True, usecols=[1, 2, 3, 4], skiprows=1)
    real = real1w + real2w
    imaginary = imaginary1w + imaginary2w
    data = real + 1j * imaginary
    return data

def load_matrix_linear(in_file):
    """ loads files into matrices and extracts columns """
    real, imaginary = loadtxt(in_file, unpack=True, usecols=[1, 2], skiprows=1)
    data = real + 1j * imaginary
    return data

def save_matrix(out_file, data):
    """ saves matrix to file """
    savetxt(out_file, data, fmt=('%5.14e'), delimiter='\t')

nonlinear_reflection()
