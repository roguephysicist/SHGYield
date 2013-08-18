"""
nrc.py is a python program designed to calculate the Nonlinear reflection
coefficient for silicon surfaces. It works in conjunction with the matrix
elements calculated using ABINIT, and open source ab initio software,
and TINIBA, our in-house optical calculation software.

The work codified in this software can be found in Phys.Rev.B66, 195329(2002).
"""

from math import sin, cos, radians
from scipy import constants, interpolate
from numpy import loadtxt, savetxt, column_stack, absolute, sqrt, linspace

########### user input ###########
OUT_PATH = "./sample_data/nrc/"
CHI1_PATH = "./sample_data/res/chi1.kk_xx_yy_zz_3107_25-nospin_scissor_0.5_Nc_8"
ZZZ_PATH = "./sample_data/res/shgC.kk_zzz_145_half-slab_10-nospin_scissor_0.5_Nc_29"
ZXX_PATH = "./sample_data/res/shgC.kk_zxx_145_half-slab_10-nospin_scissor_0.5_Nc_29"
XXZ_PATH = "./sample_data/res/shgC.kk_xxz_145_half-slab_10-nospin_scissor_0.5_Nc_29"
XXX_PATH = "./sample_data/res/shgC.kk_xxx_145_half-slab_10-nospin_scissor_0.5_Nc_29"
# Angles
THETA_DEG = 65 # angle of incidence in degrees
PHI_DEG = 30 # azimuthal angle
# Misc
ELEC_DENS = 1.5e10 # electronic density
########### end user input ###########

## Housekeeping
THETA_RAD = radians(THETA_DEG)
PHI_RAD = radians(PHI_DEG)

########### functions ###########
def nonlinear_reflection():
    """ calls the different math functions and returns matrix,
    which is written to file """
    onee = linspace(0, 20, 2001)
    twoe = 2 * onee
    polarization = [["p", "p"], ["p", "s"], ["s", "p"], ["s", "s"]]
    for state in polarization:
        nrc = rif_constants(onee) * absolute(fresnel_vs(state[1], twoe) *
              fresnel_sb(state[1], twoe) * ((fresnel_vs(state[0], onee) *
              fresnel_sb(state[0], onee)) ** 2) *
              reflection_components(state[0], state[1], onee, twoe)) ** 2
        nrc = column_stack((onee, nrc))
        out = OUT_PATH + "R" + state[0] + state[1]
        save_matrix(out, nrc)

def chi1_real(energy):
    """ creates spline from real part of chi1 matrix"""
    energies = linspace(0, 20, 2001)
    chi1 = load_matrix(CHI1_PATH)
    interpolated = interpolate.InterpolatedUnivariateSpline(energies, chi1.real)
    return interpolated(energy)

def chi1_imag(energy):
    """ creates spline from imaginary part of chi1 matrix"""
    energies = linspace(0, 20, 2001)
    chi1 = load_matrix(CHI1_PATH)
    interpolated = interpolate.InterpolatedUnivariateSpline(energies, chi1.imag)
    return interpolated(energy)

def epsilon(energy):
    """ combines splines for real and imaginary parts of chi1 """
    chi1 = chi1_real(energy) + 1j * chi1_imag(energy)
    linear = 1 + (4 * constants.pi * chi1)
    return linear

def wave_vector(energy):
    """ math for wave vectors """
    k = (energy / constants.c) * sqrt(epsilon(energy) -
        (sin(THETA_RAD) ** 2))
    return k

def rif_constants(energy):
    """ math for constant term """
    const = (32 * (constants.pi ** 3) * (energy ** 2)) / ((ELEC_DENS *
        constants.e ** 2) * (constants.c ** 3) * (cos(THETA_RAD) ** 2))
    return const

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
    xxx = load_matrix(XXX_PATH)
    if polar_in == "p" and polar_out == "p":
        zzz = load_matrix(ZZZ_PATH)
        zxx = load_matrix(ZXX_PATH)
        xxz = load_matrix(XXZ_PATH)
        r_factor = sin(THETA_RAD) * epsilon(twoenergy) * \
                (((sin(THETA_RAD) ** 2) * (epsilon(energy) ** 2) * zzz) +
                (wave_vector(energy) ** 2) * (epsilon(energy) ** 2) * zxx) \
                 + epsilon(energy) * epsilon(twoenergy) * \
                 wave_vector(energy) * wave_vector(twoenergy) * \
                 (-2 * sin(THETA_RAD) * epsilon(energy) * xxz +
                wave_vector(energy) * epsilon(energy) * xxx *
                cos(3 * PHI_RAD))
    elif polar_in == "s" and polar_out == "p":
        zxx = load_matrix(ZXX_PATH)
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
    real, imaginary = loadtxt(in_file, unpack=True, usecols=[3, 4])
    data = real + 1j * imaginary
    return data

def save_matrix(out_file, data):
    """ saves matrix to file """
    savetxt(out_file, data, fmt=('%5.14e'), delimiter='\t')

nonlinear_reflection()
