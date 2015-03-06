#!/Users/sma/anaconda/bin/python
"""
nrc.py is a python program designed to calculate the Nonlinear reflection
coefficient for silicon surfaces. It works in conjunction with the matrix
elements calculated using ABINIT, and open source ab initio software,
and TINIBA, our in-house optical calculation software.

The work codified in this software can be found in Phys.Rev.B66, 195329(2002).
Please refer to "Strain induced SHG" manuscript (BMS) for equation references
unless stated otherwise.

Experimental gap for bulk Si = 3.4 according to Landolt-Boernstein

TO-DO
- read files for batch processing from an input file.
"""

import sys
from math import sin, cos, radians
from scipy import constants, interpolate
from numpy import loadtxt, savetxt, column_stack, absolute, \
                  sqrt, linspace, ones, complex128

########### User Input ###########

#CHI = sys.argv[1] # reads layers from command line
#XXX = sys.argv[2] # reads kpoints from command line
#XXZ = sys.argv[3] # reads ecut from command line
#ZXX = sys.argv[4] # reads N_c from command line
#ZZZ = sys.argv[5] # reads scissor correction from command line
##OUT = "../calculated/nrc/"

# Angles
THETA_RAD = radians(65)
PHI_RAD = radians(30)
# Misc
ENERGY = linspace(0.01, 20, 2000) # for debugging only, don't need this motherfucking shit

########### Functions ###########

### Equations 
###
def nonlinear_reflection(state, energy):
    """
    eq. 22: R_{if}
    dependencies: polarization state arrays
    parents: control
    children: rif_constants, fresnel_vl, fresnel_lb, reflection_components
    Calls math functions and returns numpy array for each polarization
    """
    onee = energy
    twoe = 2 * onee
    nrc = rif_constants(onee) * absolute(((fresnel_vl(state[1], twoe) *
          fresnel_lb(state[1], twoe) * ((fresnel_vl(state[0], onee) *
          fresnel_lb(state[0], onee)) ** 2)) / 2j) *
          reflection_components(state[0], state[1], onee, twoe)) ** 2
    return nrc

def rif_constants(energy):
    """
    eq. 22: constant term outside absolute value
    dependencies: 1w energy array
    parents: nonlinear_reflection
    children: none
    Multiplies constants. "electrondensity" merits revision.
    """
    electrondensity = 1e-28 # electronic density and scaling factor (1e-7 * 1e-21)
    const = (32 * (constants.pi ** 3) * ((energy / constants.value("Planck constant over 2 pi in eV s")) ** 2)) / (electrondensity * ((constants.c * 100) ** 3) * (cos(THETA_RAD) ** 2))
    return const

def fresnel_vl(polarization, energy):
    """
    eq. 26: t^{vl}_{s} and t^{vl}_{p}
    dependencies: energy array
    parents: nonlinear_reflection
    children: epsilon, wave_vector
    Calculates fresnel factors for vacuum to surface
    """
    if polarization == "s":
        fresnel = (2 * cos(THETA_RAD)) / (cos(THETA_RAD) +
                   wave_vector(energy))
    elif polarization == "p":
        fresnel = (2 * cos(THETA_RAD)) / (epsilon(energy) *
                   cos(THETA_RAD) + wave_vector(energy))
    return fresnel

def fresnel_lb(polarization, energy):
    """
    eq. 27: t^{lb}_{s} and t^{lb}_{p}
    dependencies: energy array
    parents: nonlinear_reflection
    children: epsilon, wave_vector
    Calculates fresnel factors for surface to bulk
    """
    if polarization == "s":
        fresnel = ones(2000, dtype=complex128)
        # fresnel = (2 * wave_vector(energy)) / (wave_vector(energy)
        #              + wave_vector(energy))
    elif polarization == "p":
        fresnel = 1 / epsilon(energy)
        # fresnel = (2 * wave_vector(energy)) / (epsilon(energy)
        #              * wave_vector(energy) + epsilon(energy)
        #              * wave_vector(energy))
    return fresnel

def reflection_components(polar_in, polar_out, energy, twoenergy):
    """
    eqs. 34a--34d
    dependencies: polarization states, 1w and 2w energy arrays
    parents: nonlinear_reflection
    children: load_shg, epsilon, wave_vector
    Calculates r_{if} factors. Loads shg arrays and does a lot of matrix
    operations. BMS program had electrostatic units multiplying shg arrays.
    """
    # zzz = electrostatic_units(energy) * load_shg(ZZZ)
    # zxx = electrostatic_units(energy) * load_shg(ZXX)
    # xxz = electrostatic_units(energy) * load_shg(XXZ)
    # xxx = electrostatic_units(energy) * load_shg(XXX)
    zzz = load_shg(ZZZ)
    zxx = load_shg(ZXX)
    xxz = load_shg(XXZ)
    xxx = load_shg(XXX)
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

def epsilon(energy):
    """
    eq. ?? REFERENCE NEEDED
    dependencies: energy array
    parents: wave_vector, fresnel_vl, fresnel_lb, reflection_components
    children: chi_spline
    Combines splines for real and imaginary parts of Chi^(1)
    """
    chi1 = chi_spline("real", energy) + 1j * chi_spline("imag", energy)
    linear = 4 * constants.pi * chi1
    return linear

def wave_vector(energy):
    """
    eq. 5: k_{z}(\omega)
    dependencies: energy array
    parents: fresnel_vl, fresnel_lb, reflection_components
    children: epsilon
    Calculates wave vector k
    """
    k = sqrt(epsilon(energy) - (sin(THETA_RAD) ** 2))
    return k

def electrostatic_units(energy):
    """
    coefficient to convert to appropriate electrostatic units
    possibly deprecated?
    """
    area = (1 / ((2 * sqrt(2)) ** 2)) * 2 * sqrt(3)
    factor = (1j * ((2 *
              constants.value("Rydberg constant times hc in eV")) ** 5) *
              1e-5 * 2.08e-15 *
            ((constants.value("lattice parameter of silicon") / 1e-10) ** 3))\
              / (area * ((energy) ** 3))
    return factor

### Control Functions
###
def control():
    """
    CHAIN START
    parents: none
    children: nonlinear_reflection, save_matrix
    Creates final matrix and writes to file.
    """
    onee = linspace(0.01, 20, 2000)
    nrc = column_stack((onee, nonlinear_reflection(["p", "p"], onee),
                              nonlinear_reflection(["p", "s"], onee),
                              nonlinear_reflection(["s", "p"], onee),
                              nonlinear_reflection(["s", "s"], onee)))
    out  = "reflex.dat"
    save_matrix(out, nrc)

def load_chi(in_file):
    """
    dependencies: input file
    parents: chi_spline
    children: none
    Loads Chi^(1) file, unpacks columns, and combines into complex numpy array.
    """
    real, imaginary = loadtxt(in_file, unpack=True, usecols=[1, 2], skiprows=1)
    data = real + 1j * imaginary
    return data

def load_shg(in_file):
    """
    dependencies: input file
    parents: reflection_components
    children: none
    Loads shg Chi^(2) files, unpacks columns, sums 1w and 2w for real and 
    imag, and combines into complex numpy array.
    """
    real1w, imaginary1w, real2w, imaginary2w = loadtxt(in_file, unpack=True, usecols=[1, 2, 3, 4], skiprows=1)
    real = real1w + real2w
    imaginary = imaginary1w + imaginary2w
    data = real + 1j * imaginary
    return data

def save_matrix(out_file, data):
    """
    CHAIN END
    dependencies: output file, data to be written
    parents: control
    children: none
    Saves final numpy array to output file, and writes fancy header to file. 
    FMT value for energy column differs from reflection components. 
    """
    savetxt(out_file, data, fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e'),
                            delimiter='    ', 
                            header='w(eV)  Rpp                     Rps                     Rsp                     Rss')

def chi_spline(part, energy):
    """
    dependencies: real or imaginary indicator, energy
    parents: epsilon
    children: load_chi
    Creates spline from real part of Chi^(1) array
    """
    chi1 = load_chi(CHI)
    interpolated = \
    interpolate.InterpolatedUnivariateSpline(energy, getattr(chi1, part))
    return interpolated(energy)

control()
