#!/Users/sma/anaconda/bin/python
"""
nrc.py is a python program designed to calculate the Nonlinear reflection
coefficient for silicon surfaces. It works in conjunction with the matrix
elements calculated umath.sing ABINIT, and open source ab initio software,
and TINIBA, our in-house optical calculation software.

The work codified in this software can be found in Phys.Rev.B66, 195329(2002).
Please refer to "Strain induced SHG" manuscript (BMS) for equation references
unless stated otherwise.

Experimental gap for bulk Si = 3.4 according to Landolt-Boernstein
"""

import sys
import math
import numpy as np
from scipy import constants, interpolate

# Angles and energies
THETA_RAD = math.radians(65)
PHI_RAD = math.radians(30)
ONEE = np.linspace(0.01, 10, 1000)

########### Functions ###########

### Equations
###
def nonlinear_reflection(state):
    """
    eq. 22: R_{if}
    dependencies: polarization state arrays
    parents: control
    children: rif_constants, fresnel_vl, fresnel_lb, reflection_components
    Calls math functions and returns numpy array for each polarization
    """
    hbar = constants.value("Planck constant over 2 pi in eV s")
    nrc = rif_constants() * ((ONEE / hbar) ** 2) * \
          np.absolute(
            (fresnel_vl(state[1], "twoe") * fresnel_lb(state[1], "twoe") *
            ((fresnel_vl(state[0], "onee") * fresnel_lb(state[0], "onee")) ** 2)) *
            reflection_components(state[0], state[1])
                    ) ** 2
    return nrc

def rif_constants():
    """
    eq. 22: constant term outside absolute value
    dependencies: none
    parents: nonlinear_reflection
    children: none
    Multiplies constants. "elecdens" merits revision.
    """
    #elecdens = 1e-28 # electronic density and scaling factor (1e-7 * 1e-21)
    elecdens = 1 # this term is included in chi^{2}
    #const = (32 * (constants.pi ** 3)) / \
    #        ((elecdens ** 2) * \
    #        ((constants.c * 100) ** 3) * \
    #        (math.cos(THETA_RAD) ** 2))
    const = 1
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
        fresnel = (2 * math.cos(THETA_RAD)) / \
            (math.cos(THETA_RAD) + wave_vector("l", energy))
    elif polarization == "p":
        fresnel = (2 * math.cos(THETA_RAD)) / \
        (epsilon("l", energy) * math.cos(THETA_RAD) + wave_vector("l", energy))
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
        fresnel = (2 * wave_vector("l", energy)) / \
            (wave_vector("l", energy) + wave_vector("b", energy))
    elif polarization == "p":
        fresnel = (2 * wave_vector("l", energy)) / \
            (epsilon("b", energy) * wave_vector("l", energy) + 
             epsilon("l", energy) * wave_vector("b", energy))
    return fresnel

def reflection_components(polar_in, polar_out):
    """
    eqs. 34a--34d
    dependencies: polarization states, 1w and 2w energy arrays
    parents: nonlinear_reflection
    children: load_shg, epsilon, wave_vector
    Calculates r_{if} factors. Loads shg arrays and does a lot of matrix
    operations. BMS program had electrostatic units multiplying shg arrays.
    """
    zzz = load_shg(VARS['zzz']) # * electrostatic_units(energy)
    zxx = load_shg(VARS['zxx']) # * electrostatic_units(energy)
    xxz = load_shg(VARS['xxz']) # * electrostatic_units(energy)
    xxx = load_shg(VARS['xxx']) # * electrostatic_units(energy)
    if polar_in == "p" and polar_out == "p":
        r_factor = math.sin(THETA_RAD) * epsilon("b", "twoe") * (((math.sin(THETA_RAD) ** 2) * (epsilon("b", "onee") ** 2) * zzz) + (wave_vector("b", "onee") ** 2) * (epsilon("l", "onee") ** 2) * zxx) + epsilon("l", "onee") * epsilon("l", "twoe") * wave_vector("b", "onee") * wave_vector("b", "twoe") * (-2 * math.sin(THETA_RAD) * epsilon("b", "onee") * xxz + wave_vector("b", "onee") * epsilon("l", "onee") * xxx * math.cos(3 * PHI_RAD))
    elif polar_in == "p" and polar_out == "s":
        r_factor = -(wave_vector("b", "onee") ** 2) * (epsilon("l", "onee") ** 2) * xxx * math.sin(3 * PHI_RAD)
    elif polar_in == "s" and polar_out == "p":
        r_factor = math.sin(THETA_RAD) * epsilon("b", "twoe") * zxx - wave_vector("b", "twoe") * epsilon("l", "twoe") * xxx * math.cos(3 * PHI_RAD)
    elif polar_in == "s" and polar_out == "s":
        r_factor = xxx * math.sin(3 * PHI_RAD)
    return r_factor

def epsilon(interface, energy):
    """
    eq. ?? REFERENCE NEEDED
    dependencies: energy array
    parents: wave_vector, fresnel_vl, fresnel_lb, reflection_components
    children: chi_spline
    Combines splines for real and imaginary parts of Chi^(1)
    """
    if interface == "l":
        chi1 = load_chi(VARS['chil'], energy)
    elif interface == "b":
        chi1 = load_chi(VARS['chib'], energy)
    #spline = chi_spline(chi1, "real", energy) + \
    #        1j * chi_spline(chi1, "imag", energy)
    #eps = 1 + (4 * constants.pi * spline)
    eps = 1 + (4 * constants.pi * chi1)
    return eps

def wave_vector(interface, energy):
    """
    eq. 5: k_{z}(omega)
    dependencies: energy array
    parents: fresnel_vl, fresnel_lb, reflection_components
    children: epsilon
    Calculates wave vector k
    """
    kz = np.sqrt(epsilon(interface, energy) - (math.sin(THETA_RAD) ** 2))
    return kz

def electrostatic_units(energy):
    """
    coefficient to convert to appropriate electrostatic units
    deprecated
    """
    area = (1 / ((2 * np.sqrt(2)) ** 2)) * 2 * np.sqrt(3)
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
    parents: none
    children: nonlinear_reflection, save_matrix
    Creates final matrix and writes to file.
    """
    nrc = np.column_stack((2*ONEE, nonlinear_reflection(["p", "p"]),
                                   nonlinear_reflection(["p", "s"]),
                                   nonlinear_reflection(["s", "p"]),
                                   nonlinear_reflection(["s", "s"])))
    #outfile = VARS['output']
    outf = sys.argv[2]
    np.savetxt(outf, nrc, fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e'),
                            delimiter='    ')

def parse_input():
    """
    parents: none
    children: none
    Parses input file specified in command line input
    """
    infile = sys.argv[1]
    params = {}
    targetfile = open(infile)
    data = targetfile.readlines()
    for line in data:
        if '#' in line:
            line = line.split('#', 1)
        else:
            key, value = line.split(":")
            params[key.strip()] = value.strip()
    targetfile.close()
    return params

def load_chi(in_file, energy):
    """
    dependencies: input file
    parents: chi_spline
    children: none
    Loads Chi^(1) file, unpacks columns, and combines into complex numpy array.
    """
    real, imag = np.loadtxt(in_file, unpack=True, usecols=[1, 2], skiprows=1)
    data = real + 1j * imag
    if energy == "onee":
        chi = data[:1000]
    elif energy == "twoe":
        chi = data[1::2]
    return chi

def load_shg(in_file):
    """
    dependencies: input file
    parents: reflection_components
    children: none
    Loads shg Chi^(2) files, unpacks columns, sums 1w and 2w for real and
    imag, and combines into complex numpy array.
    """
    real1w, imaginary1w, real2w, imaginary2w = \
            np.loadtxt(in_file, unpack=True, usecols=[1, 2, 3, 4], skiprows=1)
    real = real1w + real2w
    imaginary = imaginary1w + imaginary2w
    data = real + 1j * imaginary
    shg = data[:1000]
    return shg

def chi_spline(chi1, part, energy):
    """
    dependencies: real or imaginary indicator, energy
    parents: epsilon
    children: load_chi
    Creates spline from real part of Chi^(1) array
    """
    interpolated = \
    interpolate.InterpolatedUnivariateSpline(energy, getattr(chi1, part))
    return interpolated(energy)

### Debug functions
###
def debug():
    zzz = load_shg(VARS['zzz'])
    zxx = load_shg(VARS['zxx'])
    xxz = load_shg(VARS['xxz'])
    czzz = (math.sin(THETA_RAD) ** 3) * epsilon("b", "twoe") * (epsilon("b", "onee") ** 2) * zzz
    czxx = math.sin(THETA_RAD) * epsilon("b", "twoe") * (epsilon("b", "onee") ** 2) * (wave_vector("b", "onee") ** 2) * zxx
    cxxz = -2 * math.sin(THETA_RAD) * epsilon("b", "twoe") * (epsilon("b", "onee") ** 2) * wave_vector("b", "onee") * wave_vector("b", "twoe") * xxz
    deb = np.column_stack((ONEE, np.absolute(czzz)**2, np.absolute(czxx)**2, np.absolute(cxxz)**2))
    np.savetxt("debug/coefs.dat", deb, delimiter='    ')

def fort_comparison():
    # fort.301
    eps = np.column_stack((ONEE, np.absolute(epsilon("b", "onee")), np.absolute(epsilon("l", "onee")), np.absolute(epsilon("b", "twoe")), np.absolute(epsilon("l", "twoe"))))
    np.savetxt("debug/epsilon.dat", eps, delimiter='    ')
    # fort.302
    kz = np.column_stack((ONEE, np.absolute(wave_vector("b", "onee")), np.absolute(wave_vector("l", "onee")), np.absolute(wave_vector("b", "twoe")), np.absolute(wave_vector("l", "twoe"))))
    np.savetxt("debug/kz.dat", kz, delimiter='    ')    
    # fort.303
    fresnel1w = np.column_stack((ONEE, np.absolute(fresnel_vl("s", "onee")), np.absolute(fresnel_vl("p", "onee")), np.absolute(fresnel_lb("s", "onee")), np.absolute(fresnel_lb("p", "onee"))))
    np.savetxt("debug/fresnel1w.dat", fresnel1w, delimiter='    ')
    # fort.304
    fresnel2w = np.column_stack((ONEE, np.absolute(fresnel_vl("s", "twoe")), np.absolute(fresnel_vl("p", "twoe")), np.absolute(fresnel_lb("s", "twoe")), np.absolute(fresnel_lb("p", "twoe"))))
    np.savetxt("debug/fresnel2w.dat", fresnel2w, delimiter='    ')
    # fort.305
    ref = np.column_stack((ONEE, np.absolute(reflection_components("p", "p")), np.absolute(reflection_components("p", "s")), np.absolute(reflection_components("s", "p")), np.absolute(reflection_components("s", "s"))))
    np.savetxt("debug/refs.dat", ref, delimiter='    ')

def fort_output():
    ### epsilons
    epsl = np.column_stack((ONEE, epsilon("l", "onee").real, epsilon("l", "onee").imag))
    np.savetxt("debug/epsl.dat", epsl, delimiter='    ')
    epsb = np.column_stack((ONEE, epsilon("b", "onee").real, epsilon("b", "onee").imag))
    np.savetxt("debug/epsb.dat", epsb, delimiter='    ')
    ### chi2 components
    comps = np.column_stack((ONEE, load_shg(VARS['zzz']).real, load_shg(VARS['zzz']).imag, load_shg(VARS['zxx']).real, load_shg(VARS['zxx']).imag, load_shg(VARS['xxz']).real, load_shg(VARS['xxz']).imag, load_shg(VARS['xxx']).real, load_shg(VARS['xxx']).imag))
    np.savetxt("debug/comps.dat", comps, delimiter='    ')
    

VARS = parse_input()
fort_output()
#fort_comparison()
#debug()
control()
