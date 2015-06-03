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
    if energy == "1w":
        chi = data[:1000]
    elif energy == "2w":
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
    real1w, imaginary1w, real2w, imaginary2w = np.loadtxt(in_file, unpack=True, usecols=[1, 2, 3, 4], skiprows=1)
    real = real1w + real2w
    imaginary = imaginary1w + imaginary2w
    data = real + 1j * imaginary
    shg = data[:1000]
    return shg

VARS = parse_input()

onee = np.linspace(0.01, 10, 1000)
chil1w = load_chi(VARS['chil'], "1w")
chil2w = load_chi(VARS['chil'], "2w")
chib1w = load_chi(VARS['chib'], "1w")
chib2w = load_chi(VARS['chib'], "2w")
zzz = load_shg(VARS['zzz']) # * electrostatic_units(energy)
zxx = load_shg(VARS['zxx']) # * electrostatic_units(energy)
xxz = load_shg(VARS['xxz']) # * electrostatic_units(energy)
xxx = load_shg(VARS['xxx']) # * electrostatic_units(energy)
hbar = constants.value("Planck constant over 2 pi in eV s")
#elecdens = 1e-28 # electronic density and scaling factor (1e-7 * 1e-21)
elecdens = 1 # this term is included in chi^{2}
#const = (32 * (constants.pi ** 3)) / ((elecdens ** 2) * ((constants.c * 100) ** 3) * (math.cos(THETA_RAD) ** 2))
const = 1
epsl1w = 1 + (4 * constants.pi * chil1w)
epsl2w = 1 + (4 * constants.pi * chil2w)
epsb1w = 1 + (4 * constants.pi * chib1w)
epsb2w = 1 + (4 * constants.pi * chib2w)
wvl1w = np.sqrt(epsl1w - (math.sin(THETA_RAD) ** 2))
wvl2w = np.sqrt(epsl2w - (math.sin(THETA_RAD) ** 2))
wvb1w = np.sqrt(epsb1w - (math.sin(THETA_RAD) ** 2))
wvb2w = np.sqrt(epsb2w - (math.sin(THETA_RAD) ** 2))
tvls = (2 * math.cos(THETA_RAD)) / (math.cos(THETA_RAD) + wvl1w)
Tvls = (2 * math.cos(THETA_RAD)) / (math.cos(THETA_RAD) + wvl2w)
tvlp = (2 * math.cos(THETA_RAD)) / (epsl1w * math.cos(THETA_RAD) + wvl1w)
Tvlp = (2 * math.cos(THETA_RAD)) / (epsl2w * math.cos(THETA_RAD) + wvl2w)
tlbs = (2 * wvl1w) / (wvl1w + wvb1w)
Tlbs = (2 * wvl2w) / (wvl2w + wvb2w)
tlbp = (2 * wvl1w) / (epsb1w * wvl1w + epsl1w * wvb1w)
Tlbp = (2 * wvl2w) / (epsb2w * wvl2w + epsl2w * wvb2w)
rpp = math.sin(THETA_RAD) * epsb2w * (((math.sin(THETA_RAD) ** 2) * (epsb1w ** 2) * zzz) + (wvb1w ** 2) * (epsl1w ** 2) * zxx) + epsl1w * epsl2w * wvb1w * wvb2w * (-2 * math.sin(THETA_RAD) * epsb1w * xxz + wvb1w * epsl1w * xxx * math.cos(3 * PHI_RAD))
rps = -(wvb1w ** 2) * (epsl1w ** 2) * xxx * math.sin(3 * PHI_RAD)
rsp = math.sin(THETA_RAD) * epsb2w * zxx - wvb2w * epsl2w * xxx * math.cos(3 * PHI_RAD)
rss = xxx * math.sin(3 * PHI_RAD)
Rpp = const * ((onee / hbar) ** 2) * np.absolute((Tvlp * Tlbp * ((tvlp * tlbp) ** 2)) * rpp) ** 2
Rps = const * ((onee / hbar) ** 2) * np.absolute((Tvls * Tlbs * ((tvlp * tlbp) ** 2)) * rps) ** 2
Rsp = const * ((onee / hbar) ** 2) * np.absolute((Tvlp * Tlbp * ((tvls * tlbs) ** 2)) * rsp) ** 2
Rss = const * ((onee / hbar) ** 2) * np.absolute((Tvls * Tlbs * ((tvls * tlbs) ** 2)) * rss) ** 2
nrc = np.column_stack((2*onee, Rpp, Rps, Rsp, Rss))
#outfile = VARS['output']
outf = sys.argv[2]
np.savetxt(outf, nrc, fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e'), delimiter='    ')    
