#!/Users/sma/anaconda/bin/python
"""
nrc.py is a python program designed to calculate the Nonlinear reflection
coefficient for silicon surfaces. It works in conjunction with the matrix
elements calculated with ABINIT, an open source ab initio software,
and TINIBA, our in-house optical calculation software.

The work codified in this software can be found in Phys.Rev.B66, 195329(2002).
Please refer to "Strain induced SHG" manuscript (BMS) for equation references
unless stated otherwise.

Experimental gap for bulk Si = 3.4 according to Landolt-Boernstein
"""

import sys
import math
import numpy as np
from scipy import constants

# Angles and energies
THETA_RAD = math.radians(65)
PHI_RAD = math.radians(30)

def parse_input():
    """
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
    Loads Chi^(1) file, unpacks columns, and combines into complex numpy array.
    """
    global MAX_E
    real, imag = np.loadtxt(in_file, unpack=True, usecols=[1, 2], skiprows=1)
    data = real + 1j * imag
    MAX_E = len(data)/2
    if energy == "1w":
        chi = data[:MAX_E]
    elif energy == "2w":
        chi = data[1::2]
    return chi

def load_shg(in_file):
    """
    Loads shg Chi^(2) files, unpacks columns, sums 1w and 2w for real and
    imag, and combines into complex numpy array.
    """
    real1w, imaginary1w, real2w, imaginary2w = np.loadtxt(in_file,
                                                          unpack=True,
                                                          usecols=[1, 2, 3, 4],
                                                          skiprows=1)
    real = real1w + real2w
    imaginary = imaginary1w + imaginary2w
    data = real + 1j * imaginary
    shg = data[:MAX_E]
    return shg

# reads input file
param = parse_input()

# loads chi1 and epsilons
chil1w = load_chi(param['chil'], "1w")
chil2w = load_chi(param['chil'], "2w")
chib1w = load_chi(param['chib'], "1w")
chib2w = load_chi(param['chib'], "2w")
epsl1w = 1 + (4 * constants.pi * chil1w)
epsl2w = 1 + (4 * constants.pi * chil2w)
epsb1w = 1 + (4 * constants.pi * chib1w)
epsb2w = 1 + (4 * constants.pi * chib2w)

# constants and numpy array for 1w energy values
onee = np.linspace(0.01, 10, MAX_E)
hbar = constants.value("Planck constant over 2 pi in eV s")
eps0 = constants.epsilon_0 / 100 # (in F/cm)
lspeed = constants.c * 100 # (in cm/s)
pico2cent2 = 1e-20
n0e = 1.11e19 # for silicon (ref. 17 from PRB 81, 3781. not needed in mks!)
const = (32 * (constants.pi ** 3) * ((onee / hbar) ** 2)) / \
        (eps0 * (lspeed ** 3) * (math.cos(THETA_RAD) ** 2))

# wave vectors for 1w and 2w
kzl1w = np.sqrt(epsl1w - (math.sin(THETA_RAD) ** 2))
kzl2w = np.sqrt(epsl2w - (math.sin(THETA_RAD) ** 2))
kzb1w = np.sqrt(epsb1w - (math.sin(THETA_RAD) ** 2))
kzb2w = np.sqrt(epsb2w - (math.sin(THETA_RAD) ** 2))

# fresnel factors for 1w and 2w, s and p polarizations
tvls = (2 * math.cos(THETA_RAD)) / (math.cos(THETA_RAD) + kzl1w)
tvlp = (2 * math.cos(THETA_RAD)) / (epsl1w * math.cos(THETA_RAD) + kzl1w)
tlbs = (2 * kzl1w) / (kzl1w + kzb1w)
tlbp = (2 * kzl1w) / (epsb1w * kzl1w + epsl1w * kzb1w)
Tvls = (2 * math.cos(THETA_RAD)) / (math.cos(THETA_RAD) + kzl2w)
Tvlp = (2 * math.cos(THETA_RAD)) / (epsl2w * math.cos(THETA_RAD) + kzl2w)
Tlbs = (2 * kzl2w) / (kzl2w + kzb2w)
Tlbp = (2 * kzl2w) / (epsb2w * kzl2w + epsl2w * kzb2w)

# loads chi2, converts to cm^2/V, and screens them with layer epsilon
zzz = pico2cent2 * load_shg(param['zzz'])/epsl1w ** 2
zxx = pico2cent2 * load_shg(param['zxx'])
xxz = pico2cent2 * load_shg(param['xxz'])/epsl1w
xxx = pico2cent2 * load_shg(param['xxx'])

# r factors for different input and output polarizations
rpp = math.sin(THETA_RAD) * epsb2w * \
      (((math.sin(THETA_RAD) ** 2) * (epsb1w ** 2) * zzz) + \
      (kzb1w ** 2) * (epsl1w ** 2) * zxx) + epsl1w * epsl2w * \
      kzb1w * kzb2w * (-2 * math.sin(THETA_RAD) * epsb1w * xxz + \
      kzb1w * epsl1w * xxx * math.cos(3 * PHI_RAD))
rps = -(kzb1w ** 2) * (epsl1w ** 2) * xxx * math.sin(3 * PHI_RAD)
rsp = math.sin(THETA_RAD) * epsb2w * zxx - \
      kzb2w * epsl2w * xxx * math.cos(3 * PHI_RAD)
rss = xxx * math.sin(3 * PHI_RAD)

# R factors for different input and output polarizations (in cm^2/W)
Rpp = const * np.absolute((Tvlp * Tlbp * ((tvlp * tlbp) ** 2)) * rpp) ** 2
Rps = const * np.absolute((Tvls * Tlbs * ((tvlp * tlbp) ** 2)) * rps) ** 2
Rsp = const * np.absolute((Tvlp * Tlbp * ((tvls * tlbs) ** 2)) * rsp) ** 2
Rss = const * np.absolute((Tvls * Tlbs * ((tvls * tlbs) ** 2)) * rss) ** 2

# creates columns for 2w and R factors and writes to file
nrc = np.column_stack((2*onee, Rpp, Rps, Rsp, Rss))
outf = param['output']
#outf = sys.argv[2]
np.savetxt(outf, nrc, fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e'),
                      delimiter='    ')
