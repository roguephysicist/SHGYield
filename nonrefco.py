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

For gnuplot:
abso(w,x,y,z)=sqrt((w+y)**2+(x+z)**2)
"""

import sys
import math
import numpy as np
from scipy import constants

# Angles and max energy
THETA_RAD = math.radians(65)
PHI_RAD = math.radians(30)
MAXE = 1000

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

def load_chi(in_file):
    """
    Loads Chi^(1) file, unpacks columns, and creates numpy array.
    """
    data = np.loadtxt(in_file, unpack=True, skiprows=1)
    return data

def load_shg(in_file):
    """
    Loads shg Chi^(2) files, unpacks columns, sums 1w and 2w for real and
    imag, and combines into complex numpy array.
    """
    data = np.loadtxt(in_file, unpack=True, skiprows=1)
    comp = (data[1] + data[3]) + 1j * (data[2] + data[4])
    shg = comp[:MAXE]
    return shg

# reads input file
param = parse_input()

# constants, conversions, and prefactor
onee = np.linspace(0.01, 10, MAXE) # 1w energy array
hbar = constants.value("Planck constant over 2 pi in eV s")
eps0 = constants.epsilon_0 # (F/m)
lspeed = constants.c # (m/s)
n0e = 1.11e19 # for silicon, from PRB 81, 3781 Ref. 17 (V/m^2)
pm2tom2 = 1e-24 # pm^2 to m^2
m2tocm2 = 1e4 # m^2 to cm^2
tinibascale = 1e6 # for scaling chi in 1e6 (pm^2/V)
scale = 1e20 # for R in 1e-20 (cm^2/W)
prefactor = 1 / (2 * eps0 * hbar**2 * lspeed**3 * math.cos(THETA_RAD)**2)

# loads chi1 and epsilons
chil = load_chi(param['chil'])
chib = load_chi(param['chib'])
chilzz = chil[5] + 1j * chil[6]
chilxx = chil[1] + 1j * chil[2]
chibxx = chib[1] + 1j * chib[2]
epslzz = 1 + (4 * constants.pi * chilzz[:MAXE])
epsl1w = 1 + (4 * constants.pi * chilxx[:MAXE])
epsl2w = 1 + (4 * constants.pi * chilxx[1::2])
epsb1w = 1 + (4 * constants.pi * chibxx[:MAXE])
epsb2w = 1 + (4 * constants.pi * chibxx[1::2])

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
zzz = (tinibascale * pm2tom2 * load_shg(param['zzz']))/(epslzz**2)
zxx = (tinibascale * pm2tom2 * load_shg(param['zxx']))
xxz = (tinibascale * pm2tom2 * load_shg(param['xxz']))/epslzz
xxx = (tinibascale * pm2tom2 * load_shg(param['xxx']))

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

# fresnel factors multiplied out for ease of debugging
fpp = Tvlp * Tlbp * (tvlp * tlbp)**2
fps = Tvls * Tlbs * (tvlp * tlbp)**2
fsp = Tvlp * Tlbp * (tvls * tlbs)**2
fss = Tvls * Tlbs * (tvls * tlbs)**2

# R factors for different input and output polarizations (in cm^2/W)
Rpp = scale * m2tocm2 * prefactor * (onee ** 2) * np.absolute(fpp * rpp)**2
Rps = scale * m2tocm2 * prefactor * (onee ** 2) * np.absolute(fps * rps)**2
Rsp = scale * m2tocm2 * prefactor * (onee ** 2) * np.absolute(fsp * rsp)**2
Rss = scale * m2tocm2 * prefactor * (onee ** 2) * np.absolute(fss * rss)**2

# creates columns for 2w and R factors and writes to file
nrc = np.column_stack((2*onee, Rpp, Rps, Rsp, Rss))
outf = param['output']
# outf = sys.argv[2]
np.savetxt(outf, nrc, fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e'),
delimiter='    ',
header='RiF in 1e-20 (cm^2/W)\n\
2w     Rpp' + 21*' ' + 'Rps' + 21*' ' + 'Rsp' + 21*' ' + 'Rss')
