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

def load_chiz(in_file, energy):
    """
    Loads Chi^(1) file, unpacks columns, and combines into complex numpy array.
    NEEDS TO BE OPTIMIZED
    """
    global MAX_E
    real, imag = np.loadtxt(in_file, unpack=True, usecols=[5, 6], skiprows=1)
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

# for screening with eps_z ## NEEDS TO BE OPTIMIZED
chil1wz = load_chiz(param['chil'], "1w")
epsl1wz = 1 + (4 * constants.pi * chil1wz)

# constants 
hbar = constants.value("Planck constant over 2 pi in eV s")
eps0 = constants.epsilon_0 / 100 # (F/cm)
lspeed = constants.c * 100 # (cm/s)
rydberg = constants.value("Rydberg constant times hc in eV") # (eV)
a0 = constants.value("lattice parameter of silicon") * 100 # (cm)
ab = constants.value("Bohr radius") * 100 # (cm)
#n0e = 1.11e19 # for silicon, from PRB 81, 3781 Ref. 17 (V/m^2)
n0esquared = 1.4399764 * 1e6 * 1e-13 * constants.e * (1.5e10)**2 # (J/cm^5)
n0e_pm = 1.11e-5 # for silicon, from PRB 81, 3781 Ref. 17 (V/pm^2)

# unit conversions and prefactors
onee = np.linspace(0.01, 10, MAX_E) # 1w energy array
pico2cent2 = 1e-20 # pm^2 to cm^2
tinibascale = 1e6
intensity = 1e7 # intensity from CGS (erg/cm^2 s) to MKS (cm^2/W)
scale = 1e21 # for R in 1e-21 (cm^2/W) for \chi^ijk in (esu cm)
mks2cgs = 1e-12 # when converting from MKS [1e6 pm^2/V] to CGS [1e-13 esu cm]
area = 2 * 3**0.5 / 8
esufactor = 1j * (2 * rydberg)**5 * (ab/a0)**5 * 2.08e-15 * (a0 / 1e-8)**3 / area
#prefactor = (32 * constants.pi**3) / (hbar**2 * lspeed**3 * math.cos(THETA_RAD)**2) * intensity
#prefactor = (32 * constants.pi**3) / (hbar**2 * n0esquared * lspeed**3 * math.cos(THETA_RAD)**2)
#prefactor = (32 * (constants.pi ** 3) * ((onee / hbar) ** 2)) / (eps0 * (lspeed ** 3) * (math.cos(THETA_RAD) ** 2))
prefactor = 1

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
zzz = load_shg(param['zzz'])/(epsl1wz**2)
zxx = load_shg(param['zxx'])
xxz = load_shg(param['xxz'])/epsl1wz
xxx = load_shg(param['xxx'])

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
Rpp = (onee ** 2) * np.absolute(fpp * rpp)**2
Rps = (onee ** 2) * np.absolute(fps * rps)**2
Rsp = (onee ** 2) * np.absolute(fsp * rsp)**2
Rss = (onee ** 2) * np.absolute(fss * rss)**2

# creates columns for 2w and R factors and writes to file
nrc = np.column_stack((2*onee, Rpp, Rps, Rsp, Rss))
eps = np.column_stack((onee, epsl1w.real, epsl1w.imag, epsl1wz.real, epsl1wz.imag))
chi = np.column_stack((onee, xxx.real, xxx.imag, xxz.real, xxz.imag, zxx.real, zxx.imag, zzz.real, zzz.imag))
rif = np.column_stack((onee, rpp.real, rpp.imag, rps.real, rps.imag, rsp.real, rsp.imag, rss.real, rss.imag))
fre = np.column_stack((onee, fpp.real, fpp.imag, fps.real, fps.imag, fsp.real, fsp.imag, fss.real, fss.imag))
outf = param['output']
# outf = sys.argv[2]
np.savetxt(outf, nrc, fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e'), delimiter='    ')
np.savetxt('eps_' + outf, eps, fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e'), delimiter='    ')
np.savetxt('chi_' + outf, chi, fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e'), delimiter='    ')
np.savetxt('rif_' + outf, rif, fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e'), delimiter='    ')
np.savetxt('fre_' + outf, fre, fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e'), delimiter='    ')
