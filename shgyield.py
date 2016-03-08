"""
shgyield.py is a python program designed to calculate the nonlinear reflection
coefficient for semiconductor surfaces. It works in conjunction with the
matrix elements calculated with ABINIT, an open source ab initio software, and
TINIBA, our in-house optical calculation software.

The work coded in this software can be found in an upcoming publication and is
explicitly derived in 'phd-thesis'.

Tested with Anaconda Python.

requirements:
sys, math, numpy, scipy

usage:
python shgyield.py <sample.in>
"""

import sys
import math
import numpy as np
from scipy import constants

# max energy, assumes range from 0 to 20 with 2001 steps
MAXE = 1000
ONEE = np.linspace(0.01, float(MAXE)/100, MAXE) # 1w energy array

CHI1NORM = 1.2659296143 # Normalization yo
SLABHEIGHT = 71.09425 # In Bohrs
HEIGHT2 = SLABHEIGHT*0.5 # Also in Bohrs biatch

# parses input file specified in command line input
def parse_input():
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

# loads chi1 file, unpacks columns, and creates numpy array
def epsilon(in_file, norm): 
    data = np.loadtxt(in_file, unpack=True, skiprows=1)
    real = (data[1] + data[3] + data[5])/3 # real average
    imag = (data[2] + data[4] + data[6])/3 # imag average
    coma = real + 1j * imag # complex average
    epsa = 1 + (4 * constants.pi * coma * norm) # epsilon average
    return epsa

# loads chi2 file, unpacks columns, sums columns, and creates numpy array
def shgcomp(in_file):
    data = np.loadtxt(in_file, unpack=True, skiprows=1)
    comp = (data[1] + data[3]) + 1j * (data[2] + data[4])
    shg = comp[:MAXE]
    return shg

# most generic fresnel factors
def fresnel(pol, i, j, freq):
    ki = eval("w" + i + freq)
    kj = eval("w" + j + freq)
    epsi = eval("eps" + i + freq)
    epsj = eval("eps" + j + freq)
    if pol == "p":
        factor = (2 * ki * np.sqrt(epsi * epsj)) / \
                 (ki * epsj + kj * epsi)
    elif pol == "s":
        factor = (2 * ki) / (ki + kj)
    return factor

# reads input file and establishes mode
PARAM = parse_input()
MODE = str(PARAM['mode'])

# angles
THETAIN = math.radians(float(PARAM['theta']))
PHI = math.radians(float(PARAM['phi']))

# constants, conversions, and prefactor
HBAR = constants.value("Planck constant over 2 pi in eV s")
EPS0 = constants.epsilon_0 # (F/m)
LSPEED = constants.c # (m/s)
N0E = 1.11e19 # for silicon, from PRB 81, 3781 Ref. 17 (V/m^2)
PM2TOM2 = 1e-24 # pm^2 to m^2
M2TOCM2 = 1e4 # m^2 to cm^2
TINIBASCALE = 1e6 # for scaling chi in 1e6 (pm^2/V)
SCALE = 1e20 # for R in 1e-20 (cm^2/W)
PREFACTOR = 1 / (2 * EPS0 * HBAR**2 * LSPEED**3 * math.cos(THETAIN)**2)

# loads chi2, converts to m^2/V
ZZZ = (TINIBASCALE * PM2TOM2 * shgcomp(PARAM['zzz']))
ZXX = (TINIBASCALE * PM2TOM2 * shgcomp(PARAM['zxx']))
XXZ = (TINIBASCALE * PM2TOM2 * shgcomp(PARAM['xxz']))
XXX = (TINIBASCALE * PM2TOM2 * shgcomp(PARAM['xxx']))

# creates epsilons from chi1 responses
epsl = epsilon(PARAM['chil'], CHI1NORM)
epsb = epsilon(PARAM['chib'], 1)
epsv1w = 1
epsv2w = 1
epsb1w = epsb[:MAXE]
epsb2w = epsb[1::2][:MAXE]
epsl1w = epsl[:MAXE]
epsl2w = epsl[1::2][:MAXE]

# mode switching, mostly for debugging
if MODE == "3layer": #case1
    ell1w = "l"
    ell2w = "l"
elif MODE == "2layer": #case2
    ell1w = "b"
    ell2w = "v"
elif MODE == "bulk":
    ell1w = "b"
    ell2w = "b"
elif MODE == "vacuum":
    ell1w = "v"
    ell2w = "v"
elif MODE == "hybrid":
    ell1w = "b"
    ell2w = "l"
epsl1w = eval("eps" + ell1w + "1w")
epsl2w = eval("eps" + ell2w + "2w")

# refraction indices
nv = np.sqrt(epsv1w)
Nv = np.sqrt(epsv2w)
nb = np.sqrt(epsb1w)
Nb = np.sqrt(epsb2w)
nl = np.sqrt(epsl1w)
Nl = np.sqrt(epsl2w)

# wave vectors for 1w and 2w
wv1w = np.sqrt(epsv1w - (math.sin(THETAIN) ** 2))
wv2w = np.sqrt(epsv2w - (math.sin(THETAIN) ** 2))
wb1w = np.sqrt(epsb1w - (math.sin(THETAIN) ** 2))
wb2w = np.sqrt(epsb2w - (math.sin(THETAIN) ** 2))
wl1w = np.sqrt(epsl1w - (math.sin(THETAIN) ** 2))
wl2w = np.sqrt(epsl2w - (math.sin(THETAIN) ** 2))

## fresnel factors for 1w and 2w, s and p polarizations
tvls = fresnel("s", "v", ell1w, "1w")
tvlp = fresnel("p", "v", ell1w, "1w")
tlbs = fresnel("s", ell1w, "b", "1w")
tlbp = fresnel("p", ell1w, "b", "1w")
Tvls = fresnel("s", "v", ell2w, "2w")
Tvlp = fresnel("p", "v", ell2w, "2w")
Tlbs = fresnel("s", ell2w, "b", "2w")
Tlbp = fresnel("p", ell2w, "b", "2w")


#### multiple reflections framework
delta0 = 8 * math.pi * \
((ONEE * HEIGHT2 * constants.value("Bohr radius"))/(HBAR * constants.c)) * \
np.sqrt(Nl**2 - math.sin(THETAIN)**2)

delta = 8 * math.pi * \
((ONEE * SLABHEIGHT * constants.value("Bohr radius"))/(HBAR * constants.c)) *\
np.sqrt(Nl**2 - math.sin(THETAIN)**2)

edelta0 = np.exp(1j * delta0)
edelta = np.exp(1j * delta)

Rvlp = ((Nl/Nv) * Tvlp) - 1
Rlbp = ((Nb/Nl) * Tlbp) - 1
Rvls = Tvls - 1
Rlbs = Tlbs - 1

RMp = (Rlbp * edelta0)/(1 + (Rvlp * Rlbp * edelta))
RMs = (Rlbs * edelta0)/(1 + (Rvls * Rlbs * edelta))

RMpplus = 1 + RMp
RMsplus = 1 + RMs
RMpminus = 1 - RMp
RMsminus = 1 - RMs

rMRpP = RMpplus * (math.sin(THETAIN)) * \
          ((nb**4 * math.sin(THETAIN)**2 * ZZZ) + (nl**4 * wb1w**2 * ZXX)) \
    - RMpminus * nl**2 * wb1w * wl2w * \
        ((2 * nb**2 * math.sin(THETAIN) * XXZ) 
                + (nl**2 * wb1w * XXX * math.cos(3 * PHI)))
rMRpS = -nl**4 * wb1w**2 * XXX * math.sin(3 * PHI)
rMRsP = (RMpplus * math.sin(THETAIN) * ZXX) + \
      (RMpminus * wl2w * XXX * math.cos(3 * PHI))
rMRsS = XXX * math.sin(3 * PHI)

GammaMRpP = (Tvlp/Nl) * ((tvlp * tlbp)/(nl**2 * nb))**2
GammaMRpS = Tvls * RMsplus * ((tvlp * tlbp)/(nl**2 * nb))**2
GammaMRsP = (Tvlp/Nl) * (tvls * tlbs)**2
GammaMRsS = Tvls * Tlbs * (tvls * tlbs)**2

RMRpP = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammaMRpP * rMRpP)**2
RMRpS = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammaMRpS * rMRpS)**2
RMRsP = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammaMRsP * rMRsP)**2
RMRsS = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammaMRsS * rMRsS)**2
####


# r factors for different input and output polarizations
rpP = Nb**2 * (math.sin(THETAIN)) * \
     ((nb**4 * math.sin(THETAIN)**2 * ZZZ) + (nl**4 * wb1w**2 * ZXX)) \
    - nl**2 * Nl**2 * wb1w * wb2w * \
       ((2 * nb**2 * math.sin(THETAIN) * XXZ)
         + (nl**2 * wb1w * XXX * math.cos(3 * PHI)))
rpS = - nl**4 * wb1w**2 * XXX * math.sin(3 * PHI)
rsP = (Nb**2 * math.sin(THETAIN) * ZXX)\
        + (Nl**2 * wb2w * XXX * math.cos(3 * PHI))
rsS = XXX * math.sin(3 * PHI)

# fresnel factors multiplied out for ease of debugging
GammapP = ((Tvlp * Tlbp)/(Nl**2 * Nb)) * ((tvlp * tlbp)/(nl**2 * nb))**2
GammapS = Tvls * Tlbs * ((tvlp * tlbp)/(nl**2 * nb))**2
GammasP = ((Tvlp * Tlbp)/(Nl**2 * Nb)) * (tvls * tlbs)**2
GammasS = Tvls * Tlbs * (tvls * tlbs)**2

# R factors for different input and output polarizations (in cm^2/W)
RpP = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammapP * rpP)**2
RpS = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammapS * rpS)**2
RsP = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammasP * rsP)**2
RsS = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammasS * rsS)**2

# creates columns for 2w and R factors and writes to file
NRC = np.column_stack((ONEE, RpP, RpS, RsP, RsS))
NRCMR = np.column_stack((ONEE, RMRpP, RMRpS, RMRsP, RMRsS))
OUTF = PARAM['output']
OUTFMR = PARAM['outputMR']
np.savetxt(OUTF, NRC,
           fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e'),
           delimiter='    ', header='RiF 1e-20 (cm^2/W)\n1w RpP RpS RsP RsS')
np.savetxt(OUTFMR, NRCMR,
           fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e'),
           delimiter='    ', header='RiF 1e-20 (cm^2/W)\n1w RpP RpS RsP RsS')
