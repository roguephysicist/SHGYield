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

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

## debugging variables
CHI1NORM = 1.2659296143 # Normalization yo
#SLAB = (142.1885172213904)*0.5 # In Bohrs
#SLAB = 3.7621659771810236 # In nanometers
SLAB = 100
HEIGHT2 = 30 # Also in nanometers biatch
##


#### Functions ####
def parse_input(infile): # parses input file from command line
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

def savefile(file, freq, val1, val2, val3, val4): # saves to file with header
    data = np.column_stack((freq, val1, val2, val3, val4))
    np.savetxt(file, data,
           fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e'),
           delimiter='    ', header='RiF 1e-20 (cm^2/W)\n1w RpP RpS RsP RsS')

def epsilon(in_file, norm): # loads chi1 file, creates epsilon numpy array
    data = np.loadtxt(in_file, unpack=True, skiprows=1)
    real = (data[1] + data[3] + data[5])/3      # real average
    imag = (data[2] + data[4] + data[6])/3      # imag average
    coma = real + 1j * imag                     # complex average
    epsa = 1 + (4 * constants.pi * coma * norm) # epsilon with normalization
    return epsa

def shgcomp(in_file): # loads chi2 file, sums Re and Im, creates numpy array
    data = np.loadtxt(in_file, unpack=True, skiprows=1)
    comp = (data[1] + data[3]) + 1j * (data[2] + data[4])
    shg = comp[:MAXE]
    return shg

def fresnel(pol, i, j, freq): # generic fresnel factors (see Mizrahi and Sipe)
    wi = eval("w" + i + freq)
    wj = eval("w" + j + freq)
    epsi = eval("eps" + i + freq)
    epsj = eval("eps" + j + freq)
    if pol == "p":
        factor = (2 * wi * np.sqrt(epsi * epsj)) / \
                 (wi * epsj + wj * epsi)
    elif pol == "s":
        factor = (2 * wi) / (wi + wj)
    return factor


#### Init ####
PARAM = parse_input(sys.argv[1])    # parses input file
MODE = str(PARAM['mode'])           # establishes mode


#### Energy ####
## assumes range from 0 to 20 with 2001 steps
MAXE = 1000
ONEE = np.linspace(0.01, float(MAXE)/100, MAXE) # 1w energy array


#### Constants ####
HBAR = constants.value("Planck constant over 2 pi in eV s")
PLANCK = constants.value("Planck constant in eV s")
EPS0 = constants.epsilon_0              # (F/m)
LSPEED = constants.c                    # (m/s)
N0E = 1.11e19                           # for Si, PRB81.3781 (Ref. 17) (V/m^2)
PM2TOM2 = 1e-24                         # pm^2 to m^2
M2TOCM2 = 1e4                           # m^2 to cm^2
TINIBASCALE = 1e6                       # for scaling chi2 in 1e6 (pm^2/V)
SCALE = 1e20                            # for R in 1e-20 (cm^2/W)
THETA0 = math.radians(float(PARAM['theta'])) # converts theta to radians
PHI = math.radians(float(PARAM['phi'])) # converts phi to radians
PREFACTOR = 1 / (2 * EPS0 * HBAR**2 * LSPEED**3 * math.cos(THETA0)**2)


#### Math ####
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
wv1w = np.sqrt(epsv1w - (math.sin(THETA0) ** 2))
wv2w = np.sqrt(epsv2w - (math.sin(THETA0) ** 2))
wb1w = np.sqrt(epsb1w - (math.sin(THETA0) ** 2))
wb2w = np.sqrt(epsb2w - (math.sin(THETA0) ** 2))
wl1w = np.sqrt(epsl1w - (math.sin(THETA0) ** 2))
wl2w = np.sqrt(epsl2w - (math.sin(THETA0) ** 2))

## fresnel factors for 1w and 2w, s and p polarizations
tvls = fresnel("s", "v", ell1w, "1w")
tvlp = fresnel("p", "v", ell1w, "1w")
tlbs = fresnel("s", ell1w, "b", "1w")
tlbp = fresnel("p", ell1w, "b", "1w")
Tvls = fresnel("s", "v", ell2w, "2w")
Tvlp = fresnel("p", "v", ell2w, "2w")
Tlbs = fresnel("s", ell2w, "b", "2w")
Tlbp = fresnel("p", ell2w, "b", "2w")
Rvlp = ((Nl/Nv) * Tvlp) - 1
Rlbp = ((Nb/Nl) * Tlbp) - 1
Rvls = Tvls - 1
Rlbs = Tlbs - 1

#### multiple reflections framework
lambda0 = (PLANCK * LSPEED * 1e9)/ONEE # In nanometers

delta0 = 8 * math.pi * ((ONEE * HEIGHT2 * 1e-9)/(PLANCK * LSPEED)) * \
                        np.sqrt(Nl**2 - math.sin(THETA0)**2)
delta = 8 * math.pi * ((ONEE * SLAB * 1e-9)/(PLANCK * LSPEED)) *\
                        np.sqrt(Nl**2 - math.sin(THETA0)**2)

edelta0 = np.exp(1j * delta0)
edelta = np.exp(1j * delta)

RMp = (Rlbp * edelta0)/(1 + (Rvlp * Rlbp * edelta))
RMs = (Rlbs * edelta0)/(1 + (Rvls * Rlbs * edelta))
RMpplus = 1 + RMp
RMsplus = 1 + RMs
RMpminus = 1 - RMp
RMsminus = 1 - RMs
####

# r factors for different input and output polarizations
## with multiple reflections
rMRpP = RMpplus * (math.sin(THETA0)) * \
          ((nb**4 * math.sin(THETA0)**2 * ZZZ) + (nl**4 * wb1w**2 * ZXX)) \
      - RMpminus * nl**2 * wb1w * wl2w * \
        ((2 * nb**2 * math.sin(THETA0) * XXZ) 
                + (nl**2 * wb1w * XXX * math.cos(3 * PHI)))
rMRpS = -nl**4 * wb1w**2 * XXX * math.sin(3 * PHI)
rMRsP = (RMpplus * math.sin(THETA0) * ZXX) + \
      (RMpminus * wl2w * XXX * math.cos(3 * PHI))
rMRsS = XXX * math.sin(3 * PHI)
## without multiple reflections
rpP = Nb**2 * (math.sin(THETA0)) * \
     ((nb**4 * math.sin(THETA0)**2 * ZZZ) + (nl**4 * wb1w**2 * ZXX)) \
    - nl**2 * Nl**2 * wb1w * wb2w * \
       ((2 * nb**2 * math.sin(THETA0) * XXZ)
         + (nl**2 * wb1w * XXX * math.cos(3 * PHI)))
rpS = - nl**4 * wb1w**2 * XXX * math.sin(3 * PHI)
rsP = (Nb**2 * math.sin(THETA0) * ZXX)\
        + (Nl**2 * wb2w * XXX * math.cos(3 * PHI))
rsS = XXX * math.sin(3 * PHI)

# fresnel factors multiplied out for ease of debugging
## with multiple reflections
GammaMRpP = (Tvlp/Nl) * ((tvlp * tlbp)/(nl**2 * nb))**2
GammaMRpS = Tvls * RMsplus * ((tvlp * tlbp)/(nl**2 * nb))**2
GammaMRsP = (Tvlp/Nl) * (tvls * tlbs)**2
GammaMRsS = Tvls * Tlbs * (tvls * tlbs)**2
## without multiple reflections
GammapP = ((Tvlp * Tlbp)/(Nl**2 * Nb)) * ((tvlp * tlbp)/(nl**2 * nb))**2
GammapS = Tvls * Tlbs * ((tvlp * tlbp)/(nl**2 * nb))**2
GammasP = ((Tvlp * Tlbp)/(Nl**2 * Nb)) * (tvls * tlbs)**2
GammasS = Tvls * Tlbs * (tvls * tlbs)**2

# R factors for different input and output polarizations (in cm^2/W)
#### CONSIDER FUNCTIONALIZING THIS
## with multiple reflections
RMRpP = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammaMRpP * rMRpP)**2
RMRpS = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammaMRpS * rMRpS)**2
RMRsP = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammaMRsP * rMRsP)**2
RMRsS = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammaMRsS * rMRsS)**2
## without multiple reflections
RpP = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammapP * rpP)**2
RpS = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammapS * rpS)**2
RsP = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammasP * rsP)**2
RsS = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((1/nl) * GammasS * rsS)**2

test1 = np.absolute(Rlbp)
test2 = np.absolute(Rlbp/(1 + Rvlp*Rlbp))
#### Output ####
#savefile(PARAM['output'], ONEE, RpP, RpS, RsP, RsS)
#savefile(PARAM['outputMR'], ONEE, RMRpP, RMRpS, RMRsP, RMRsS)
savefile("test", ONEE, test1.real, test1.imag, test2.real, test2.imag)
