"""
shgyield.py is a python program designed to calculate the nonlinear reflection
coefficient for semiconductor surfaces. It works in conjunction with the matrix
elements calculated with ABINIT, an open source ab initio software, and TINIBA,
our in-house optical calculation software.

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
from scipy import constants, ndimage

## debugging variables
CHI1NORM = 1.2659296143 # Normalization yo
THICKNESS = 3.7621659771810236 # In nanometers
D2 = 3.7621659771810236 # Also in nanometers biatch


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

def broad(target, sigma): 
    data = ndimage.filters.gaussian_filter(target, sigma)
    return data

def savefile(file, freq, val1, val2, val3, val4): # saves to file with header
    data = np.column_stack((freq, val1, val2, val3, val4))
    np.savetxt(file, data,
            fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e'),
            delimiter='    ',
            header='RiF 1e-20 (cm^2/W)\n1w(eV) RpP'+21*" "+\
                   'RpS'+21*" "+'RsP'+21*" "+'RsS')

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

def fresnel(i, j, pol, freq): # generic fresnel factors (see Mizrahi)
    wi = eval("w" + i + freq)
    wj = eval("w" + j + freq)
    epsi = eval("eps" + i + freq)
    epsj = eval("eps" + j + freq)
    if pol == "p":
        factor = (2 * wi * np.sqrt(epsi * epsj))/(wi * epsj + wj * epsi)
    elif pol == "s":
        factor = (2 * wi)/(wi + wj)
    return factor

def shgyield(gamma, riF): # function for the final yield
    RiF = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
                    np.absolute((1/nl) * gamma * riF)**2
    broadened = broad(RiF, SIGMA)
    return broadened


#### Init ####
PARAM = parse_input(sys.argv[1]) # parses input file
MODE = str(PARAM['mode']) # establishes mode


#### Energy ####
## assumes range from 0 to 20 with 2001 steps
MAXE = 1000
ONEE = np.linspace(0.01, float(MAXE)/100, MAXE) # 1w energy array
SIGMA = 4.266666666 # broadening (sigma = 0.10 eV)


#### Constants ####
HBAR = constants.value("Planck constant over 2 pi in eV s")
PLANCK = constants.value("Planck constant in eV s")
EPS0 = constants.epsilon_0 # (F/m)
LSPEED = constants.c # (m/s)
N0E = 1.11e19 # for Si, PRB 81.3781 (Ref. 17) (V/m^2)
PM2TOM2 = 1e-24 # pm^2 to m^2
M2TOCM2 = 1e4 # m^2 to cm^2
TINIBASCALE = 1e6 # for scaling chi2 in 1e6 (pm^2/V)
SCALE = 1e20 # for R in 1e-20 (cm^2/W)
THETA0 = math.radians(float(PARAM['theta'])) # converts theta to radians
PHI = math.radians(float(PARAM['phi'])) # converts phi to radians
LAMBDA0 = (PLANCK * LSPEED * 1e9)/ONEE # In nanometers
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
wv1w = np.sqrt(epsv1w - (math.sin(THETA0)**2))
wv2w = np.sqrt(epsv2w - (math.sin(THETA0)**2))
wb1w = np.sqrt(epsb1w - (math.sin(THETA0)**2))
wb2w = np.sqrt(epsb2w - (math.sin(THETA0)**2))
wl1w = np.sqrt(epsl1w - (math.sin(THETA0)**2))
wl2w = np.sqrt(epsl2w - (math.sin(THETA0)**2))

## fresnel factors for 1w and 2w, s and p polarizations
tvls = fresnel("v", ell1w, "s", "1w")
tvlp = fresnel("v", ell1w, "p", "1w")
tlbs = fresnel(ell1w, "b", "s", "1w")
tlbp = fresnel(ell1w, "b", "p", "1w")
Tvls = fresnel("v", ell2w, "s", "2w")
Tvlp = fresnel("v", ell2w, "p", "2w")
Tlbs = fresnel(ell2w, "b", "s", "2w")
Tlbp = fresnel(ell2w, "b", "p", "2w")

# r factors for different input and output polarizations
rpP = Nb**2 * (math.sin(THETA0)) * \
        ((nb**4 * math.sin(THETA0)**2 * ZZZ) + (nl**4 * wb1w**2 * ZXX)) \
    - nl**2 * Nl**2 * wb1w * wb2w * \
        ((2 * nb**2 * math.sin(THETA0) * XXZ)
       + (nl**2 * wb1w * XXX * math.cos(3 * PHI)))
rpS = -nl**4 * wb1w**2 * XXX * math.sin(3 * PHI)
rsP = (Nb**2 * math.sin(THETA0) * ZXX)\
        + (Nl**2 * wb2w * XXX * math.cos(3 * PHI))
rsS = XXX * math.sin(3 * PHI)

# fresnel factors multiplied out for ease of debugging
GammapP = ((Tvlp * Tlbp)/(Nl**2 * Nb)) * ((tvlp * tlbp)/(nl**2 * nb))**2
GammapS = Tvls * Tlbs * ((tvlp * tlbp)/(nl**2 * nb))**2
GammasP = ((Tvlp * Tlbp)/(Nl**2 * Nb)) * (tvls * tlbs)**2
GammasS = Tvls * Tlbs * (tvls * tlbs)**2

# R factors for different input and output polarizations (in cm^2/W)
RpP = shgyield(GammapP, rpP)
RpS = shgyield(GammapS, rpS)
RsP = shgyield(GammasP, rsP)
RsS = shgyield(GammasS, rsS)


#### Output ####
savefile(PARAM['output'], ONEE, RpP, RpS, RsP, RsS)
