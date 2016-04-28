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

## some variables
CHI1NORM = 1.2659296143 # Normalization yo
#SIGMA = 4.266666666 # broadening (sigma = 0.10 eV)
SIGMA = 3.19999995 # broadening (sigma = 0.075 eV)
#SIGMA = 2.133333333 # broadening (sigma = 0.05 eV)


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

def fresnel(kind, i, j, pol, freq): # generic fresnel factors (see Mizrahi)
    wi = eval("w" + i + freq)
    wj = eval("w" + j + freq)
    epsi = eval("eps" + i + freq)
    epsj = eval("eps" + j + freq)
    if kind == "t":
        if pol == "p":
            factor = (2 * wi * np.sqrt(epsi * epsj))/(wi * epsj + wj * epsi)
        elif pol == "s":
            factor = (2 * wi)/(wi + wj)
    elif kind == "r":
        if pol == "p":
            factor = ((wi * epsj) - (wj * epsi))/((wi * epsj) + (wj * epsi))
        elif pol == "s":
            factor =  (wi - wj)/(wi + wj)
    return factor

def shgyield(gamma, riF): # function for the final yield
    RiF = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
                    np.absolute((1/nl) * gamma * riF)**2
    broadened = broad(RiF, SIGMA)
    return broadened


#### Init ####
PARAM = parse_input(sys.argv[1]) # parses input file
MODE = str(PARAM['mode']) # establishes mode
MULTIREF = str(PARAM['multiref']) # if multiple reflections are considered

#### Energy ####
## assumes range from 0 to 20 with 2001 steps
MAXE = 1000
ONEE = np.linspace(0.01, float(MAXE)/100, MAXE) # 1w energy array


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
THICKNESS = float(PARAM['thickness'])
D2 = float(PARAM['d2'])
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
tvls = fresnel("t", "v", ell1w, "s", "1w")
tvlp = fresnel("t", "v", ell1w, "p", "1w")
tlbs = fresnel("t", ell1w, "b", "s", "1w")
tlbp = fresnel("t", ell1w, "b", "p", "1w")
Tvls = fresnel("t", "v", ell2w, "s", "2w")
Tvlp = fresnel("t", "v", ell2w, "p", "2w")
Tlbs = fresnel("t", ell2w, "b", "s", "2w")
Tlbp = fresnel("t", ell2w, "b", "p", "2w")
rvlp = fresnel("r", "v", ell1w, "p", "1w")
rvls = fresnel("r", "v", ell1w, "s", "1w")
rlbp = fresnel("r", ell1w, "b", "p", "1w")
rlbs = fresnel("r", ell1w, "b", "s", "1w")
Rvlp = fresnel("r", "v", ell2w, "p", "2w")
Rvls = fresnel("r", "v", ell2w, "s", "2w")
Rlbp = fresnel("r", ell2w, "b", "p", "2w")
Rlbs = fresnel("r", ell2w, "b", "s", "2w")

#### multiple reflections framework
varphi = 4 * math.pi * ((ONEE * THICKNESS * 1e-9)/(PLANCK * LSPEED)) * wl1w
delta = 8 * math.pi * ((ONEE * THICKNESS * 1e-9)/(PLANCK * LSPEED)) * wl2w
delta0 = 8 * math.pi * ((ONEE * D2 * 1e-9)/(PLANCK * LSPEED)) * wl2w

if MULTIREF == "yes":
    rMp = (rlbp * np.exp(1j * varphi))/(1 + (rvlp * rlbp * np.exp(1j * varphi)))
    rMs = (rlbs * np.exp(1j * varphi))/(1 + (rvls * rlbs * np.exp(1j * varphi)))
    RMp = (Rlbp * np.exp(1j * delta0))/(1 + (Rvlp * Rlbp * np.exp(1j * delta)))
    RMs = (Rlbs * np.exp(1j * delta0))/(1 + (Rvls * Rlbs * np.exp(1j * delta)))
    RMpav = (Rlbp * np.exp(1j * delta/2))/\
            (1 + (Rvlp * Rlbp * np.exp(1j * delta)))\
            * np.sin(delta/2)/(delta/2)
    RMsav = (Rlbs * np.exp(1j * delta/2))/\
            (1 + (Rvls * Rlbs * np.exp(1j * delta)))\
            * np.sin(delta/2)/(delta/2)
elif MULTIREF == "no":
    rMp = rlbp
    rMs = rlbs
    RMp = Rlbp
    RMs = Rlbs
    RMpav = Rlbp
    RMsav = Rlbs

rMplusp = 1 + rMp
rMpluss = 1 + rMs
RMplusp = 1 + RMp
RMpluss = 1 + RMs
rMminusp = 1 - rMp
rMminuss = 1 - rMs
RMminusp = 1 - RMp
RMminuss = 1 - RMs
####

# r factors for different input and output polarizations
rMRpP = RMplusp * (math.sin(THETA0)) * \
            ((rMplusp**2 * math.sin(THETA0)**2 * ZZZ) \
           + (rMminusp**2 * wl1w**2 * ZXX)) \
      - RMminusp * wl1w * wl2w * \
            ((2 * rMplusp * rMminusp * math.sin(THETA0) * XXZ) 
           + (rMminusp**2 * wl1w * XXX * math.cos(3 * PHI)))
rMRpS = -rMminusp * wl1w**2 * XXX * math.sin(3 * PHI)
rMRsP = (RMplusp * math.sin(THETA0) * ZXX) + \
      (RMminusp * wl2w * XXX * math.cos(3 * PHI))
rMRsS = XXX * math.sin(3 * PHI)

# fresnel factors multiplied out for ease of debugging
GammaMRpP = (Tvlp/Nl) * (tvlp/nl)**2
GammaMRpS = Tvls * RMpluss * (tvlp/nl)**2
GammaMRsP = (Tvlp/Nl) * (tvls * rMpluss)**2
GammaMRsS = Tvls * RMpluss * (tvls * rMpluss)**2

# R factors for different input and output polarizations (in cm^2/W)
RMRpP = shgyield(GammaMRpP, rMRpP)
RMRpS = shgyield(GammaMRpS, rMRpS)
RMRsP = shgyield(GammaMRsP, rMRsP)
RMRsS = shgyield(GammaMRsS, rMRsS)


#### Output ####
savefile(PARAM['output'], ONEE, RMRpP, RMRpS, RMRsP, RMRsS)
