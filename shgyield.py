"""
shgyield.py is a python program designed to calculate the nonlinear reflection
coefficient for semiconductor surfaces. It works in conjunction with the matrix
elements calculated with ABINIT, an open source ab initio software, and TINIBA,
our in-house optical calculation software.

Tested with Anaconda Python 4.0.0.

requirements:
sys, math, numpy, scipy

usage:
python shgyield.py <sample.in>
"""

# import sys
import math
import numpy as np
from scipy import constants, ndimage


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
    value = sigma * 42.666666666
    data = ndimage.filters.gaussian_filter(target, value)
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

def layer_depth(in_file):
    data = np.loadtxt(in_file, dtype=float, unpack=True, usecols=(1,), skiprows=1)
    depth = (data[0] - data[data>=0]) * constants.value("Bohr radius") * 1e9 # nanometers
    return depth


#### Input section ####
MODE = "3-layer"
theta = 65
phi = 30
norm = 1.2659296143
sigma = 0.075
# Multiple reflections (values in nanometers)
MULTIREF = "yes"
thickness = 10

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
THETA0 = math.radians(float(theta)) # converts theta to radians
PHI = math.radians(float(phi)) # converts phi to radians
THICKNESS = float(thickness) # thickness of the thin layer for in nanometers multiref
# DEPTH = str(depth) # depth at which we place the polarization sheet
SIGMA = float(sigma) # standard deviation for gaussian broadening
CHI1NORM = float(norm) # DEBUG: normalization factor for chi1
LAMBDA0 = (PLANCK * LSPEED * 1e9)/ONEE # In nanometers
PREFACTOR = 1 / (2 * EPS0 * HBAR**2 * LSPEED**3 * math.cos(THETA0)**2)


#### Math ####

# creates epsilons from chi1 responses
chib = "res/chi1.sm_0.075_xx_yy_zz_3107_25-nospin_scissor_0.70_Nc_26"
chil = "res/calChi1.sm_0.075_xx_yy_zz_576_half-slab_12-nospin_scissor_0.70_Nc_103"
epsl = epsilon(chil, CHI1NORM)
epsb = epsilon(chib, 1)
epsv1w = 1
epsv2w = 1
epsb1w = epsb[:MAXE]
epsb2w = epsb[1::2][:MAXE]
epsl1w = epsl[:MAXE]
epsl2w = epsl[1::2][:MAXE]

# mode switching, mostly for debugging
if MODE == "3-layer":
    ell1w = "l"
    ell2w = "l"
elif MODE == "2-layer-fresnel":
    ell1w = "b"
    ell2w = "v"
elif MODE == "2-layer-bulk":
    ell1w = "b"
    ell2w = "b"
elif MODE == "2-layer-vacuum":
    ell1w = "v"
    ell2w = "v"
elif MODE == "3-layer-hybrid":
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


sumGammaMRpP = 0
sumGammaMRpS = 0
sumGammaMRsP = 0
sumGammaMRsS = 0
sumrMRpP = 0
sumrMRpS = 0
sumrMRsP = 0
sumrMRsS = 0

LAYERS = layer_depth("front.layers.xy")

for idx, val in enumerate(LAYERS, start=1):
# for idx in range(1, 26):
    
    DEPTH = float(val)
    # DEPTH = "average"
    # DEPTH = LAYERS[idx-1]
    
    xxx = "res-layers/shgC.sm_0.075_xxx_576_" + str(idx) + "_12-nospin_scissor_0.70_Nc_103"
    xxz = "res-layers/shgC.sm_0.075_xxz_576_" + str(idx) + "_12-nospin_scissor_0.70_Nc_103"
    zxx = "res-layers/shgC.sm_0.075_zxx_576_" + str(idx) + "_12-nospin_scissor_0.70_Nc_103"
    zzz = "res-layers/shgC.sm_0.075_zzz_576_" + str(idx) + "_12-nospin_scissor_0.70_Nc_103"
    
    XXX = (TINIBASCALE * PM2TOM2 * shgcomp(xxx))
    XYY = - XXX
    YYX = - XXX
    XXZ = (TINIBASCALE * PM2TOM2 * shgcomp(xxz))
    YYZ = XXZ
    ZXX = (TINIBASCALE * PM2TOM2 * shgcomp(zxx))
    ZYY = ZXX
    ZZZ = (TINIBASCALE * PM2TOM2 * shgcomp(zzz))
    XXY = 0
    XYZ = 0
    XZZ = 0
    YXX = 0
    YXZ = 0
    YYY = 0
    YZZ = 0
    ZXY = 0
    ZXZ = 0
    ZZY = 0
    
    #### multiple reflections framework
    if MULTIREF == "yes":
        varphi = 4 * math.pi * ((ONEE * THICKNESS * 1e-9)/(PLANCK * LSPEED)) * wl1w
        delta = 8 * math.pi * ((ONEE * THICKNESS * 1e-9)/(PLANCK * LSPEED)) * wl2w
        if DEPTH == "average":
            RMpav = (Rlbp * np.exp(1j * delta/2))/(1 + (Rvlp * Rlbp * np.exp(1j * delta))) * np.sin(delta/2)/(delta/2)
            RMsav = (Rlbs * np.exp(1j * delta/2))/(1 + (Rvls * Rlbs * np.exp(1j * delta))) * np.sin(delta/2)/(delta/2)
            RMplusp = 1 + RMpav
            RMpluss = 1 + RMsav
            RMminusp = 1 - RMpav
            RMminuss = 1 - RMsav
        else:
            D2 = float(DEPTH)
            delta0 = 8 * math.pi * ((ONEE * D2 * 1e-9)/(PLANCK * LSPEED)) * wl2w
            RMp = (Rlbp * np.exp(1j * delta0))/(1 + (Rvlp * Rlbp * np.exp(1j * delta)))
            RMs = (Rlbs * np.exp(1j * delta0))/(1 + (Rvls * Rlbs * np.exp(1j * delta)))
            RMplusp = 1 + RMp
            RMpluss = 1 + RMs
            RMminusp = 1 - RMp
            RMminuss = 1 - RMs    
        rMp = (rlbp * np.exp(1j * varphi))/(1 + (rvlp * rlbp * np.exp(1j * varphi)))
        rMs = (rlbs * np.exp(1j * varphi))/(1 + (rvls * rlbs * np.exp(1j * varphi)))
        rMplusp = 1 + rMp
        rMpluss = 1 + rMs
        rMminusp = 1 - rMp
        rMminuss = 1 - rMs
    elif MULTIREF == "no":
        RMplusp = 1 + Rlbp
        RMpluss = 1 + Rlbs
        RMminusp = 1 - Rlbp
        RMminuss = 1 - Rlbs
        rMplusp = 1 + rlbp
        rMpluss = 1 + rlbs
        rMminusp = 1 - rlbp
        rMminuss = 1 - rlbs 
    ####
    
    # r factors for different input and output polarizations
    rMRpP = - (RMminusp * rMminusp**2 * wl1w**2 * wl2w * math.cos(PHI)**3 * XXX) \
            - (2 * RMminusp * rMminusp**2 * wl1w**2 * wl2w * math.sin(PHI) * math.cos(PHI)**2 * XXY) \
            - (2 * RMminusp * rMplusp * rMminusp * wl1w * wl2w * math.sin(THETA0) * math.cos(PHI)**2 * XXZ) \
            - (RMminusp * rMminusp**2 * wl1w**2 * wl2w * math.sin(PHI)**2 * math.cos(PHI) * XYY) \
            - (2 * RMminusp * rMplusp * rMminusp * wl1w * wl2w * math.sin(THETA0) * math.sin(PHI) * math.cos(PHI) * XYZ) \
            - (RMminusp * rMplusp**2 * wl2w * math.sin(THETA0)**2 * math.cos(PHI) * XZZ) \
            - (RMminusp * rMminusp**2 * wl1w**2 * wl2w * math.sin(PHI) * math.cos(PHI)**2 * YXX) \
            - (2 * RMminusp * rMminusp**2 * wl1w**2 * wl2w * math.sin(PHI)**2 * math.cos(PHI) * YYX) \
            - (2 * RMminusp * rMplusp * rMminusp * wl1w * wl2w * math.sin(THETA0) * math.sin(PHI) * math.cos(PHI) * YXZ) \
            - (RMminusp * rMminusp**2 * wl1w**2 * wl2w * math.sin(PHI)**3 * YYY) \
            - (2 * RMminusp * rMplusp * rMminusp * wl1w * wl2w * math.sin(THETA0) * math.sin(PHI)**2 * YYZ) \
            - (RMminusp * rMplusp**2 * wl2w * math.sin(THETA0)**2 * math.sin(PHI) * YZZ) \
            + (RMplusp * rMminusp**2 * wl1w**2 * math.sin(THETA0) * math.cos(PHI)**2 * ZXX) \
            + (2 * RMplusp * rMplusp * rMminusp * wl1w * math.sin(THETA0)**2 * math.cos(PHI) * ZXZ) \
            + (2 * RMplusp * rMminusp**2 * wl1w**2 * math.sin(THETA0) * math.sin(PHI) * math.cos(PHI) * ZXY) \
            + (RMplusp * rMminusp**2 * wl1w**2 * math.sin(THETA0) * math.sin(PHI)**2 * ZYY) \
            + (2 * RMplusp * rMplusp * rMminusp * wl1w * math.sin(THETA0)**2 * math.sin(PHI) * ZZY) \
            + (RMplusp * rMplusp**2 * math.sin(PHI)**3 * ZZZ)
    rMRpS = - (rMminusp**2 * wl1w**2 * math.sin(PHI) * math.cos(PHI)**2 * XXX) \
            - (2 * rMminusp**2 * wl1w**2 * math.sin(PHI)**2 * math.cos(PHI) * XXY) \
            - (2 * rMplusp * rMminusp * wl1w * math.sin(THETA0) * math.sin(PHI) * math.cos(PHI) * XXZ) \
            - (rMminusp**2 * wl1w**2 * math.sin(PHI)**3 * XYY) \
            - (2 * rMplusp * rMminusp * wl1w * math.sin(THETA0) * math.sin(PHI)**2 * XYZ) \
            - (rMplusp**2 * math.sin(THETA0)**2 * math.sin(PHI) * XZZ) \
            + (rMminusp**2 * wl1w**2 * math.cos(PHI)**3 * YXX) \
            + (2 * rMminusp**2 * wl1w**2 * math.sin(PHI) * math.cos(PHI)**2 * YYX) \
            + (2 * rMplusp * rMminusp * wl1w * math.sin(THETA0) * math.cos(PHI)**2 * YXZ) \
            + (rMminusp**2 * wl1w**2 * math.sin(PHI)**2 * math.cos(PHI) * YYY) \
            + (2 * rMplusp * rMminusp * wl1w * math.sin(THETA0) * math.sin(PHI) * math.cos(PHI) * YYZ) \
            + (rMplusp**2 * math.sin(THETA0)**2 * math.cos(PHI) * YZZ)
    rMRsP = - (RMminusp * wl2w * math.sin(PHI)**2 * math.cos(PHI) * XXX) \
            + (RMminusp * wl2w * 2 * math.sin(PHI) * math.cos(PHI)**2 * XXY) \
            - (RMminusp * wl2w * math.cos(PHI)**3 * XYY) \
            - (RMminusp * wl2w * math.sin(PHI)**3 * YXX) \
            + (RMminusp * wl2w * 2 * math.sin(PHI)**2 * math.cos(PHI) * YYX) \
            - (RMminusp * wl2w * math.sin(PHI) * math.cos(PHI)**2 * YYY) \
            + (RMplusp * math.sin(THETA0) * math.sin(PHI)**2 * ZXX) \
            - (RMplusp * math.sin(THETA0) * 2 * math.sin(PHI) * math.cos(PHI) * ZXY) \
            + (RMplusp * math.sin(THETA0) * math.cos(PHI)**2 * ZYY)
    rMRsS = - (math.sin(PHI)**3 * XXX) \
            + (2 * math.sin(PHI)**2 * math.cos(PHI) * XXY) \
            - (math.sin(PHI) * math.cos(PHI)**2 * XYY) \
            + (math.sin(PHI)**2 * math.cos(PHI) * YXX) \
            + (math.cos(PHI)**3 * YYY) \
            - (2 * math.sin(PHI) * math.cos(PHI)**2 * YYX)
    
    # fresnel factors multiplied out for ease of debugging
    GammaMRpP = (Tvlp/Nl) * (tvlp/nl)**2
    GammaMRpS = Tvls * RMpluss * (tvlp/nl)**2
    GammaMRsP = (Tvlp/Nl) * (tvls * rMpluss)**2
    GammaMRsS = Tvls * RMpluss * (tvls * rMpluss)**2
    
    sumGammaMRpP = sumGammaMRpP + GammaMRpP
    sumGammaMRpS = sumGammaMRpS + GammaMRpS
    sumGammaMRsP = sumGammaMRsP + GammaMRsP
    sumGammaMRsS = sumGammaMRsS + GammaMRsS
    sumrMRpP = sumrMRpP + rMRpP
    sumrMRpS = sumrMRpS + rMRpS
    sumrMRsP = sumrMRsP + rMRsP
    sumrMRsS = sumrMRsS + rMRsS

# R factors for different input and output polarizations (in cm^2/W)
RMRpP = shgyield(sumGammaMRpP, sumrMRpP)
RMRpS = shgyield(sumGammaMRpS, sumrMRpS)
RMRsP = shgyield(sumGammaMRsP, sumrMRsP)
RMRsS = shgyield(sumGammaMRsS, sumrMRsS)


#### Output ####
output = "layered.dat"
savefile(output, ONEE, RMRpP, RMRpS, RMRsP, RMRsS)
