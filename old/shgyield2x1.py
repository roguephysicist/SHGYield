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

# input
MODE = "3layer"
theta = 65
phi = 00
chib = "data/bulk/chi1-vnl.sm_0.075_xx_yy_zz_3107_25-nospin_scissor_0.50_Nc_26"
chil = "data/2x1/calChi1-vnl.sm_0.075_xx_yy_zz_244_half-slab_10-nospin_scissor_0.50_Nc_130"
xxx = "data/2x1/shgC.vnl.sm_0.075_xxx_244_half-slab_10-nospin_scissor_0.50_Nc_130"
xxy = "data/2x1/shgC.vnl.sm_0.075_xxy_244_half-slab_10-nospin_scissor_0.50_Nc_130"
xxz = "data/2x1/shgC.vnl.sm_0.075_xxz_244_half-slab_10-nospin_scissor_0.50_Nc_130"
xyx = "data/2x1/shgC.vnl.sm_0.075_xyx_244_half-slab_10-nospin_scissor_0.50_Nc_130"
xyy = "data/2x1/shgC.vnl.sm_0.075_xyy_244_half-slab_10-nospin_scissor_0.50_Nc_130"
xyz = "data/2x1/shgC.vnl.sm_0.075_xyz_244_half-slab_10-nospin_scissor_0.50_Nc_130"
xzx = "data/2x1/shgC.vnl.sm_0.075_xzx_244_half-slab_10-nospin_scissor_0.50_Nc_130"
xzy = "data/2x1/shgC.vnl.sm_0.075_xzy_244_half-slab_10-nospin_scissor_0.50_Nc_130"
xzz = "data/2x1/shgC.vnl.sm_0.075_xzz_244_half-slab_10-nospin_scissor_0.50_Nc_130"
yxx = "data/2x1/shgC.vnl.sm_0.075_yxx_244_half-slab_10-nospin_scissor_0.50_Nc_130"
yxy = "data/2x1/shgC.vnl.sm_0.075_yxy_244_half-slab_10-nospin_scissor_0.50_Nc_130"
yxz = "data/2x1/shgC.vnl.sm_0.075_yxz_244_half-slab_10-nospin_scissor_0.50_Nc_130"
yyx = "data/2x1/shgC.vnl.sm_0.075_yyx_244_half-slab_10-nospin_scissor_0.50_Nc_130"
yyy = "data/2x1/shgC.vnl.sm_0.075_yyy_244_half-slab_10-nospin_scissor_0.50_Nc_130"
yyz = "data/2x1/shgC.vnl.sm_0.075_yyz_244_half-slab_10-nospin_scissor_0.50_Nc_130"
yzx = "data/2x1/shgC.vnl.sm_0.075_yzx_244_half-slab_10-nospin_scissor_0.50_Nc_130"
yzy = "data/2x1/shgC.vnl.sm_0.075_yzy_244_half-slab_10-nospin_scissor_0.50_Nc_130"
yzz = "data/2x1/shgC.vnl.sm_0.075_yzz_244_half-slab_10-nospin_scissor_0.50_Nc_130"
zxx = "data/2x1/shgC.vnl.sm_0.075_zxx_244_half-slab_10-nospin_scissor_0.50_Nc_130"
zxy = "data/2x1/shgC.vnl.sm_0.075_zxy_244_half-slab_10-nospin_scissor_0.50_Nc_130"
zxz = "data/2x1/shgC.vnl.sm_0.075_zxz_244_half-slab_10-nospin_scissor_0.50_Nc_130"
zyx = "data/2x1/shgC.vnl.sm_0.075_zyx_244_half-slab_10-nospin_scissor_0.50_Nc_130"
zyy = "data/2x1/shgC.vnl.sm_0.075_zyy_244_half-slab_10-nospin_scissor_0.50_Nc_130"
zyz = "data/2x1/shgC.vnl.sm_0.075_zyz_244_half-slab_10-nospin_scissor_0.50_Nc_130"
zzx = "data/2x1/shgC.vnl.sm_0.075_zzx_244_half-slab_10-nospin_scissor_0.50_Nc_130"
zzy = "data/2x1/shgC.vnl.sm_0.075_zzy_244_half-slab_10-nospin_scissor_0.50_Nc_130"
zzz = "data/2x1/shgC.vnl.sm_0.075_zzz_244_half-slab_10-nospin_scissor_0.50_Nc_130"
output = "R2x1." + str(phi) + "." + MODE + ".dat"


## some variables
CHI1NORM = 1.2535831585 # Normalization yo
TESTSIG = 3.19999995 # broadening (sigma = 0.075 eV)


#### Functions ####
def broad(target, sigma): 
    data = ndimage.filters.gaussian_filter(target, sigma)
    return data

def savefile(file, freq, val1, val2): # saves to file with header
    data = np.column_stack((freq, val1, val2))
    np.savetxt(file, data,
            fmt=('%05.2f', '%.14e', '%.14e'),
            delimiter='    ')

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
    broadened = broad(RiF, TESTSIG)
    return broadened


#### Energy ####
## assumes range from 0 to 20 with 2001 steps
MAXE = 1000
ONEE = np.linspace(0.01, float(MAXE)/100, MAXE) # 1w energy array

#### Constants ####
HBAR = constants.value("Planck constant over 2 pi in eV s")
PLANCK = constants.value("Planck constant in eV s")
EPS0 = constants.epsilon_0 # (F/m)
LSPEED = constants.c # (m/s)
PM2TOM2 = 1e-24 # pm^2 to m^2
M2TOCM2 = 1e4 # m^2 to cm^2
TINIBASCALE = 1e6 # for scaling chi2 in 1e6 (pm^2/V)
SCALE = 1e20 # for R in 1e-20 (cm^2/W)
THETA0 = math.radians(float(theta)) # converts theta to radians
PHI = math.radians(float(phi)) # converts phi to radians
LAMBDA0 = (PLANCK * LSPEED * 1e9)/ONEE # In nanometers
PREFACTOR = 1 / (2 * EPS0 * HBAR**2 * LSPEED**3 * math.cos(THETA0)**2)


#### Math ####
# loads chi2, converts to m^2/V
XXX = (TINIBASCALE * PM2TOM2 * shgcomp(xxx))
XXY = (TINIBASCALE * PM2TOM2 * shgcomp(xxy))
XXZ = (TINIBASCALE * PM2TOM2 * shgcomp(xxz))
XYX = (TINIBASCALE * PM2TOM2 * shgcomp(xyx))
XYY = (TINIBASCALE * PM2TOM2 * shgcomp(xyy))
XYZ = (TINIBASCALE * PM2TOM2 * shgcomp(xyz))
XZX = (TINIBASCALE * PM2TOM2 * shgcomp(xzx))
XZY = (TINIBASCALE * PM2TOM2 * shgcomp(xzy))
XZZ = (TINIBASCALE * PM2TOM2 * shgcomp(xzz))
YXX = (TINIBASCALE * PM2TOM2 * shgcomp(yxx))
YXY = (TINIBASCALE * PM2TOM2 * shgcomp(yxy))
YXZ = (TINIBASCALE * PM2TOM2 * shgcomp(yxz))
YYX = (TINIBASCALE * PM2TOM2 * shgcomp(yyx))
YYY = (TINIBASCALE * PM2TOM2 * shgcomp(yyy))
YYZ = (TINIBASCALE * PM2TOM2 * shgcomp(yyz))
YZX = (TINIBASCALE * PM2TOM2 * shgcomp(yzx))
YZY = (TINIBASCALE * PM2TOM2 * shgcomp(yzy))
YZZ = (TINIBASCALE * PM2TOM2 * shgcomp(yzz))
ZXX = (TINIBASCALE * PM2TOM2 * shgcomp(zxx))
ZXY = (TINIBASCALE * PM2TOM2 * shgcomp(zxy))
ZXZ = (TINIBASCALE * PM2TOM2 * shgcomp(zxz))
ZYX = (TINIBASCALE * PM2TOM2 * shgcomp(zyx))
ZYY = (TINIBASCALE * PM2TOM2 * shgcomp(zyy))
ZYZ = (TINIBASCALE * PM2TOM2 * shgcomp(zyz))
ZZX = (TINIBASCALE * PM2TOM2 * shgcomp(zzx))
ZZY = (TINIBASCALE * PM2TOM2 * shgcomp(zzy))
ZZZ = (TINIBASCALE * PM2TOM2 * shgcomp(zzz))

# creates epsilons from chi1 responses
epsl = epsilon(chil, CHI1NORM)
epsb = epsilon(chib, 1)
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

# shg radiation
rsP = Nl**2 * wb2w * (
        - math.sin(PHI)**2 * math.cos(PHI) * XXX
        + 2 * math.sin(PHI) * math.cos(PHI)**2 * XXY
        - math.cos(PHI)**3 * XYY) \
    + Nl**2 * wb2w * (
        - math.sin(PHI)**3 * YXX
        + 2 * math.sin(PHI)**2 * math.cos(PHI) * YXY
        - math.sin(PHI) * math.cos(PHI)**2 * YYY) \
    + Nb**2 * math.sin(THETA0) * (
          math.sin(PHI)**2 * ZXX
        - 2 * math.sin(PHI) * math.cos(PHI) * ZXY
        + math.cos(PHI)**2 * ZYY)
GammasP = ((Tvlp * Tlbp)/(Nl**2 * Nb)) * (tvls * tlbs)**2
RsP = shgyield(GammasP, rsP)

rsS = (- math.sin(PHI)**3 * XXX) + (2 * math.sin(PHI)**2 * math.cos(PHI) * XXY) - (math.sin(PHI) * math.cos(PHI)**2 * XYY) + (math.sin(PHI)**2 * math.cos(PHI) * YXX) - (2 * math.sin(PHI) * math.cos(PHI)**2 * YXY) + (math.cos(PHI)**3 * YYY)
GammasS = Tvls * Tlbs * (tvls * tlbs)**2
RsS = shgyield(GammasS, rsS)


#### Output ####
savefile(output, ONEE, RsP, RsS)
