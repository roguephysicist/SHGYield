#!/usr/bin/env python
"""
shgyield.py is a python program designed to calculate the nonlinear reflection
coefficient for semiconductor surfaces. It works in conjunction with the matrix
elements calculated with ABINIT, an open source ab initio software, and TINIBA,
our in-house optical calculation software. For a complete overview of the
theory, see PRB 94, 115314 (2016).

TODO:
* Make an epsilon spline so there is no need to halve energy range
* Dict comprehension to reduce array creation
* Consider switching over to matrix expressions for 'rpP', 'rpS', 'rsP', and 'rsS'

required packages:
`sys, yaml, numpy, scipy`

usage:
`python shgyield.py input.yml`
"""

import sys
import yaml
import numpy as np
from scipy import constants, ndimage


def parse_input(infile):
    '''
    Parses the YAML input file specified on the command line. Returns dictionary
    with variables and values from the input file.
    '''
    with open(infile) as data_file:
        params = yaml.load(data_file)
    return params

def broad(target, sigma):
    '''
    A function for applying Gaussian broadening on the final output data.
    '''
    data = ndimage.filters.gaussian_filter(target, sigma)
    return data

def savefile(file, freq, val1, val2, val3, val4):
    '''
    Saves specified dataset to file with the following columns:
    Energy(1w)    R_{pP}    R_{pS}    R_{sP}    R_{sS}
    '''
    data = np.column_stack((freq, val1, val2, val3, val4))
    np.savetxt(file, data, fmt=('%05.2f', '%.8e', '%.8e', '%.8e', '%.8e'),
               delimiter='    ',
               header='RiF 1e-20 (cm^2/W)\n1w(eV) RpP'+15*" "+\
                      'RpS'+15*" "+'RsP'+15*" "+'RsS')

def epsload(in_file, norm):
    '''
    Reads calculated chi1 file, that is organized as
    Energy(1w) Re[chi_xx] Im[chi_xx] Re[chi_yy] Im[chi_yy] Re[chi_zz] Im[chi_zz].
    Converts to epsilon = 1 + 4*pi*chi1, and normalized chi1 according to
    PRB 92, 245308 (2015). Returns the averaged epsilon
    epsilon_{avg} = (epsilon^{xx} + epsilon^{yy} + epsilon^{zz})/3
    in a numpy array.
    '''
    data = np.loadtxt(in_file, unpack=True, skiprows=1)
    real = (data[1] + data[3] + data[5])/3      # real average
    imag = (data[2] + data[4] + data[6])/3      # imag average
    coma = real + 1j * imag                     # complex average
    epsa = 1 + (4 * np.pi * coma * norm) # epsilon with normalization
    return epsa

def shgload(infile):
    '''
    Loads chi^{abc} listed in the input file, with the following columns:
    Energy(1w) Re[1w] Im[1w] Re[2w] Im[2w]. Sums 1w and 2w real and imaginary
    parts, and handles scale and units appropriately. Returns numpy array.
    '''
    data = np.loadtxt(infile, unpack=True, skiprows=1)
    comp = (data[1] + data[3]) + 1j * (data[2] + data[4])
    chi2 = TINIBASCALE * PM2TOM2 * comp[:MAXE]
    return chi2

def freflc(pol, veci, vecj, epsi, epsj):
    '''
    Generic reflection fresnel factors, see Eq. (13) of PRB 94, 115314 (2016).
    '''
    if pol == "p":
        factor = ((veci * epsj) - (vecj * epsi))/((veci * epsj) + (vecj * epsi))
    elif pol == "s":
        factor = (veci - vecj)/(veci + vecj)
    return factor

def ftrans(pol, veci, vecj, epsi, epsj):
    '''
    Generic transmission fresnel factors, see Eq. (13) of PRB 94, 115314 (2016).
    '''
    if pol == "p":
        factor = (2 * veci * np.sqrt(epsi * epsj))/(veci * epsj + vecj * epsi)
    elif pol == "s":
        factor = (2 * veci)/(veci + vecj)
    return factor

def shgyield(gamma, riF):
    '''
    Calculates the final broadened SHG yield, ready to be written to file.
    See Eq. (38) of PRB 94, 115314 (2016).
    '''
    RiF = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * np.absolute((1/nl) * gamma * riF)**2
    broadened = broad(RiF, SIGMA)
    return broadened


## Initialization: Parse input file, establish relevant modes.
PARAM = parse_input(sys.argv[1])  # Parses input file
MODE = PARAM['mode']              # Establishes the 'layer model' to be used;
                                  # see PRB 93, 235304 (2016).


## Energy: establishes energy range
## Assumes chi1/2 components are taken over 0-20 eV with 2001 (0.01 eV) steps
MAXE = 1000 # Upper energy limit, usually 1000 = 10 eV
ONEE = np.linspace(0.01, float(MAXE)/100, MAXE) # 1w energy array, 0.01-10 eV


## Constants and parameters
HBAR = constants.value("Planck constant over 2 pi in eV s")
PLANCK = constants.value("Planck constant in eV s")
EPS0 = constants.epsilon_0                   # \epsilon_{0} in F/m
LSPEED = constants.c                         # Speed of light in m/s
PM2TOM2 = 1e-24                              # Convert from pm^2 to m^2
M2TOCM2 = 1e4                                # Convert from m^2 to cm^2
TINIBASCALE = 1e6                            # Scaling chi2 in 1e6 (pm^2/V)
SCALE = 1e20                                 # Final yield in 1e-20 (cm^2/W)
THETA0 = np.radians(PARAM['parameters']['theta']) # Converts theta to radians
PHI = np.radians(PARAM['parameters']['phi']) # Converts phi to radians
SIGMA = PARAM['parameters']['sigma']         # Std. dev. for gaussian broadening
PREFACTOR = 1/(2 * EPS0 * HBAR**2 * LSPEED**3 * np.cos(THETA0)**2)


## Linear responses: chi1 and epsilons
CHI1NORM = PARAM['chi1']['norm']                # Normalization for layered chi1
epsb = epsload(PARAM['chi1']['chib'], 1)        # Epsilon from chi1, bulk
epsv1w = 1                              # Epsilon for vacuum = 1
epsv2w = 1                              # Epsilon for vacuum = 1
epsb1w = epsb[:MAXE]                    # Epsilon for bulk, 1w
epsb2w = epsb[1::2][:MAXE]              # Epsilon for bulk, 2w


## Reflection model, see PRB 93, 235304 (2016).
if MODE == "3-layer": # The incident fields and SHG both occur in the thin layer (l)
    epsl = epsload(PARAM['chi1']['chil'], CHI1NORM) # Epsilon from chi1, layered, normalized
    epsl1w = epsl[:MAXE] # Epsilon for layer, 1w
    epsl2w = epsl[1::2][:MAXE] # Epsilon for layer, 2w
elif MODE == "2-layer-fresnel": # The incident fields in bulk, SHG in vacuum
    epsl1w = epsb1w
    epsl2w = epsv2w
elif MODE == "2-layer-bulk": # Both incident fields and SHG in bulk
    epsl1w = epsb1w
    epsl2w = epsb2w
elif MODE == "2-layer-vacuum": # Both incident fields and SHG in vacuum
    epsl1w = epsv2w
    epsl2w = epsv2w
elif MODE == "3-layer-hybrid": # Incident field in bulk, SHG in thin layer (l)
    epsl = epsload(PARAM['chi1']['chil'], CHI1NORM) # Epsilon from chi1, layered, normalized
    epsl1w = epsb1w
    epsl2w = epsl[1::2][:MAXE] # Epsilon for layer, 2w


## Nonlinear responses: chi2 components to be used in the formulas below.
## There are 18 possible components for SHG. We have 4 scenarios regarding the
## input file and how it is interpreted:
## 1. Component not listed (or commented)                  (chi2^{abc} = 0)
## 2. Component specified as 0                             (chi2^{abc} = 0)
## 3. Component specified with path to file                (chi2^{abc} loaded from file)
## 4. Component has symmetry relation, such as 'yyx: -xxx' (chi2^{yyx} = -chi2^{xxx})
##
## The end result is a dictionary (CHI2) with the component name as the 'key',
## and a numpy array with the component data as the 'value'.
CHI2 = {}
CHI2_equiv = {}
all_components = ['xxx', 'xyy', 'xzz', 'xyz', 'xxz', 'xxy',
                  'yxx', 'yyy', 'yzz', 'yyz', 'yxz', 'yxy',
                  'zxx', 'zyy', 'zzz', 'zyz', 'zxz', 'zxy']
for component in all_components:
    if component in PARAM['chi2']:
        value = PARAM['chi2'][component]
        try:
            shg = shgload(value)
        except (ValueError, OSError, IOError):
            if value == 0:
                shg = 0
            else:
                CHI2_equiv[component] = value
                continue
    elif component not in PARAM['chi2']:
        shg = 0
    CHI2[component] = shg
for key, value in CHI2_equiv.items():
    equivalence = value.split('-')
    if equivalence[0]:
        CHI2[key] = CHI2[equivalence[0]]
    elif not equivalence[0]:
        CHI2[key] = -1 * CHI2[equivalence[1]]


## Indices of refraction for 1w (n) and 2w (N), where n = sqrt{epsilon}
nv = np.sqrt(epsv1w) # Index of refraction for vacuum, 1w
Nv = np.sqrt(epsv2w) # Index of refraction for vacuum, 2w
nb = np.sqrt(epsb1w) # Index of refraction for bulk, 1w
Nb = np.sqrt(epsb2w) # Index of refraction for bulk, 2w
nl = np.sqrt(epsl1w) # Index of refraction for layer, 1w
Nl = np.sqrt(epsl2w) # Index of refraction for layer, 2w


## Wave vectors for 1w and 2w, see Eq. (6) of PRB 94, 115314 (2016).
wv1w = np.sqrt(epsv1w - (np.sin(THETA0)**2)) # Wave vector for vacuum, 1w
wv2w = np.sqrt(epsv2w - (np.sin(THETA0)**2)) # Wave vector for vacuum, 2w
wb1w = np.sqrt(epsb1w - (np.sin(THETA0)**2)) # Wave vector for bulk, 1w
wb2w = np.sqrt(epsb2w - (np.sin(THETA0)**2)) # Wave vector for bulk, 2w
wl1w = np.sqrt(epsl1w - (np.sin(THETA0)**2)) # Wave vector for layer, 1w
wl2w = np.sqrt(epsl2w - (np.sin(THETA0)**2)) # Wave vector for layer, 2w


## Fresnel factors for 1w and 2w, s and p polarizations. See Eqs. (13) and (14)
## of PRB 94, 115314 (2016).
tvls = ftrans("s", wv1w, wl1w, epsv1w, epsl1w) # Transmission, 1w, vacuum-layer, s
tvlp = ftrans("p", wv1w, wl1w, epsv1w, epsl1w) # Transmission, 1w, vacuum-layer, p
tlbs = ftrans("s", wl1w, wb1w, epsl1w, epsb1w) # Transmission, 1w, layer-bulk,   s
tlbp = ftrans("p", wl1w, wb1w, epsl1w, epsb1w) # Transmission, 1w, layer-bulk,   p
Tvls = ftrans("s", wv2w, wl2w, epsv2w, epsl2w) # Transmission, 2w, vacuum-layer, s
Tvlp = ftrans("p", wv2w, wl2w, epsv2w, epsl2w) # Transmission, 2w, vacuum-layer, p
Tlbs = ftrans("s", wl2w, wb2w, epsl2w, epsb2w) # Transmission, 2w, layer-bulk,   s
Tlbp = ftrans("p", wl2w, wb2w, epsl2w, epsb2w) # Transmission, 2w, layer-bulk,   p
rvls = freflc("s", wv1w, wl1w, epsv1w, epsl1w) # Reflection,   1w, vacuum-layer, s
rvlp = freflc("p", wv1w, wl1w, epsv1w, epsl1w) # Reflection,   1w, vacuum-layer, p
rlbs = freflc("s", wl1w, wb1w, epsl1w, epsb1w) # Reflection,   1w, layer-bulk,   s
rlbp = freflc("p", wl1w, wb1w, epsl1w, epsb1w) # Reflection,   1w, layer-bulk,   p
Rvls = freflc("s", wv2w, wl2w, epsv2w, epsl2w) # Reflection,   2w, vacuum-layer, s
Rvlp = freflc("p", wv2w, wl2w, epsv2w, epsl2w) # Reflection,   2w, vacuum-layer, p
Rlbs = freflc("s", wl2w, wb2w, epsl2w, epsb2w) # Reflection,   2w, layer-bulk,   s
Rlbp = freflc("p", wl2w, wb2w, epsl2w, epsb2w) # Reflection,   2w, layer-bulk,   p


## Multiple reflections framework. See Eqs. (16), (17), (21), (22), (26),
## and (30) of PRB 94, 115314 (2016).
if PARAM['multiref']['enable'] and MODE == '3-layer':
    thickness = PARAM['multiref']['thickness'] # Thickness d of the thin layer \ell
    depth = PARAM['multiref']['depth']         # Depth d2 of the polarization sheet
    varphi = 4 * np.pi * ((ONEE * thickness * 1e-9)/(PLANCK * LSPEED)) * wl1w
    delta = 8 * np.pi * ((ONEE * thickness * 1e-9)/(PLANCK * LSPEED)) * wl2w
    if depth == "average":
        RMpav = (Rlbp * np.exp(1j * delta/2))/\
                (1 + (Rvlp * Rlbp * np.exp(1j * delta))) * np.sin(delta/2)/(delta/2)
        RMsav = (Rlbs * np.exp(1j * delta/2))/\
                (1 + (Rvls * Rlbs * np.exp(1j * delta))) * np.sin(delta/2)/(delta/2)
        RMplusp = 1 + RMpav
        RMpluss = 1 + RMsav
        RMminusp = 1 - RMpav
        RMminuss = 1 - RMsav
    else:
        D2 = float(depth)
        delta0 = 8 * np.pi * ((ONEE * D2 * 1e-9)/(PLANCK * LSPEED)) * wl2w
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
elif not PARAM['multiref']['enable'] or MODE != '3-layer':
    RMplusp = 1 + Rlbp
    RMpluss = 1 + Rlbs
    RMminusp = 1 - Rlbp
    RMminuss = 1 - Rlbs
    rMplusp = 1 + rlbp
    rMpluss = 1 + rlbs
    rMminusp = 1 - rlbp
    rMminuss = 1 - rlbs


## r factors for different input and output polarizations. See Eqs. (50), (55),
## (60), and (65) of PRB 94, 115314 (2016).
### r_{pP}
rpP = - (RMminusp * rMminusp**2 * wl1w**2 * wl2w * np.cos(PHI)**3 * CHI2['xxx']) \
      - (2 * RMminusp * rMminusp**2 * wl1w**2 * wl2w * np.sin(PHI) * np.cos(PHI)**2 * CHI2['xxy']) \
      - (2 * RMminusp * rMplusp * rMminusp * wl1w * wl2w * np.sin(THETA0) * np.cos(PHI)**2 * CHI2['xxz']) \
      - (RMminusp * rMminusp**2 * wl1w**2 * wl2w * np.sin(PHI)**2 * np.cos(PHI) * CHI2['xyy']) \
      - (2 * RMminusp * rMplusp * rMminusp * wl1w * wl2w * np.sin(THETA0) * np.sin(PHI) * np.cos(PHI) * CHI2['xyz']) \
      - (RMminusp * rMplusp**2 * wl2w * np.sin(THETA0)**2 * np.cos(PHI) * CHI2['xzz']) \
      - (RMminusp * rMminusp**2 * wl1w**2 * wl2w * np.sin(PHI) * np.cos(PHI)**2 * CHI2['yxx']) \
      - (2 * RMminusp * rMminusp**2 * wl1w**2 * wl2w * np.sin(PHI)**2 * np.cos(PHI) * CHI2['yxy']) \
      - (2 * RMminusp * rMplusp * rMminusp * wl1w * wl2w * np.sin(THETA0) * np.sin(PHI) * np.cos(PHI) * CHI2['yxz']) \
      - (RMminusp * rMminusp**2 * wl1w**2 * wl2w * np.sin(PHI)**3 * CHI2['yyy']) \
      - (2 * RMminusp * rMplusp * rMminusp * wl1w * wl2w * np.sin(THETA0) * np.sin(PHI)**2 * CHI2['yyz']) \
      - (RMminusp * rMplusp**2 * wl2w * np.sin(THETA0)**2 * np.sin(PHI) * CHI2['yzz']) \
      + (RMplusp * rMminusp**2 * wl1w**2 * np.sin(THETA0) * np.cos(PHI)**2 * CHI2['zxx']) \
      + (2 * RMplusp * rMplusp * rMminusp * wl1w * np.sin(THETA0)**2 * np.cos(PHI) * CHI2['zxz']) \
      + (2 * RMplusp * rMminusp**2 * wl1w**2 * np.sin(THETA0) * np.sin(PHI) * np.cos(PHI) * CHI2['zxy']) \
      + (RMplusp * rMminusp**2 * wl1w**2 * np.sin(THETA0) * np.sin(PHI)**2 * CHI2['zyy']) \
      + (2 * RMplusp * rMplusp * rMminusp * wl1w * np.sin(THETA0)**2 * np.sin(PHI) * CHI2['zyz']) \
      + (RMplusp * rMplusp**2 * np.sin(PHI)**3 * CHI2['zzz'])
### r_{pS}
rpS = - (rMminusp**2 * wl1w**2 * np.sin(PHI) * np.cos(PHI)**2 * CHI2['xxx']) \
      - (2 * rMminusp**2 * wl1w**2 * np.sin(PHI)**2 * np.cos(PHI) * CHI2['xxy']) \
      - (2 * rMplusp * rMminusp * wl1w * np.sin(THETA0) * np.sin(PHI) * np.cos(PHI) * CHI2['xxz']) \
      - (rMminusp**2 * wl1w**2 * np.sin(PHI)**3 * CHI2['xyy']) \
      - (2 * rMplusp * rMminusp * wl1w * np.sin(THETA0) * np.sin(PHI)**2 * CHI2['xyz']) \
      - (rMplusp**2 * np.sin(THETA0)**2 * np.sin(PHI) * CHI2['xzz']) \
      + (rMminusp**2 * wl1w**2 * np.cos(PHI)**3 * CHI2['yxx']) \
      + (2 * rMminusp**2 * wl1w**2 * np.sin(PHI) * np.cos(PHI)**2 * CHI2['yxy']) \
      + (2 * rMplusp * rMminusp * wl1w * np.sin(THETA0) * np.cos(PHI)**2 * CHI2['yxz']) \
      + (rMminusp**2 * wl1w**2 * np.sin(PHI)**2 * np.cos(PHI) * CHI2['yyy']) \
      + (2 * rMplusp * rMminusp * wl1w * np.sin(THETA0) * np.sin(PHI) * np.cos(PHI) * CHI2['yyz']) \
      + (rMplusp**2 * np.sin(THETA0)**2 * np.cos(PHI) * CHI2['yzz'])
### r_{sP}
rsP = - (RMminusp * wl2w * np.sin(PHI)**2 * np.cos(PHI) * CHI2['xxx']) \
      + (RMminusp * wl2w * 2 * np.sin(PHI) * np.cos(PHI)**2 * CHI2['xxy']) \
      - (RMminusp * wl2w * np.cos(PHI)**3 * CHI2['xyy']) \
      - (RMminusp * wl2w * np.sin(PHI)**3 * CHI2['yxx']) \
      + (RMminusp * wl2w * 2 * np.sin(PHI)**2 * np.cos(PHI) * CHI2['yxy']) \
      - (RMminusp * wl2w * np.sin(PHI) * np.cos(PHI)**2 * CHI2['yyy']) \
      + (RMplusp * np.sin(THETA0) * np.sin(PHI)**2 * CHI2['zxx']) \
      - (RMplusp * np.sin(THETA0) * 2 * np.sin(PHI) * np.cos(PHI) * CHI2['zxy']) \
      + (RMplusp * np.sin(THETA0) * np.cos(PHI)**2 * CHI2['zyy'])
### r_{sS}
rsS = - (np.sin(PHI)**3 * CHI2['xxx']) \
      + (2 * np.sin(PHI)**2 * np.cos(PHI) * CHI2['xxy']) \
      - (np.sin(PHI) * np.cos(PHI)**2 * CHI2['xyy']) \
      + (np.sin(PHI)**2 * np.cos(PHI) * CHI2['yxx']) \
      + (np.cos(PHI)**3 * CHI2['yyy']) \
      - (2 * np.sin(PHI) * np.cos(PHI)**2 * CHI2['yxy'])


## Gamma prefactor for different polarizations. See Eqs. (49), (54), (59), (64)
## of PRB 94, 115314 (2016).
GammapP = (Tvlp/Nl) * (tvlp/nl)**2             # p-in, P-out
GammapS = Tvls * RMpluss * (tvlp/nl)**2        # p-in, S-out
GammasP = (Tvlp/Nl) * (tvls * rMpluss)**2      # s-in, P-out
GammasS = Tvls * RMpluss * (tvls * rMpluss)**2 # s-in, S-out


## Final SHG yield for different input and output polarizations (in cm^2/W).
## See Eqs. (44) and (38) of PRB 94, 115314 (2016).
RpP = shgyield(GammapP, rpP) # p-in, P-out
RpS = shgyield(GammapS, rpS) # p-in, S-out
RsP = shgyield(GammasP, rsP) # s-in, P-out
RsS = shgyield(GammasS, rsS) # s-in, S-out


## Write output to file.
savefile(PARAM['output'], ONEE, RpP, RpS, RsP, RsS)
