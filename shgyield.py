"""
shgyield.py is a python program designed to calculate the nonlinear reflection
coefficient for semiconductor surfaces. It works in conjunction with the matrix
elements calculated with ABINIT, an open source ab initio software, and TINIBA,
our in-house optical calculation software.

For a complete overview of the theory, see PRB 94, 115314 (2016).

Tested with Anaconda Python 4.0.0.

requirements:
sys, math, numpy, scipy

usage:
python shgyield.py <sample.in>
"""

import sys
import math
import numpy as np
from scipy import constants, ndimage


def parse_input(infile):
    '''
    Parses the input file specified on the command line. Returns dictionary
    with variables and values from the input file.
    '''
    params = {}
    targetfile = open(infile)
    data = targetfile.readlines()
    for line in data:
        if line.strip().startswith('#') or line == '\n':
            pass
        else:
            key, value = line.partition('#')[0].split(':')
            params[key.strip()] = value.strip()
    targetfile.close()
    return params

def broad(target, sigma):
    '''
    A function for applying Gaussian broadening on the final output data.
    '''
    value = sigma * 42.666666666
    data = ndimage.filters.gaussian_filter(target, value)
    return data

def savefile(file, freq, val1, val2, val3, val4):
    '''
    Saves specified dataset to file with the following columns:
    Energy(1w)    R_{pP}    R_{pS}    R_{sP}    R_{sS}
    '''
    data = np.column_stack((freq, val1, val2, val3, val4))
    np.savetxt(file, data, fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e'),
               delimiter='    ',
               header='RiF 1e-20 (cm^2/W)\n1w(eV) RpP'+21*" "+\
                      'RpS'+21*" "+'RsP'+21*" "+'RsS')

def epsilon(in_file, norm):
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
    epsa = 1 + (4 * constants.pi * coma * norm) # epsilon with normalization
    return epsa

def shgcomp():
    '''
    Reads each chi^{abc} listed in the input file, with the following columns:
    Energy(1w) Re[1w] Im[1w] Re[2w] Im[2w].
    Checks for all 18 components possible for SHG; those not listed in the input
    file are assumed to be 0. Sums 1w and 2w real and imaginary parts, and
    returns a dictionary with the component name as the 'key', and a numpy array
    with the component data as the 'value'.
    '''
    chi2 = {}
    components = ['xxx', 'xxy', 'xxz', 'xyy', 'xyz', 'xzz',
                  'yxx', 'yyx', 'yxz', 'yyy', 'yyz', 'yzz',
                  'zxx', 'zxy', 'zxz', 'zyy', 'zzy', 'zzz']
    for response in components:
        if response in PARAM:
            try:
                factor = float(PARAM[response])
                shg = TINIBASCALE * PM2TOM2 * factor
            except ValueError:
                data = np.loadtxt(PARAM[response], unpack=True, skiprows=1)
                comp = (data[1] + data[3]) + 1j * (data[2] + data[4])
                shg = TINIBASCALE * PM2TOM2 * comp[:MAXE]
        elif response not in PARAM:
            shg = 0
        chi2[response] = shg
    return chi2

def fresnel(kind, i, j, pol, freq):
    '''
    Generic fresnel factors, see Eq. (13) of PRB 94, 115314 (2016).
    '''
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
            factor = (wi - wj)/(wi + wj)
    return factor

def shgyield(gamma, riF):
    '''
    Calculates the final broadened SHG yield, ready to be written to file.
    See Eq. (38) of PRB 94, 115314 (2016).
    '''
    RiF = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
                    np.absolute((1/nl) * gamma * riF)**2
    broadened = broad(RiF, SIGMA)
    return broadened


## Initialization: Parse input file, establish relevant modes.
PARAM = parse_input(sys.argv[1])  # Parses input file
MODE = str(PARAM['mode'])         # Establishes the 'layer model' to be used;
                                  # see PRB 93, 235304 (2016). Mostly deprecated.
MULTIREF = str(PARAM['multiref']) # Whether or not multiple reflections are
                                  # considered.

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
THETA0 = math.radians(float(PARAM['theta'])) # Converts theta to radians
PHI = math.radians(float(PARAM['phi']))      # Converts phi to radians
THICKNESS = float(PARAM['thickness'])        # Thickness d of the thin layer \ell
DEPTH = PARAM['depth']                       # Depth d2 of the polarization sheet
SIGMA = float(PARAM['sigma'])                # Std. dev. for gaussian broadening
LAMBDA0 = (PLANCK * LSPEED * 1e9)/ONEE       # Energy range expressed in nm
CHI1NORM = float(PARAM['norm'])              # DEBUG: normalization for chi1
PREFACTOR = 1/(2 * EPS0 * HBAR**2 * LSPEED**3 * math.cos(THETA0)**2)

## Components and pre-calculation
CHI2 = shgcomp()                        # Dictionary with all chi2 components
epsl = epsilon(PARAM['chil'], CHI1NORM) # Epsilon from chi1, layered, normalized
epsb = epsilon(PARAM['chib'], 1)        # Epsilon from chi1, bulk
epsv1w = 1                              # Epsilon for vacuum = 1
epsv2w = 1                              # Epsilon for vacuum = 1
epsb1w = epsb[:MAXE]                    # Epsilon for bulk, 1w
epsb2w = epsb[1::2][:MAXE]              # Epsilon for bulk, 2w
epsl1w = epsl[:MAXE]                    # Epsilon for layer, 1w
epsl2w = epsl[1::2][:MAXE]              # Epsilon for layer, 2w

## Reflection model, see PRB 93, 235304 (2016). This is mostly deprecated.
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

## Indices of refraction for 1w (n) and 2w (N), where n = sqrt{epsilon}
nv = np.sqrt(epsv1w) # Index of refraction for vacuum, 1w
Nv = np.sqrt(epsv2w) # Index of refraction for vacuum, 2w
nb = np.sqrt(epsb1w) # Index of refraction for bulk, 1w
Nb = np.sqrt(epsb2w) # Index of refraction for bulk, 2w
nl = np.sqrt(epsl1w) # Index of refraction for layer, 1w
Nl = np.sqrt(epsl2w) # Index of refraction for layer, 2w

## Wave vectors for 1w and 2w, see Eq. (6) of PRB 94, 115314 (2016).
wv1w = np.sqrt(epsv1w - (math.sin(THETA0)**2)) # Wave vector for vacuum, 1w
wv2w = np.sqrt(epsv2w - (math.sin(THETA0)**2)) # Wave vector for vacuum, 2w
wb1w = np.sqrt(epsb1w - (math.sin(THETA0)**2)) # Wave vector for bulk, 1w
wb2w = np.sqrt(epsb2w - (math.sin(THETA0)**2)) # Wave vector for bulk, 2w
wl1w = np.sqrt(epsl1w - (math.sin(THETA0)**2)) # Wave vector for layer, 1w
wl2w = np.sqrt(epsl2w - (math.sin(THETA0)**2)) # Wave vector for layer, 2w

## Fresnel factors for 1w and 2w, s and p polarizations. See Eqs. (13) and (14)
## of PRB 94, 115314 (2016).
tvls = fresnel("t", "v", ell1w, "s", "1w") # Transmission, 1w, vacuum-layer, s
tvlp = fresnel("t", "v", ell1w, "p", "1w") # Transmission, 1w, vacuum-layer, p
tlbs = fresnel("t", ell1w, "b", "s", "1w") # Transmission, 1w, layer-bulk,   s
tlbp = fresnel("t", ell1w, "b", "p", "1w") # Transmission, 1w, layer-bulk,   p
Tvls = fresnel("t", "v", ell2w, "s", "2w") # Transmission, 2w, vacuum-layer, s
Tvlp = fresnel("t", "v", ell2w, "p", "2w") # Transmission, 2w, vacuum-layer, p
Tlbs = fresnel("t", ell2w, "b", "s", "2w") # Transmission, 2w, layer-bulk,   s
Tlbp = fresnel("t", ell2w, "b", "p", "2w") # Transmission, 2w, layer-bulk,   p
rvls = fresnel("r", "v", ell1w, "s", "1w") # Reflection,   1w, vacuum-layer, s
rvlp = fresnel("r", "v", ell1w, "p", "1w") # Reflection,   1w, vacuum-layer, p
rlbs = fresnel("r", ell1w, "b", "s", "1w") # Reflection,   1w, layer-bulk,   s
rlbp = fresnel("r", ell1w, "b", "p", "1w") # Reflection,   1w, layer-bulk,   p
Rvls = fresnel("r", "v", ell2w, "s", "2w") # Reflection,   2w, vacuum-layer, s
Rvlp = fresnel("r", "v", ell2w, "p", "2w") # Reflection,   2w, vacuum-layer, p
Rlbs = fresnel("r", ell2w, "b", "s", "2w") # Reflection,   2w, layer-bulk,   s
Rlbp = fresnel("r", ell2w, "b", "p", "2w") # Reflection,   2w, layer-bulk,   p

## Multiple reflections framework. See Eqs. (16), (17), (21), (22), (26),
## and (30) of PRB 94, 115314 (2016).
if MULTIREF == "yes":
    varphi = 4 * math.pi * ((ONEE * THICKNESS * 1e-9)/(PLANCK * LSPEED)) * wl1w
    delta = 8 * math.pi * ((ONEE * THICKNESS * 1e-9)/(PLANCK * LSPEED)) * wl2w
    if DEPTH == "average":
        RMpav = (Rlbp * np.exp(1j * delta/2))/\
                (1 + (Rvlp * Rlbp * np.exp(1j * delta))) * np.sin(delta/2)/(delta/2)
        RMsav = (Rlbs * np.exp(1j * delta/2))/\
                (1 + (Rvls * Rlbs * np.exp(1j * delta))) * np.sin(delta/2)/(delta/2)
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

## r factors for different input and output polarizations. See Eqs. (50), (55),
## (60), and (65) of PRB 94, 115314 (2016).
### r_{pP}
rMRpP = - (RMminusp * rMminusp**2 * wl1w**2 * wl2w * math.cos(PHI)**3 * CHI2['xxx']) \
        - (2 * RMminusp * rMminusp**2 * wl1w**2 * wl2w * math.sin(PHI) * math.cos(PHI)**2 * CHI2['xxy']) \
        - (2 * RMminusp * rMplusp * rMminusp * wl1w * wl2w * math.sin(THETA0) * math.cos(PHI)**2 * CHI2['xxz']) \
        - (RMminusp * rMminusp**2 * wl1w**2 * wl2w * math.sin(PHI)**2 * math.cos(PHI) * CHI2['xyy']) \
        - (2 * RMminusp * rMplusp * rMminusp * wl1w * wl2w * math.sin(THETA0) * math.sin(PHI) * math.cos(PHI) * CHI2['xyz']) \
        - (RMminusp * rMplusp**2 * wl2w * math.sin(THETA0)**2 * math.cos(PHI) * CHI2['xzz']) \
        - (RMminusp * rMminusp**2 * wl1w**2 * wl2w * math.sin(PHI) * math.cos(PHI)**2 * CHI2['yxx']) \
        - (2 * RMminusp * rMminusp**2 * wl1w**2 * wl2w * math.sin(PHI)**2 * math.cos(PHI) * CHI2['yyx']) \
        - (2 * RMminusp * rMplusp * rMminusp * wl1w * wl2w * math.sin(THETA0) * math.sin(PHI) * math.cos(PHI) * CHI2['yxz']) \
        - (RMminusp * rMminusp**2 * wl1w**2 * wl2w * math.sin(PHI)**3 * CHI2['yyy']) \
        - (2 * RMminusp * rMplusp * rMminusp * wl1w * wl2w * math.sin(THETA0) * math.sin(PHI)**2 * CHI2['yyz']) \
        - (RMminusp * rMplusp**2 * wl2w * math.sin(THETA0)**2 * math.sin(PHI) * CHI2['yzz']) \
        + (RMplusp * rMminusp**2 * wl1w**2 * math.sin(THETA0) * math.cos(PHI)**2 * CHI2['zxx']) \
        + (2 * RMplusp * rMplusp * rMminusp * wl1w * math.sin(THETA0)**2 * math.cos(PHI) * CHI2['zxz']) \
        + (2 * RMplusp * rMminusp**2 * wl1w**2 * math.sin(THETA0) * math.sin(PHI) * math.cos(PHI) * CHI2['zxy']) \
        + (RMplusp * rMminusp**2 * wl1w**2 * math.sin(THETA0) * math.sin(PHI)**2 * CHI2['zyy']) \
        + (2 * RMplusp * rMplusp * rMminusp * wl1w * math.sin(THETA0)**2 * math.sin(PHI) * CHI2['zzy']) \
        + (RMplusp * rMplusp**2 * math.sin(PHI)**3 * CHI2['zzz'])
### r_{pS}
rMRpS = - (rMminusp**2 * wl1w**2 * math.sin(PHI) * math.cos(PHI)**2 * CHI2['xxx']) \
        - (2 * rMminusp**2 * wl1w**2 * math.sin(PHI)**2 * math.cos(PHI) * CHI2['xxy']) \
        - (2 * rMplusp * rMminusp * wl1w * math.sin(THETA0) * math.sin(PHI) * math.cos(PHI) * CHI2['xxz']) \
        - (rMminusp**2 * wl1w**2 * math.sin(PHI)**3 * CHI2['xyy']) \
        - (2 * rMplusp * rMminusp * wl1w * math.sin(THETA0) * math.sin(PHI)**2 * CHI2['xyz']) \
        - (rMplusp**2 * math.sin(THETA0)**2 * math.sin(PHI) * CHI2['xzz']) \
        + (rMminusp**2 * wl1w**2 * math.cos(PHI)**3 * CHI2['yxx']) \
        + (2 * rMminusp**2 * wl1w**2 * math.sin(PHI) * math.cos(PHI)**2 * CHI2['yyx']) \
        + (2 * rMplusp * rMminusp * wl1w * math.sin(THETA0) * math.cos(PHI)**2 * CHI2['yxz']) \
        + (rMminusp**2 * wl1w**2 * math.sin(PHI)**2 * math.cos(PHI) * CHI2['yyy']) \
        + (2 * rMplusp * rMminusp * wl1w * math.sin(THETA0) * math.sin(PHI) * math.cos(PHI) * CHI2['yyz']) \
        + (rMplusp**2 * math.sin(THETA0)**2 * math.cos(PHI) * CHI2['yzz'])
### r_{sP}
rMRsP = - (RMminusp * wl2w * math.sin(PHI)**2 * math.cos(PHI) * CHI2['xxx']) \
        + (RMminusp * wl2w * 2 * math.sin(PHI) * math.cos(PHI)**2 * CHI2['xxy']) \
        - (RMminusp * wl2w * math.cos(PHI)**3 * CHI2['xyy']) \
        - (RMminusp * wl2w * math.sin(PHI)**3 * CHI2['yxx']) \
        + (RMminusp * wl2w * 2 * math.sin(PHI)**2 * math.cos(PHI) * CHI2['yyx']) \
        - (RMminusp * wl2w * math.sin(PHI) * math.cos(PHI)**2 * CHI2['yyy']) \
        + (RMplusp * math.sin(THETA0) * math.sin(PHI)**2 * CHI2['zxx']) \
        - (RMplusp * math.sin(THETA0) * 2 * math.sin(PHI) * math.cos(PHI) * CHI2['zxy']) \
        + (RMplusp * math.sin(THETA0) * math.cos(PHI)**2 * CHI2['zyy'])
### r_{sS}
rMRsS = - (math.sin(PHI)**3 * CHI2['xxx']) \
        + (2 * math.sin(PHI)**2 * math.cos(PHI) * CHI2['xxy']) \
        - (math.sin(PHI) * math.cos(PHI)**2 * CHI2['xyy']) \
        + (math.sin(PHI)**2 * math.cos(PHI) * CHI2['yxx']) \
        + (math.cos(PHI)**3 * CHI2['yyy']) \
        - (2 * math.sin(PHI) * math.cos(PHI)**2 * CHI2['yyx'])

## Gamma prefactor for different polarizations. See Eqs. (49), (54), (59), (64)
## of PRB 94, 115314 (2016).
GammaMRpP = (Tvlp/Nl) * (tvlp/nl)**2             # p-in, P-out
GammaMRpS = Tvls * RMpluss * (tvlp/nl)**2        # p-in, S-out
GammaMRsP = (Tvlp/Nl) * (tvls * rMpluss)**2      # s-in, P-out
GammaMRsS = Tvls * RMpluss * (tvls * rMpluss)**2 # s-in, S-out

## Final SHG yield for different input and output polarizations (in cm^2/W).
## See Eqs. (44) and (38) of PRB 94, 115314 (2016).
RMRpP = shgyield(GammaMRpP, rMRpP) # p-in, P-out
RMRpS = shgyield(GammaMRpS, rMRpS) # p-in, S-out
RMRsP = shgyield(GammaMRsP, rMRsP) # s-in, P-out
RMRsS = shgyield(GammaMRsS, rMRsS) # s-in, S-out


## Write output to file.
savefile(PARAM['output'], ONEE, RMRpP, RMRpS, RMRsP, RMRsS)
