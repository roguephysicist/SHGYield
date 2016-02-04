#!/Users/sma/anaconda/bin/python
"""
nonrefco.py is a python program designed to calculate the Nonlinear reflection
coefficient for silicon surfaces. It works in conjunction with the matrix
elements calculated with ABINIT, an open source ab initio software, and
TINIBA, our in-house optical calculation software.

The work codified in this software can be found in Phys.Rev.B66, 195329(2002).
Please refer to "Strain induced SHG" manuscript (BMS) for equation references
unless stated otherwise.

Experimental gap for bulk Si = 3.4 according to Landolt-Boernstein.

For gnuplot:
abso(w,x,y,z)=sqrt((w+y)**2+(x+z)**2)
"""

import sys
import math
import numpy as np
from scipy import constants

# max energy
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

def epsilon(in_file):
    """
    Loads Chi^(1) file, unpacks columns, and creates numpy array.
    """
    data = np.loadtxt(in_file, unpack=True, skiprows=1)
    real = (data[1] + data[3] + data[5])/3 # real average
    imag = (data[2] + data[4] + data[6])/3 # imag average
    coma = real + 1j * imag # complex average
    comx = data[1] + 1j * data[2] # complex x component
    comy = data[3] + 1j * data[4] # complex y component
    comz = data[5] + 1j * data[6] # complex z component
    epsa = 1 + (4 * constants.pi * coma) # epsilon average
    epsx = 1 + (4 * constants.pi * comx) # epsilon x component
    epsy = 1 + (4 * constants.pi * comy) # epsilon y component
    epsz = 1 + (4 * constants.pi * comz) # epsilon z component
    eps = [epsa, epsx, epsy, epsz]
    return eps

def shgcomp(in_file):
    """
    Loads shg Chi^(2) files, unpacks columns, sums 1w and 2w for real and
    imag, and combines into complex numpy array.
    """
    data = np.loadtxt(in_file, unpack=True, skiprows=1)
    comp = (data[1] + data[3]) + 1j * (data[2] + data[4])
    shg = comp[:MAXE]
    return shg

def fresnel(pol, i, j, freq):
    """
    Generic Fresnel factors
    """
    ki = eval("k" + i + freq)
    kj = eval("k" + j + freq)
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
ONEE = np.linspace(0.01, 10, MAXE) # 1w energy array
HBAR = constants.value("Planck constant over 2 pi in eV s")
EPS0 = constants.epsilon_0 # (F/m)
LSPEED = constants.c # (m/s)
N0E = 1.11e19 # for silicon, from PRB 81, 3781 Ref. 17 (V/m^2)
PM2TOM2 = 1e-24 # pm^2 to m^2
M2TOCM2 = 1e4 # m^2 to cm^2
TINIBASCALE = 1e6 # for scaling chi in 1e6 (pm^2/V)
SCALE = 1e20 # for R in 1e-20 (cm^2/W)
PREFACTOR = 1 / (2 * EPS0 * HBAR**2 * LSPEED**3 * math.cos(THETAIN)**2)

# creates epsilons from chi1 responses
epsl = epsilon(PARAM['chil'])
epsb = epsilon(PARAM['chib'])

# epsilons
epsv1w = 1
epsv2w = 1
epsb1w = epsb[0][:MAXE]
epsb2w = epsb[0][1::2]
if MODE == "2layer":
    epsl1w = epsb[0][:MAXE]
    epsl2w = 1
elif MODE == "3layer":
    epsl1w = epsl[0][:MAXE]
    epsl2w = epsl[0][1::2]

# refraction indices
nb = np.sqrt(epsb1w)
Nb = np.sqrt(epsb2w)
nl = np.sqrt(epsl1w)
Nl = np.sqrt(epsl2w)

# wave vectors for 1w and 2w
kv1w = np.sqrt(epsv1w - (math.sin(THETAIN) ** 2))
kv2w = np.sqrt(epsv2w - (math.sin(THETAIN) ** 2))
kb1w = np.sqrt(epsb1w - (math.sin(THETAIN) ** 2))
kb2w = np.sqrt(epsb2w - (math.sin(THETAIN) ** 2))
kl1w = np.sqrt(epsl1w - (math.sin(THETAIN) ** 2))
kl2w = np.sqrt(epsl2w - (math.sin(THETAIN) ** 2))

## fresnel factors for 1w and 2w, s and p polarizations
if MODE == "2layer":
    ell1w = "b"
    ell2w = "v"
elif MODE == "3layer":
    ell1w = "l"
    ell2w = "l"
tvls = fresnel("s", "v", ell1w, "1w")
tvlp = fresnel("p", "v", ell1w, "1w")
tlbs = fresnel("s", ell1w, "b", "1w")
tlbp = fresnel("p", ell1w, "b", "1w")
Tvls = fresnel("s", "v", ell2w, "2w")
Tvlp = fresnel("p", "v", ell2w, "2w")
Tlbs = fresnel("s", ell2w, "b", "2w")
Tlbp = fresnel("p", ell2w, "b", "2w")

# loads chi2, converts to m^2/V
ZZZ = (TINIBASCALE * PM2TOM2 * shgcomp(PARAM['zzz']))
ZXX = (TINIBASCALE * PM2TOM2 * shgcomp(PARAM['zxx']))
XXZ = (TINIBASCALE * PM2TOM2 * shgcomp(PARAM['xxz']))
XXX = (TINIBASCALE * PM2TOM2 * shgcomp(PARAM['xxx']))

# r factors for different input and output polarizations
rpP = epsb2w * (math.sin(THETAIN)) * \
     (epsb1w**2 * (math.sin(THETAIN))**2 * ZZZ + \
                     epsl1w**2 * kb1w**2 * ZXX) - \
     epsl1w * epsl2w * kb1w * kb2w * \
       (2 * epsb1w * (math.sin(THETAIN)) * XXZ + \
                           epsl1w * kb1w * XXX * math.cos(3 * PHI))
rpS = -(kb1w ** 2) * (epsl1w ** 2) * XXX * math.sin(3 * PHI)
rsP = (epsb2w * (math.sin(THETAIN)) * ZXX) - \
      (epsl2w * kb2w * XXX * math.cos(3 * PHI))
rsS = XXX * math.sin(3 * PHI)

# fresnel factors multiplied out for ease of debugging
GammapP = ((Tvlp * Tlbp)/(epsl2w * Nb)) * \
          ((tvlp * tlbp)/(epsl1w * nb))**2
GammapS = Tvls * Tlbs * ((tvlp * tlbp)/(epsl1w * nb))**2
GammasP = ((Tvlp * Tlbp)/(epsl2w * Nb)) * (tvls * tlbs)**2
GammasS = Tvls * Tlbs * (tvls * tlbs)**2

# R factors for different input and output polarizations (in cm^2/W)
RpP = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((np.sqrt(Nl)/nl) * GammapP * rpP)**2
RpS = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((np.sqrt(Nl)/nl) * GammapS * rpS)**2
RsP = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((np.sqrt(Nl)/nl) * GammasP * rsP)**2
RsS = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
      np.absolute((np.sqrt(Nl)/nl) * GammasS * rsS)**2

# RUN THESE IN 3 LAYER MODE
# These are RpP tests that evaluate
tvbp = fresnel("p", "v", "b", "1w")
Tvbp = fresnel("p", "v", "b", "2w")
# P(2w) and Ew in the bulk (Case 5)
rpPbulk = ((math.sin(THETAIN)**3) * ZZZ) \
        + ((kb1w**2) * math.sin(THETAIN) * ZXX) \
        - (2 * kb1w * kb2w * math.sin(THETAIN) * XXZ) \
        - ((kb1w**2) * kb2w * XXX * math.cos(3 * PHI))
GammapPbulk = (Tvbp * (tvbp**2))/(epsb1w * Nb)
RpPbulk = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
          np.absolute((np.sqrt(Nb)/nb) * GammapPbulk * rpPbulk)**2
# P(2w) and Ew in the vacuum (Case 3)
rpPvacuum = ((epsb1w**2) * epsb2w * (math.sin(THETAIN)**3) * ZZZ) \
          + (epsb2w * (kb1w**2) * math.sin(THETAIN) * ZXX) \
          - (2 * epsb1w * kb1w * kb2w * math.sin(THETAIN) * XXZ) \
          - ((kb1w**2) * kb2w * XXX * math.cos(3 * PHI))
GammapPvacuum = (Tvbp * (tvbp**2))/(epsb1w * Nb)
RpPvacuum = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
            np.absolute(1 * GammapPvacuum * rpPvacuum)**2
# P(2w) in l and Ew in the bulk (Case 4)
rpPlb = (epsb2w * (math.sin(THETAIN)**3) * ZZZ) \
      + (epsb2w * (kb1w**2) * math.sin(THETAIN) * ZXX) \
      - (2 * epsl2w * kb1w * kb2w * math.sin(THETAIN) * XXZ) \
      - (epsl2w * (kb1w**2) * kb2w * XXX * math.cos(3 * PHI))
GammapPlb = (Tvlp * Tlbp * (tvbp**2))/(epsl2w * epsb1w * Nb)
RpPlb = SCALE * M2TOCM2 * PREFACTOR * (ONEE ** 2) * \
            np.absolute((np.sqrt(Nl)/nb) * GammapPvacuum * rpPvacuum)**2


# creates columns for 2w and R factors and writes to file
NRC = np.column_stack((ONEE, RpP, RpS, RsP, RsS, RpPbulk, RpPvacuum, RpPlb))
OUTF = PARAM['output']
np.savetxt(OUTF, NRC,
           fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e'),
           delimiter='    ', header='RiF 1e-20 (cm^2/W) 1w RpP RpS RsP RsS')
