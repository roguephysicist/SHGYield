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
        factor = (2 * ki * np.sqrt(epsi * epsj)) / (ki * epsj + kj * epsi)
    elif pol == "s":
        factor = (2 * ki) / (ki + kj)
    return factor

# reads input file and establishes mode
param = parse_input()
mode = str(param['mode'])

# angles
thetain = math.radians(float(param['theta']))
phi = math.radians(float(param['phi']))

# constants, conversions, and prefactor
onee = np.linspace(0.01, 10, MAXE) # 1w energy array
hbar = constants.value("Planck constant over 2 pi in eV s")
eps0 = constants.epsilon_0 # (F/m)
lspeed = constants.c # (m/s)
n0e = 1.11e19 # for silicon, from PRB 81, 3781 Ref. 17 (V/m^2)
pm2tom2 = 1e-24 # pm^2 to m^2
m2tocm2 = 1e4 # m^2 to cm^2
tinibascale = 1e6 # for scaling chi in 1e6 (pm^2/V)
scale = 1e20 # for R in 1e-20 (cm^2/W)
prefactor = 1 / (2 * eps0 * hbar**2 * lspeed**3 * math.cos(thetain)**2)

# creates epsilons from chi1 responses
epsl = epsilon(param['chil'])
epsb = epsilon(param['chib'])

# epsilons
epsv1w = 1
epsv2w = 1
epsb1w = epsb[0][:MAXE]
epsb2w = epsb[0][1::2]
if mode == "2layer":
    epsl1w = epsb[0][:MAXE]
    epsl2w = 1
elif mode == "3layer":
    epsl1w = epsl[0][:MAXE]
    epsl2w = epsl[0][1::2]

# refraction indices
nb = np.sqrt(epsb1w)
Nb = np.sqrt(epsb2w)
nl = np.sqrt(epsl1w)
Nl = np.sqrt(epsl2w)

# wave vectors for 1w and 2w
kv1w = np.sqrt(epsv1w - (math.sin(thetain) ** 2))
kv2w = np.sqrt(epsv2w - (math.sin(thetain) ** 2))
kb1w = np.sqrt(epsb1w - (math.sin(thetain) ** 2))
kb2w = np.sqrt(epsb2w - (math.sin(thetain) ** 2))
kl1w = np.sqrt(epsl1w - (math.sin(thetain) ** 2))
kl2w = np.sqrt(epsl2w - (math.sin(thetain) ** 2))

## fresnel factors for 1w and 2w, s and p polarizations
if mode == "2layer":
    ell1w = "b"
    ell2w = "v"
elif mode == "3layer":
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
tvbp = fresnel("p", "v", "b", "1w")
Tvbp = fresnel("p", "v", "b", "2w")

# loads chi2, converts to m^2/V
zzz = (tinibascale * pm2tom2 * shgcomp(param['zzz']))
zxx = (tinibascale * pm2tom2 * shgcomp(param['zxx']))
xxz = (tinibascale * pm2tom2 * shgcomp(param['xxz']))
xxx = (tinibascale * pm2tom2 * shgcomp(param['xxx']))

# r factors for different input and output polarizations
rpp = epsb2w * (math.sin(thetain)) * \
     (epsb1w**2 * (math.sin(thetain))**2 * zzz + \
                     epsl1w**2 * kb1w**2 * zxx) - \
     epsl1w * epsl2w * kb1w * kb2w * \
       (2 * epsb1w * (math.sin(thetain)) * xxz + \
                           epsl1w * kb1w * xxx * math.cos(3 * phi))
rps = -(kb1w ** 2) * (epsl1w ** 2) * xxx * math.sin(3 * phi)
rsp = (epsb2w * (math.sin(thetain)) * zxx) - \
      (epsl2w * kb2w * xxx * math.cos(3 * phi))
rss = xxx * math.sin(3 * phi)
## 2w in the bulk
brpp = ((math.sin(thetain)**3) * zzz) \
     + ((kb1w**2) * math.sin(thetain) * zxx) \
     - (2 * kb1w * kb2w * math.sin(thetain) * xxz) \
     - ((kb1w**2) * kb2w * xxx * math.cos(3 * phi))
## 1w in the vacuum
vrpp = ((epsb1w**2) * epsb2w * (math.sin(thetain)**3) * zzz) \
     + (epsb2w * (kb1w**2) * math.sin(thetain) * zxx) \
     - (2 * epsb1w * kb1w * kb2w * math.sin(thetain) * xxz) \
     - ((kb1w**2) * kb2w * xxx * math.cos(3 * phi))

# fresnel factors multiplied out for ease of debugging
gammapp = ((Tvlp * Tlbp)/(epsl2w * Nb)) * \
          ((tvlp * tlbp)/(epsl1w * nb))**2
gammaps = Tvls * Tlbs * ((tvlp * tlbp)/(epsl1w * nb))**2
gammasp = ((Tvlp * Tlbp)/(epsl2w * Nb)) * (tvls * tlbs)**2
gammass = Tvls * Tlbs * (tvls * tlbs)**2
## 2w in the bulk
bgammapp = (Tvbp * (tvbp**2))/(epsb1w * Nb)
## 1w in the vacuum
vgammapp = (Tvbp/Nb) * (tvbp/nb)**2

# R factors for different input and output polarizations (in cm^2/W)
Rpp = scale * m2tocm2 * prefactor * (onee ** 2) * np.absolute(gammapp * rpp)**2
Rps = scale * m2tocm2 * prefactor * (onee ** 2) * np.absolute(gammaps * rps)**2
Rsp = scale * m2tocm2 * prefactor * (onee ** 2) * np.absolute(gammasp * rsp)**2
Rss = scale * m2tocm2 * prefactor * (onee ** 2) * np.absolute(gammass * rss)**2
bRpp = scale * m2tocm2 * prefactor * (onee ** 2) * np.absolute(bgammapp * brpp)**2
vRpp = scale * m2tocm2 * prefactor * (onee ** 2) * np.absolute(vgammapp * vrpp)**2

# creates columns for 2w and R factors and writes to file
nrc = np.column_stack((onee, Rpp, Rps, Rsp, Rss, bRpp, vRpp))
outf = param['output']
np.savetxt(outf, nrc, fmt=('%05.2f', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e', '%.14e'),
           delimiter='    ', header='RiF in 1e-20 (cm^2/W)\n\
           2w     Rpp' + 21*' ' + 'Rps' + 21*' ' + 'Rsp' + 21*' ' + 'Rss')
