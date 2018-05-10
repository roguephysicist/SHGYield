#!/usr/bin/env python
"""
shgyield.py is a python program designed to calculate the nonlinear reflection
coefficient for semiconductor surfaces. It works in conjunction with the matrix
elements calculated with ABINIT, an open source ab initio software, and TINIBA,
our in-house optical calculation software. For a complete overview of the
theory, see PRB 94, 115314 (2016).

TODO:
* Dict comprehension to reduce array creation
* Switch to matrix expressions
* Add broadening for input components

required packages:
`sys, yaml, numpy, scipy`

usage:
`python shgyield.py input.yml`
"""

import sys
import yaml
import numpy as np
from scipy import constants, ndimage
from scipy.interpolate import InterpolatedUnivariateSpline

np.seterr(divide='ignore', invalid='ignore', over='ignore') # ignores overflow and divide-by-zero


def parse_input(infile):
    '''
    Parses the YAML input file specified on the command line. Returns dictionary
    with variables and values from the input file.
    '''
    with open(infile) as data_file:
        params = yaml.load(data_file)
    return params


def epsload(infile, norm):
    '''
    Reads calculated chi1 file, that is organized as Energy(1w) Re[chi_xx]
    Im[chi_xx] Re[chi_yy] Im[chi_yy] Re[chi_zz] Im[chi_zz]. Normalizes chi1
    according to PRB 92, 245308 (2015), and averages chi1 = (chi1^{xx} +
    chi1^{yy} + chi1^{zz})/3. Creates two 'InterpolatedUnivariateSpline' objects
    (splines that pass through EVERY point of the original data), one each for
    the real and imaginary parts. It then generates new arrays based on the
    input energy range or value. Finally, converts to epsilon = 1 + 4*pi*chi1,
    and returns numpy arrays with the 1w and 2w epsilons.
    '''
    freq, rexx, imxx, reyy, imyy, rezz, imzz = np.loadtxt(infile, unpack=True)
    real = broad(((rexx + reyy + rezz)/3) * norm, PARAM['chi1']['sigma']) # real average
    imag = broad(((imxx + imyy + imzz)/3) * norm, PARAM['chi1']['sigma']) # imag average
    respl = InterpolatedUnivariateSpline(freq, real, ext=2)
    imspl = InterpolatedUnivariateSpline(freq, imag, ext=2)
    chi1w = respl(ENERGY) + 1j * imspl(ENERGY)  # complex average, 1w
    chi2w = respl(2*ENERGY) + 1j * imspl(2*ENERGY)  # complex average, 2w
    eps1w = 1 + (4 * np.pi * chi1w)
    eps2w = 1 + (4 * np.pi * chi2w)
    return eps1w, eps2w


def shgload(infile, norm):
    '''
    Loads chi^{abc} listed in the input file, with the following columns:
    Energy(1w) Re[1w] Im[1w] Re[2w] Im[2w]. Sums 1w and 2w real and imaginary
    parts, and then reates two 'InterpolatedUnivariateSpline' objects (splines
    that pass through EVERY point of the original data), one each for the real
    and imaginary parts. It then generates new arrays based on the input energy
    range or value. Returns numpy array with the appropriate scale and units.
    '''
    freq, re1w, im1w, re2w, im2w = np.loadtxt(infile, unpack=True)
    real = broad((re1w + re2w) * norm, PARAM['chi2']['sigma'])
    imag = broad((im1w + im2w) * norm, PARAM['chi2']['sigma'])
    respl = InterpolatedUnivariateSpline(freq, real, ext=2)
    imspl = InterpolatedUnivariateSpline(freq, imag, ext=2)
    comp = respl(ENERGY) + 1j * imspl(ENERGY)
    chi2 = SCALE * PM2TOM2 * comp
    return chi2


def wvec(eps):
    '''
    Wave vector, where w = sqrt[epsilon - sin(theta)^2].
    '''
    wavevector = np.sqrt(eps - (np.sin(THETA0)**2))
    return wavevector


def frefs(epsi, epsj):
    '''
    Generic reflection fresnel factors, see Eqs. (13) and (14) of PRB 94, 115314 (2016).
    '''
    factor = (wvec(epsi) - wvec(epsj))/(wvec(epsi) + wvec(epsj))
    return factor


def frefp(epsi, epsj):
    '''
    Generic reflection fresnel factors, see Eqs. (13) and (14) of PRB 94, 115314 (2016).
    '''
    factor = ((wvec(epsi) * epsj) - (wvec(epsj) * epsi))/\
             ((wvec(epsi) * epsj) + (wvec(epsj) * epsi))
    return factor


def ftrans(epsi, epsj):
    '''
    s-polarized transmission fresnel factors , see Eqs. (13) and (14) of PRB 94, 115314 (2016).
    '''
    factor = (2 * wvec(epsi))/(wvec(epsi) + wvec(epsj))
    return factor


def ftranp(epsi, epsj):
    '''
    p-polarized transmission fresnel factors, see Eqs. (13) and (14) of PRB 94, 115314 (2016).
    '''
    factor = (2 * wvec(epsi) * np.sqrt(epsi * epsj))/\
             (wvec(epsi) * epsj + wvec(epsj) * epsi)
    return factor


def mrc(fres1, fres2, phase1, phase2, rem):
    '''
    Multiple reflection coefficient, see Eqs. (18) and (21) of PRB 94, 115314 (2016).
    '''
    coeff = (fres1 * np.exp(1j * phase1)) / (1 + (fres2 * fres1 * np.exp(1j * phase2))) * rem
    return coeff


def rad_pp(phi):
    '''
    r factors for different input and output polarizations. See Eqs. (50), (55),
    (60), and (65) of PRB 94, 115314 (2016). Returns rpP, rpS, rsP, and rsS.
    '''
    pre = (ftranp(EPS['2w']['v'], EPS['2w']['l'])/INDEX['Nl']) * \
          (ftranp(EPS['1w']['v'], EPS['1w']['l'])/INDEX['nl'])**2
    ### r_{pP}
    rpp = - (FRES['2w']['p-'] * FRES['1w']['p-']**2 \
            * wvec(EPS['1w']['l'])**2 * wvec(EPS['2w']['l']) \
            * np.cos(phi)**3 * CHI2['xxx']) \
          - (2 * FRES['2w']['p-'] * FRES['1w']['p-']**2 \
            * wvec(EPS['1w']['l'])**2 * wvec(EPS['2w']['l']) \
            * np.sin(phi) * np.cos(phi)**2 * CHI2['xxy']) \
          - (2 * FRES['2w']['p-'] * FRES['1w']['p+'] * FRES['1w']['p-'] \
            * wvec(EPS['1w']['l']) * wvec(EPS['2w']['l']) \
            * np.sin(THETA0) * np.cos(phi)**2 * CHI2['xxz']) \
          - (FRES['2w']['p-'] * FRES['1w']['p-']**2 \
            * wvec(EPS['1w']['l'])**2 * wvec(EPS['2w']['l']) \
            * np.sin(phi)**2 * np.cos(phi) * CHI2['xyy']) \
          - (2 * FRES['2w']['p-'] * FRES['1w']['p+'] * FRES['1w']['p-'] \
            * wvec(EPS['1w']['l']) * wvec(EPS['2w']['l']) \
            * np.sin(THETA0) * np.sin(phi) * np.cos(phi) * CHI2['xyz']) \
          - (FRES['2w']['p-'] * FRES['1w']['p+']**2 \
            * wvec(EPS['2w']['l']) \
            * np.sin(THETA0)**2 * np.cos(phi) * CHI2['xzz']) \
          - (FRES['2w']['p-'] * FRES['1w']['p-']**2 \
            * wvec(EPS['1w']['l'])**2 * wvec(EPS['2w']['l']) \
            * np.sin(phi) * np.cos(phi)**2 * CHI2['yxx']) \
          - (2 * FRES['2w']['p-'] * FRES['1w']['p-']**2 \
            * wvec(EPS['1w']['l'])**2 * wvec(EPS['2w']['l']) \
            * np.sin(phi)**2 * np.cos(phi) * CHI2['yxy']) \
          - (2 * FRES['2w']['p-'] * FRES['1w']['p+'] * FRES['1w']['p-'] \
            * wvec(EPS['1w']['l']) * wvec(EPS['2w']['l']) \
            * np.sin(THETA0) * np.sin(phi) * np.cos(phi) * CHI2['yxz']) \
          - (FRES['2w']['p-'] * FRES['1w']['p-']**2 \
            * wvec(EPS['1w']['l'])**2 * wvec(EPS['2w']['l']) \
            * np.sin(phi)**3 * CHI2['yyy']) \
          - (2 * FRES['2w']['p-'] * FRES['1w']['p+'] * FRES['1w']['p-'] \
            * wvec(EPS['1w']['l']) * wvec(EPS['2w']['l']) \
            * np.sin(THETA0) * np.sin(phi)**2 * CHI2['yyz']) \
          - (FRES['2w']['p-'] * FRES['1w']['p+']**2 \
            * wvec(EPS['2w']['l']) \
            * np.sin(THETA0)**2 * np.sin(phi) * CHI2['yzz']) \
          + (FRES['2w']['p+'] * FRES['1w']['p-']**2 \
            * wvec(EPS['1w']['l'])**2 \
            * np.sin(THETA0) * np.cos(phi)**2 * CHI2['zxx']) \
          + (2 * FRES['2w']['p+'] * FRES['1w']['p+'] * FRES['1w']['p-'] \
            * wvec(EPS['1w']['l']) \
            * np.sin(THETA0)**2 * np.cos(phi) * CHI2['zxz']) \
          + (2 * FRES['2w']['p+'] * FRES['1w']['p-']**2 \
            * wvec(EPS['1w']['l'])**2 \
            * np.sin(THETA0) * np.sin(phi) * np.cos(phi) * CHI2['zxy']) \
          + (FRES['2w']['p+'] * FRES['1w']['p-']**2 \
            * wvec(EPS['1w']['l'])**2 \
            * np.sin(THETA0) * np.sin(phi)**2 * CHI2['zyy']) \
          + (2 * FRES['2w']['p+'] * FRES['1w']['p+'] * FRES['1w']['p-'] \
            * wvec(EPS['1w']['l']) \
            * np.sin(THETA0)**2 * np.sin(phi) * CHI2['zyz']) \
          + (FRES['2w']['p+'] * FRES['1w']['p+']**2 * np.sin(phi)**3 * CHI2['zzz'])
    return pre*rpp


def rad_ps(phi):
    '''
    r factors for different input and output polarizations. See Eqs. (50), (55),
    (60), and (65) of PRB 94, 115314 (2016). Returns rpP, rpS, rsP, and rsS.
    '''
    pre = (ftrans(EPS['2w']['v'], EPS['2w']['l']) * FRES['2w']['s+']) * \
          (ftranp(EPS['1w']['v'], EPS['1w']['l'])/INDEX['nl'])**2
    ### r_{pS}
    rps = - (FRES['1w']['p-']**2 * wvec(EPS['1w']['l'])**2 \
            * np.sin(phi) * np.cos(phi)**2 * CHI2['xxx']) \
          - (2 * FRES['1w']['p-']**2 * wvec(EPS['1w']['l'])**2 \
            * np.sin(phi)**2 * np.cos(phi) * CHI2['xxy']) \
          - (2 * FRES['1w']['p+'] * FRES['1w']['p-'] * wvec(EPS['1w']['l']) \
            * np.sin(THETA0) * np.sin(phi) * np.cos(phi) * CHI2['xxz']) \
          - (FRES['1w']['p-']**2 * wvec(EPS['1w']['l'])**2 \
            * np.sin(phi)**3 * CHI2['xyy']) \
          - (2 * FRES['1w']['p+'] * FRES['1w']['p-'] * wvec(EPS['1w']['l']) \
            * np.sin(THETA0) * np.sin(phi)**2 * CHI2['xyz']) \
          - (FRES['1w']['p+']**2 \
            * np.sin(THETA0)**2 * np.sin(phi) * CHI2['xzz']) \
          + (FRES['1w']['p-']**2 * wvec(EPS['1w']['l'])**2 \
            * np.cos(phi)**3 * CHI2['yxx']) \
          + (2 * FRES['1w']['p-']**2 * wvec(EPS['1w']['l'])**2 \
            * np.sin(phi) * np.cos(phi)**2 * CHI2['yxy']) \
          + (2 * FRES['1w']['p+'] * FRES['1w']['p-'] * wvec(EPS['1w']['l']) \
            * np.sin(THETA0) * np.cos(phi)**2 * CHI2['yxz']) \
          + (FRES['1w']['p-']**2 * wvec(EPS['1w']['l'])**2 \
            * np.sin(phi)**2 * np.cos(phi) * CHI2['yyy']) \
          + (2 * FRES['1w']['p+'] * FRES['1w']['p-'] * wvec(EPS['1w']['l']) \
            * np.sin(THETA0) * np.sin(phi) * np.cos(phi) * CHI2['yyz']) \
          + (FRES['1w']['p+']**2 \
            * np.sin(THETA0)**2 * np.cos(phi) * CHI2['yzz'])
    return pre*rps


def rad_sp(phi):
    '''
    r factors for different input and output polarizations. See Eqs. (50), (55),
    (60), and (65) of PRB 94, 115314 (2016). Returns rpP, rpS, rsP, and rsS.
    '''
    pre = (ftranp(EPS['2w']['v'], EPS['2w']['l'])/INDEX['Nl']) * \
          (ftrans(EPS['1w']['v'], EPS['1w']['l']) * FRES['1w']['s+'])**2
    ### r_{sP}
    rsp = - (FRES['2w']['p-'] * wvec(EPS['2w']['l']) \
            * np.sin(phi)**2 * np.cos(phi) * CHI2['xxx']) \
          + (FRES['2w']['p-'] * wvec(EPS['2w']['l']) \
            * 2 * np.sin(phi) * np.cos(phi)**2 * CHI2['xxy']) \
          - (FRES['2w']['p-'] * wvec(EPS['2w']['l']) \
            * np.cos(phi)**3 * CHI2['xyy']) \
          - (FRES['2w']['p-'] * wvec(EPS['2w']['l']) \
            * np.sin(phi)**3 * CHI2['yxx']) \
          + (FRES['2w']['p-'] * wvec(EPS['2w']['l']) \
            * 2 * np.sin(phi)**2 * np.cos(phi) * CHI2['yxy']) \
          - (FRES['2w']['p-'] * wvec(EPS['2w']['l']) \
            * np.sin(phi) * np.cos(phi)**2 * CHI2['yyy']) \
          + (FRES['2w']['p+'] \
            * np.sin(THETA0) * np.sin(phi)**2 * CHI2['zxx']) \
          - (FRES['2w']['p+'] \
            * np.sin(THETA0) * 2 * np.sin(phi) * np.cos(phi) * CHI2['zxy']) \
          + (FRES['2w']['p+'] \
            * np.sin(THETA0) * np.cos(phi)**2 * CHI2['zyy'])
    return pre*rsp


def rad_ss(phi):
    '''
    r factors for different input and output polarizations. See Eqs. (50), (55),
    (60), and (65) of PRB 94, 115314 (2016). Returns rpP, rpS, rsP, and rsS.
    '''
    pre = (ftrans(EPS['2w']['v'], EPS['2w']['l']) * FRES['2w']['s+']) * \
          (ftrans(EPS['1w']['v'], EPS['1w']['l']) * FRES['1w']['s+'])**2
    ### r_{sS}
    rss = - (np.sin(phi)**3 * CHI2['xxx']) \
          + (2 * np.sin(phi)**2 * np.cos(phi) * CHI2['xxy']) \
          - (np.sin(phi) * np.cos(phi)**2 * CHI2['xyy']) \
          + (np.sin(phi)**2 * np.cos(phi) * CHI2['yxx']) \
          + (np.cos(phi)**3 * CHI2['yyy']) \
          - (2 * np.sin(phi) * np.cos(phi)**2 * CHI2['yxy'])
    return pre*rss


def shgyield(factor):
    '''
    Calculates the final broadened SHG yield, ready to be written to file.
    See Eq. (38) of PRB 94, 115314 (2016).
    '''
    rif = M2TOCM2 * PREFACTOR * (ENERGY ** 2) * np.absolute((1/INDEX['nl']) * factor)**2
    broadened = broad(rif, SIGMA)
    return broadened


def broad(target, sigma):
    '''
    A function for applying Gaussian broadening on the final output data.
    '''
    data = ndimage.filters.gaussian_filter(target, sigma)
    return data


def savefile(file, data):
    '''
    Saves specified dataset to file with the following columns:
    Energy(1w)    R_{pP}    R_{pS}    R_{sP}    R_{sS}
    '''
    np.savetxt(file, data, fmt=('%07.4f', '%.8e', '%.8e', '%.8e', '%.8e', '%05.1f'),
               delimiter='    ',
               header='RiF (cm^2/W)\n1w(eV)   RpP'+15*" "+\
                      'RpS'+15*" "+'RsP'+15*" "+'RsS'+15*" "+'phi(deg)')





################################################################################
########################### BEGIN INITIALIZATION ###############################
################################################################################

## Initialization: Parse input file, establish relevant modes.
PARAM = parse_input(sys.argv[1])    # Parses input file
MODE = PARAM['mode']                # Establishes the 'layer model' to be used;
                                    # see PRB 93, 235304 (2016).
ENINPUT = PARAM['energy']           # Establishes energy range over which the
                                    # yield will be calculated.
if isinstance(ENINPUT, list):
    ENERGY = np.linspace(*ENINPUT)
elif isinstance(ENINPUT, (int, float)):
    ENERGY = float(ENINPUT)

PHIINPUT = PARAM['parameters']['phi'] # Converts phi to radians
if isinstance(PHIINPUT, list):
    PHI = np.arange(*PHIINPUT)
elif isinstance(PHIINPUT, (int, float)):
    PHI = [float(PHIINPUT)]

## Constants and parameters
HBAR = constants.value("Planck constant over 2 pi in eV s")
PLANCK = constants.value("Planck constant in eV s")
EPS0 = constants.epsilon_0                   # \epsilon_{0} in F/m
LSPEED = constants.c                         # Speed of light in m/s
PM2TOM2 = 1e-24                              # Convert from pm^2 to m^2
M2TOCM2 = 1e4                                # Convert from m^2 to cm^2
SCALE = float(PARAM['chi2']['scale'])        # if a scaling factor is included (pm^2/V)
THETA0 = np.radians(PARAM['parameters']['theta']) # Converts theta to radians
SIGMA = PARAM['parameters']['sigma']         # Std. dev. for gaussian broadening
PREFACTOR = 1/(2 * EPS0 * HBAR**2 * LSPEED**3 * np.cos(THETA0)**2)

################################################################################
############################# END INITIALIZATION ###############################
################################################################################


################################################################################
######################### BEGIN LINEAR RESPONSES ###############################
################################################################################

## Linear responses: chi1 and epsilons
CHI1NORM = PARAM['chi1']['norm']            # Normalization for layered chi1
try:
    EPSB = epsload(PARAM['chi1']['chib'], 1)    # Epsilon from chi1, bulk
except (ValueError, OSError, IOError):
    if isinstance(PARAM['chi1']['chib'], (float, int)):
        EPSB = (PARAM['chi1']['chib'].real, PARAM['chi1']['chib'].imag)
EPS = {'1w': {'v': 1, 'b': EPSB[0]}, '2w': {'v': 1, 'b': EPSB[1]}}

## Reflection model, see PRB 93, 235304 (2016).
if MODE == "3-layer": # The incident fields and SHG both occur in the thin layer (l)
    EPSL = epsload(PARAM['chi1']['chil'], CHI1NORM) # Epsilon from chi1, layered, normalized
    EPS['1w']['l'] = EPSL[0]
    EPS['2w']['l'] = EPSL[1]
elif MODE == "2-layer-fresnel": # The incident fields in bulk, SHG in vacuum
    EPS['1w']['l'] = EPS['1w']['b']
    EPS['2w']['l'] = EPS['2w']['v']
elif MODE == "2-layer-bulk": # Both incident fields and SHG in bulk
    EPS['1w']['l'] = EPS['1w']['b']
    EPS['2w']['l'] = EPS['2w']['b']
elif MODE == "2-layer-vacuum": # Both incident fields and SHG in vacuum
    EPS['1w']['l'] = EPS['2w']['v']
    EPS['2w']['l'] = EPS['2w']['v']
elif MODE == "3-layer-hybrid": # Incident field in bulk, SHG in thin layer (l)
    EPSL = epsload(PARAM['chi1']['chil'], CHI1NORM) # Epsilon from chi1, layered, normalized
    EPS['1w']['l'] = EPS['1w']['b']
    EPS['2w']['l'] = EPSL[1]                        # Epsilon for layer, 2w

## Indices of refraction for 1w (n) and 2w (N), where n = sqrt{epsilon}
# INDEX = {key: np.sqrt(value) for key, value in EPS.items()}
INDEX = {'nl': np.sqrt(EPS['1w']['l']), 'Nl': np.sqrt(EPS['2w']['l'])}

################################################################################
########################### END LINEAR RESPONSES ###############################
################################################################################


################################################################################
########################## BEGIN FRESNEL COEFFS ################################
################################################################################

## Fresnel factors and multiple reflections framework. See Eqs. (16), (17),
## (21), (22), (26), and (30) of PRB 94, 115314 (2016).
FRES = {'2w' : {'p' : {'lb' : frefp(EPS['2w']['l'], EPS['2w']['b']),
                       'vl' : frefp(EPS['2w']['v'], EPS['2w']['l'])},
                's' : {'lb' : frefs(EPS['2w']['l'], EPS['2w']['b']),
                       'vl' : frefs(EPS['2w']['v'], EPS['2w']['l'])}},
        '1w' : {'p' : {'lb' : frefp(EPS['1w']['l'], EPS['1w']['b']),
                       'vl' : frefp(EPS['1w']['v'], EPS['1w']['l'])},
                's' : {'lb' : frefs(EPS['1w']['l'], EPS['1w']['b']),
                       'vl' : frefs(EPS['1w']['v'], EPS['1w']['l'])}}}

## Multiple reflections
if PARAM['multiref']['enable'] and MODE == '3-layer':
    THICKNESS = PARAM['multiref']['thickness'] # Thickness d of the thin layer \ell
    DEPTH = PARAM['multiref']['depth']         # Depth d2 of the polarization sheet
    if DEPTH == "average":
        DELTA = 8*np.pi * ((ENERGY * THICKNESS * 1e-9)/(PLANCK * LSPEED)) * wvec(EPS['2w']['l'])
        FRES['2w']['p'] = mrc(FRES['2w']['p']['lb'], FRES['2w']['p']['vl'], DELTA/2, DELTA, np.sin(DELTA/2)/(DELTA/2))
        FRES['2w']['s'] = mrc(FRES['2w']['s']['lb'], FRES['2w']['s']['vl'], DELTA/2, DELTA, np.sin(DELTA/2)/(DELTA/2))
    else:
        DELTA = 8*np.pi * ((ENERGY * float(DEPTH) * 1e-9)/(PLANCK * LSPEED)) * wvec(EPS['2w']['l'])
        FRES['2w']['p'] = mrc(FRES['2w']['p']['lb'], FRES['2w']['p']['vl'], DELTA, DELTA, 1)
        FRES['2w']['s'] = mrc(FRES['2w']['s']['lb'], FRES['2w']['s']['vl'], DELTA, DELTA, 1)
    VARPHI = 4 * np.pi * ((ENERGY * THICKNESS * 1e-9)/(PLANCK * LSPEED)) * wvec(EPS['1w']['l'])
    FRES['1w']['p'] = mrc(FRES['1w']['p']['lb'], FRES['1w']['p']['vl'], VARPHI, VARPHI, 1)
    FRES['1w']['s'] = mrc(FRES['1w']['s']['lb'], FRES['1w']['s']['vl'], VARPHI, VARPHI, 1)
elif not PARAM['multiref']['enable'] or MODE != '3-layer':
    FRES['2w']['p'] = FRES['2w']['p']['lb']
    FRES['2w']['s'] = FRES['2w']['s']['lb']
    FRES['1w']['p'] = FRES['1w']['p']['lb']
    FRES['1w']['s'] = FRES['1w']['s']['lb']
## Final Fresnel factors with multiple reflection
FRES = {'2w' : {'p+' : 1 + FRES['2w']['p'],
                's+' : 1 + FRES['2w']['s'],
                'p-' : 1 - FRES['2w']['p'],
                's-' : 1 - FRES['2w']['s']},
        '1w' : {'p+' : 1 + FRES['1w']['p'],
                's+' : 1 + FRES['1w']['s'],
                'p-' : 1 - FRES['1w']['p'],
                's-' : 1 - FRES['1w']['s']}}

################################################################################
############################ END FRESNEL COEFFS ################################
################################################################################


################################################################################
####################### BEGIN NONLINEAR RESPONSES ##############################
################################################################################

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
CHI2NORM = PARAM['chi2']['norm']            # Normalization for chi2
CHI2 = {}
CHI2_EQUIV = {}
ALL_COMPONENTS = ['xxx', 'xyy', 'xzz', 'xyz', 'xxz', 'xxy',
                  'yxx', 'yyy', 'yzz', 'yyz', 'yxz', 'yxy',
                  'zxx', 'zyy', 'zzz', 'zyz', 'zxz', 'zxy']
for component in ALL_COMPONENTS:
    if component in PARAM['chi2']:
        value = PARAM['chi2'][component]
        try:
            shg = shgload(value, CHI2NORM)
        except (ValueError, OSError, IOError):
            if value == 0:
                shg = 0
            else:
                CHI2_EQUIV[component] = value
                continue
    elif component not in PARAM['chi2']:
        shg = 0
    CHI2[component] = shg
for key, value in CHI2_EQUIV.items():
    equivalence = value.split('-')
    if equivalence[0]:
        CHI2[key] = CHI2[equivalence[0]]
    elif not equivalence[0]:
        CHI2[key] = -1 * CHI2[equivalence[1]]
### CHANGE EQUIVALENCE STUFF TO A FUNCTION, FOR DICT COMPREHENSION DIRECTLY

################################################################################
######################### END NONLINEAR RESPONSES ##############################
################################################################################


################################################################################
######################## BEGIN FINAL SSHG YIELD ################################
################################################################################

## Final SHG yield for different input and output polarizations (in cm^2/W).
## See Eqs. (44) and (38) of PRB 94, 115314 (2016).
DATA = []
for azimuth in PHI:
    angle = np.radians(azimuth) # Converts phi to radians
    rfactors = {'pP': rad_pp(angle), 'pS': rad_ps(angle), 'sP': rad_sp(angle), 'sS': rad_ss(angle)}
    shg = {key: shgyield(val) for key, val in rfactors.items()}
    DATA.append(np.column_stack((ENERGY, shg['pP'], shg['pS'], shg['sP'], shg['sS'], np.full_like(ENERGY, angle))))

FINAL = np.concatenate(DATA)
savefile(PARAM['output'], FINAL)

################################################################################
########################## END FINAL SSHG YIELD ################################
################################################################################
