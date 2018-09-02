import sys
import numpy as np
from scipy import constants, ndimage
from scipy.interpolate import InterpolatedUnivariateSpline

np.seterr(divide='ignore', invalid='ignore', over='ignore') # ignores overflow and divide-by-zero

def geneps(real, imag, sigma):
    redata = broad(real, sigma)
    imdata = broad(imag, sigma)
    chi = (redata + 1j * imdata) * NORM  # complex 1w
    eps = 1 + (4 * np.pi * chi)
    return eps

def spline(data, energy_old, energy_new):
    respl = InterpolatedUnivariateSpline(energy, data.real, ext=2)
    imspl = InterpolatedUnivariateSpline(energy, data.imag, ext=2)
    return respl(energy) + 1j * imspl(energy)

def broad(data, sigma):
    broadened = ndimage.filters.gaussian_filter(data, sigma)
    return broadened

INFILE = sys.argv[1]
OUTFILE = sys.argv[2]
# NORM = 18.27647196771197
NORM = 7.5182259762807755
SIGMA = 0

energy, rexx, imxx, reyy, imyy, rezz, imzz = np.loadtxt(INFILE, unpack=True)
wavelength = constants.c*1e9*constants.value("Planck constant in eV s")/energy

epsxx = geneps(rexx, imxx, SIGMA)
epsyy = geneps(reyy, imyy, SIGMA)
epszz = geneps(rezz, imzz, SIGMA)
epsavg = np.mean((epsxx, epsyy), axis=0)

index = np.sqrt(epsavg)
alpha = (2/(constants.c*100*constants.value("Planck constant over 2 pi in eV s"))) * energy * index.imag
reflectivity = np.abs((1 - index)/(1 + index))**2
transmission = 1 - reflectivity

np.savetxt(OUTFILE,
           np.column_stack((energy, wavelength, epsavg.real, epsavg.imag, np.abs(epsavg), index.real, index.imag, alpha, reflectivity, transmission)),
           fmt='%05.2f  % 14.12e  % 14.12e  % 14.12e  % 14.12e  % 14.12e  % 14.12e  % 14.12e  % 14.12e  % 14.12e',
           header='E     lambda               re[eps]              im[eps]              abs[eps]             n                    k                    alpha                R                    T')
