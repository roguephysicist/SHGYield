import numpy as np
from scipy import constants, ndimage
from scipy.interpolate import InterpolatedUnivariateSpline
np.seterr(divide='ignore', invalid='ignore', over='ignore') # ignores overflow and divide-by-zero

def geneps(real, imag, sigma):
    redata = broad(real, sigma)
    imdata = broad(imag, sigma)
    chi = redata + 1j * imdata  # complex 1w
    eps = 1 + (4 * np.pi * chi)
    return eps

def spline(data, energy_old, energy_new):
    respl = InterpolatedUnivariateSpline(energy, data.real, ext=2)
    imspl = InterpolatedUnivariateSpline(energy, data.imag, ext=2)
    return respl(energy) + 1j * imspl(energy)

def broad(data, sigma):
    broadened = ndimage.filters.gaussian_filter(data, sigma)
    return broadened

# INFILE = '/Users/sma/Developer/ferroelectrics/In2Se3/same/fe-zbp/1ql/tiniba/res/chi1-vnl.kk_xx_yy_zz_3364_40-nospin_scissor_0.00_Nc_28_wi_0.00_wf_10.00_nw_1000_g_0.01'
INFILE = '/Users/sma/Developer/SHGYield/example/chi1-linear/SiBulk-chi1-xx_yy_zz'
SIGMA = 5

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

np.savetxt('linear.dat',
           np.column_stack((energy, wavelength, epsavg.real, epsavg.imag, np.abs(epsavg), index.real, index.imag, alpha, reflectivity, transmission)),
           fmt='%05.2f  % 14.12e  % 14.12e  % 14.12e  % 14.12e  % 14.12e  % 14.12e  % 14.12e  % 14.12e  % 14.12e',
           header='E     lambda               re[eps]              im[eps]              abs[eps]             n                    k                    alpha                R                    T')
