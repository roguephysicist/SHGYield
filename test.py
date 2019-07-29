import sys
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

ENERGY = 3.02

def chi(real, imag, ener):
    return 1 + 4*np.pi*(real(ener) + 1j*imag(ener))


freq, chi_real, chi_imag = np.loadtxt('example/chi1-linear/SiH1x1-chi1-xx', unpack=True)

respl = InterpolatedUnivariateSpline(freq, chi_real)
imspl = InterpolatedUnivariateSpline(freq, chi_imag)

VAR = chi(respl, imspl, 3.0255)

print(VAR)
