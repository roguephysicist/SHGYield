#!/usr/bin/env python
"""
shgyield.py is a python module for exploring the SHG optical response of
materials. It is well suited for 2D-materials, surfaces, bulks, and
metamaterials. For a complete overview of the theory, see PRB 94, 115314 (2016).

required packages:
`numpy, scipy`
"""

import numpy as np
from scipy import constants, ndimage
from scipy.interpolate import InterpolatedUnivariateSpline

np.seterr(divide='ignore', invalid='ignore', over='ignore') # ignores overflow and divide-by-zero


def broad(data, sigma):
    ''' applies Gaussian broadening to real number '''
    return ndimage.filters.gaussian_filter(data, sigma)


def broadC(data, sigma):
    ''' applies Gaussian broadening to complex and returns complex '''
    real = ndimage.filters.gaussian_filter(data.real, sigma)
    imag = ndimage.filters.gaussian_filter(data.imag, sigma)
    return real + 1j*imag


def spline(data, energy):
    ''' returns spline object '''
    return InterpolatedUnivariateSpline(energy, data, ext=2)


def splineC(data, energy):
    ''' creates a spline for complex and returns tuple with spline objects '''
    real = InterpolatedUnivariateSpline(energy, data.real, ext=2)
    imag = InterpolatedUnivariateSpline(energy, data.imag, ext=2)
    return (real, imag)


def splineEPS(data, energy):
    ''' IMPROVE: creates a spline for EPS, returns 1w and 2w '''
    splines = {key: splineC(val['data'], val['energy']) for key, val in data.items()}
    new1w = {key: val[0](energy) + 1j*val[1](energy) for key, val in splines.items()}
    new2w = {key: val[0](2*energy) + 1j*val[1](2*energy) for key, val in splines.items()}
    return (new1w, new2w)


def avgEPS(data):
    ''' IMPROVE: averages over all values in a dict '''
    return np.mean(list(data.values()), axis=0)


def rotate(chi2, gamma):
    ## in-plane rotatation for chi2 tensor components
    chi2rot = {
        'xxx' : + np.sin(gamma)**3*chi2['xxx'] \
                + np.sin(gamma)*np.cos(gamma)**2*chi2['xyy'] \
                - 2*np.sin(gamma)**2*np.cos(gamma)*chi2['xxy'] \
                - np.sin(gamma)**2*np.cos(gamma)*chi2['yxx'] \
                - np.cos(gamma)**3*chi2['yyy'] \
                + 2*np.sin(gamma)*np.cos(gamma)**2*chi2['yxy'],
        'xyy' : + np.sin(gamma)*np.cos(gamma)**2*chi2['xxx'] \
                + np.sin(gamma)**3*chi2['xyy'] \
                + 2*np.sin(gamma)**2*np.cos(gamma)*chi2['xxy'] \
                - np.cos(gamma)**3*chi2['yxx'] \
                - np.sin(gamma)**2*np.cos(gamma)*chi2['yyy'] \
                - 2*np.sin(gamma)*np.cos(gamma)**2*chi2['yxy'],
        'xzz' : + np.sin(gamma)*chi2['xzz'] - np.cos(gamma)*chi2['yzz'],
        'xyz' : + np.sin(gamma)**2*chi2['xyz'] \
                + np.sin(gamma)*np.cos(gamma)*chi2['xxz'] \
                - np.sin(gamma)*np.cos(gamma)*chi2['yyz'] \
                - np.cos(gamma)**2*chi2['yxz'],
        'xxz' : - np.sin(gamma)*np.cos(gamma)*chi2['xyz'] \
                + np.sin(gamma)**2*chi2['xxz'] \
                + np.cos(gamma)**2*chi2['yyz'] \
                - np.sin(gamma)*np.cos(gamma)*chi2['yxz'],
        'xxy' : + np.sin(gamma)**2*np.cos(gamma)*chi2['xxx'] \
                - np.sin(gamma)**2*np.cos(gamma)*chi2['xyy'] \
                + (np.sin(gamma)**3 - np.sin(gamma)*np.cos(gamma)**2)*chi2['xxy'] \
                - np.sin(gamma)*np.cos(gamma)**2*chi2['yxx'] \
                - np.sin(gamma)*np.cos(gamma)**2*chi2['yyy'] \
                + (np.cos(gamma)**3 - np.sin(gamma)**2*np.cos(gamma))*chi2['yxy'],
        'yxx' : + np.sin(gamma)**2*np.cos(gamma)*chi2['xxx'] \
                + np.cos(gamma)**3*chi2['xyy'] \
                - 2*np.sin(gamma)*np.cos(gamma)**2*chi2['xxy'] \
                + np.sin(gamma)**3*chi2['yxx'] \
                + np.sin(gamma)*np.cos(gamma)**2*chi2['yyy'] \
                - 2*np.sin(gamma)**2*np.cos(gamma)*chi2['yxy'],
        'yyy' : + np.cos(gamma)**3*chi2['xxx'] \
                + np.sin(gamma)**2*np.cos(gamma)*chi2['xyy'] \
                + 2*np.sin(gamma)*np.cos(gamma)**2*chi2['xxy'] \
                + np.sin(gamma)*np.cos(gamma)**2*chi2['yxx'] \
                + np.sin(gamma)**3*chi2['yyy'] \
                + 2*np.sin(gamma)**2*np.cos(gamma)*chi2['yxy'],
        'yzz' : + np.cos(gamma)*chi2['xzz'] + np.sin(gamma)*chi2['yzz'],
        'yyz' : + np.sin(gamma)*np.cos(gamma)*chi2['xyz'] \
                + np.cos(gamma)**2*chi2['xxz'] \
                + np.sin(gamma)**2*chi2['yyz'] \
                + np.sin(gamma)*np.cos(gamma)*chi2['yxz'],
        'yxz' : - np.cos(gamma)**2*chi2['xyz'] \
                + np.sin(gamma)*np.cos(gamma)*chi2['xxz'] \
                - np.sin(gamma)*np.cos(gamma)*chi2['yyz'] \
                + np.sin(gamma)**2*chi2['yxz'],
        'yxy' : + np.sin(gamma)*np.cos(gamma)**2*chi2['xxx'] \
                - np.sin(gamma)*np.cos(gamma)**2*chi2['xyy'] \
                - (np.cos(gamma)**3 - np.sin(gamma)**2*np.cos(gamma))*chi2['xxy'] \
                + np.sin(gamma)**2*np.cos(gamma)*chi2['yxx'] \
                - np.sin(gamma)**2*np.cos(gamma)*chi2['yyy'] \
                + (np.sin(gamma)**3 - np.sin(gamma)*np.cos(gamma)**2)*chi2['yxy'],
        'zxx' : + np.sin(gamma)**2*chi2['zxx'] \
                + np.cos(gamma)**2*chi2['zyy'] \
                - 2*np.sin(gamma)*np.cos(gamma)*chi2['zxy'],
        'zyy' : + np.cos(gamma)**2*chi2['zxx'] \
                + np.sin(gamma)**2*chi2['zyy'] \
                + 2*np.sin(gamma)*np.cos(gamma)*chi2['zxy'],
        'zzz' : + chi2['zzz'],
        'zyz' : + np.sin(gamma)*chi2['zyz'] + np.cos(gamma)*chi2['zxz'],
        'zxz' : - np.cos(gamma)*chi2['zyz'] + np.sin(gamma)*chi2['zxz'],
        'zxy' : + np.sin(gamma)*np.cos(gamma)*chi2['zxx'] \
                - np.sin(gamma)*np.cos(gamma)*chi2['zyy'] \
                - np.cos(2*gamma)*chi2['zxy']
    }
    return chi2rot


def wvec(eps, theta):
    '''
    Wave vector, where w = sqrt[epsilon - sin(theta)^2].
    '''
    return np.sqrt(eps - (np.sin(theta)**2))


def frefs(epsi, epsj, theta):
    '''
    Generic reflection fresnel factors, see Eqs. (13) and (14) of PRB 94, 115314 (2016).
    '''
    return (wvec(epsi, theta) - wvec(epsj, theta))/(wvec(epsi, theta) + wvec(epsj, theta))


def frefp(epsi, epsj, theta):
    '''
    Generic reflection fresnel factors, see Eqs. (13) and (14) of PRB 94, 115314 (2016).
    '''
    return ((wvec(epsi, theta) * epsj) - (wvec(epsj, theta) * epsi))/\
           ((wvec(epsi, theta) * epsj) + (wvec(epsj, theta) * epsi))


def ftrans(epsi, epsj, theta):
    '''
    s-polarized transmission fresnel factors , see Eqs. (13) and (14) of PRB 94, 115314 (2016).
    '''
    return (2 * wvec(epsi, theta))/(wvec(epsi, theta) + wvec(epsj, theta))


def ftranp(epsi, epsj, theta):
    '''
    p-polarized transmission fresnel factors, see Eqs. (13) and (14) of PRB 94, 115314 (2016).
    '''
    return (2 * wvec(epsi, theta) * np.sqrt(epsi * epsj))/\
           (wvec(epsi, theta) * epsj + wvec(epsj, theta) * epsi)


def mrc2w(energy, eps0, fres1, fres2, theta, thick, mref):
    '''
    2w multiple reflection coefficient, see Eq. (18) of PRB 94, 115314 (2016).
    '''
    if mref:
        delta = 8*np.pi * ((energy * thick * 1e-9)/(constants.value("Planck constant in eV s") * constants.c)) * wvec(eps0, theta) 
        return (fres1 * np.exp(1j * (delta/2)))/(1 + (fres2 * fres1 * np.exp(1j * delta))) * np.sinc(delta/2)
    elif not mref:
        return fres1


def mrc1w(energy, eps0, fres1, fres2, theta, thick, mref):
    '''
    1w multiple reflection coefficient, see Eq. (21) of PRB 94, 115314 (2016).
    '''
    if mref:
        varphi = 4*np.pi * ((energy * thick * 1e-9)/(constants.value("Planck constant in eV s") * constants.c)) * wvec(eps0, theta)
        return (fres1 * np.exp(1j * varphi))/(1 + (fres2 * fres1 * np.exp(1j * varphi)))
    elif not mref:
        return fres1


def rad_pp(energy, eps1w, eps2w, chi2, theta, phi, thick, mref):
    '''
    rpP, see Eq. (50) of PRB 94, 115314 (2016).
    '''
    fres2p = mrc2w(energy,
                   eps2w['M2'],
                   frefp(eps2w['M2'], eps2w['M3'], theta),
                   frefp(eps2w['M1'], eps2w['M2'], theta),
                   theta, thick, mref)
    fres1p = mrc1w(energy,
                   eps1w['M2'],
                   frefp(eps1w['M2'], eps1w['M3'], theta),
                   frefp(eps1w['M1'], eps1w['M2'], theta),
                   theta, thick, mref)
    pre = (ftranp(eps2w['M1'], eps2w['M2'], theta)/np.sqrt(eps2w['M2'])) * \
          (ftranp(eps1w['M1'], eps1w['M2'], theta)/np.sqrt(eps1w['M2']))**2
    ### r_{pP}
    rpp = - ((1 - fres2p) * (1 - fres1p)**2 \
            * wvec(eps1w['M2'], theta)**2 * wvec(eps2w['M2'], theta) \
            * np.cos(phi)**3 * chi2['xxx']) \
          - (2 * (1 - fres2p) * (1 - fres1p)**2 \
            * wvec(eps1w['M2'], theta)**2 * wvec(eps2w['M2'], theta) \
            * np.sin(phi) * np.cos(phi)**2 * chi2['xxy']) \
          - (2 * (1 - fres2p) * (1 + fres1p) * (1 - fres1p) \
            * wvec(eps1w['M2'], theta) * wvec(eps2w['M2'], theta) \
            * np.sin(theta) * np.cos(phi)**2 * chi2['xxz']) \
          - ((1 - fres2p) * (1 - fres1p)**2 \
            * wvec(eps1w['M2'], theta)**2 * wvec(eps2w['M2'], theta) \
            * np.sin(phi)**2 * np.cos(phi) * chi2['xyy']) \
          - (2 * (1 - fres2p) * (1 + fres1p) * (1 - fres1p) \
            * wvec(eps1w['M2'], theta) * wvec(eps2w['M2'], theta) \
            * np.sin(theta) * np.sin(phi) * np.cos(phi) * chi2['xyz']) \
          - ((1 - fres2p) * (1 + fres1p)**2 \
            * wvec(eps2w['M2'], theta) \
            * np.sin(theta)**2 * np.cos(phi) * chi2['xzz']) \
          - ((1 - fres2p) * (1 - fres1p)**2 \
            * wvec(eps1w['M2'], theta)**2 * wvec(eps2w['M2'], theta) \
            * np.sin(phi) * np.cos(phi)**2 * chi2['yxx']) \
          - (2 * (1 - fres2p) * (1 - fres1p)**2 \
            * wvec(eps1w['M2'], theta)**2 * wvec(eps2w['M2'], theta) \
            * np.sin(phi)**2 * np.cos(phi) * chi2['yxy']) \
          - (2 * (1 - fres2p) * (1 + fres1p) * (1 - fres1p) \
            * wvec(eps1w['M2'], theta) * wvec(eps2w['M2'], theta) \
            * np.sin(theta) * np.sin(phi) * np.cos(phi) * chi2['yxz']) \
          - ((1 - fres2p) * (1 - fres1p)**2 \
            * wvec(eps1w['M2'], theta)**2 * wvec(eps2w['M2'], theta) \
            * np.sin(phi)**3 * chi2['yyy']) \
          - (2 * (1 - fres2p) * (1 + fres1p) * (1 - fres1p) \
            * wvec(eps1w['M2'], theta) * wvec(eps2w['M2'], theta) \
            * np.sin(theta) * np.sin(phi)**2 * chi2['yyz']) \
          - ((1 - fres2p) * (1 + fres1p)**2 \
            * wvec(eps2w['M2'], theta) \
            * np.sin(theta)**2 * np.sin(phi) * chi2['yzz']) \
          + ((1 + fres2p) * (1 - fres1p)**2 \
            * wvec(eps1w['M2'], theta)**2 \
            * np.sin(theta) * np.cos(phi)**2 * chi2['zxx']) \
          + (2 * (1 + fres2p) * (1 + fres1p) * (1 - fres1p) \
            * wvec(eps1w['M2'], theta) \
            * np.sin(theta)**2 * np.cos(phi) * chi2['zxz']) \
          + (2 * (1 + fres2p) * (1 - fres1p)**2 \
            * wvec(eps1w['M2'], theta)**2 \
            * np.sin(theta) * np.sin(phi) * np.cos(phi) * chi2['zxy']) \
          + ((1 + fres2p) * (1 - fres1p)**2 \
            * wvec(eps1w['M2'], theta)**2 \
            * np.sin(theta) * np.sin(phi)**2 * chi2['zyy']) \
          + (2 * (1 + fres2p) * (1 + fres1p) * (1 - fres1p) \
            * wvec(eps1w['M2'], theta) \
            * np.sin(theta)**2 * np.sin(phi) * chi2['zyz']) \
          + ((1 + fres2p) * (1 + fres1p)**2 * np.sin(phi)**3 * chi2['zzz'])
    return pre*rpp


def rad_ps(energy, eps1w, eps2w, chi2, theta, phi, thick, mref):
    '''
    rpS, see Eq. (55) of PRB 94, 115314 (2016).
    '''
    fres2s = mrc2w(energy,
                   eps2w['M2'],
                   frefs(eps2w['M2'], eps2w['M3'], theta),
                   frefs(eps2w['M1'], eps2w['M2'], theta),
                   theta, thick, mref)
    fres1p = mrc1w(energy,
                   eps1w['M2'],
                   frefp(eps1w['M2'], eps1w['M3'], theta),
                   frefp(eps1w['M1'], eps1w['M2'], theta),
                   theta, thick, mref)
    pre = (ftrans(eps2w['M1'], eps2w['M2'], theta) * (1 + fres2s)) * \
          (ftranp(eps1w['M1'], eps1w['M2'], theta)/np.sqrt(eps1w['M2']))**2
    ### r_{pS}
    rps = - ((1 - fres1p)**2 * wvec(eps1w['M2'], theta)**2 \
            * np.sin(phi) * np.cos(phi)**2 * chi2['xxx']) \
          - (2 * (1 - fres1p)**2 * wvec(eps1w['M2'], theta)**2 \
            * np.sin(phi)**2 * np.cos(phi) * chi2['xxy']) \
          - (2 * (1 + fres1p) * (1 - fres1p) * wvec(eps1w['M2'], theta) \
            * np.sin(theta) * np.sin(phi) * np.cos(phi) * chi2['xxz']) \
          - ((1 - fres1p)**2 * wvec(eps1w['M2'], theta)**2 \
            * np.sin(phi)**3 * chi2['xyy']) \
          - (2 * (1 + fres1p) * (1 - fres1p) * wvec(eps1w['M2'], theta) \
            * np.sin(theta) * np.sin(phi)**2 * chi2['xyz']) \
          - ((1 + fres1p)**2 \
            * np.sin(theta)**2 * np.sin(phi) * chi2['xzz']) \
          + ((1 - fres1p)**2 * wvec(eps1w['M2'], theta)**2 \
            * np.cos(phi)**3 * chi2['yxx']) \
          + (2 * (1 - fres1p)**2 * wvec(eps1w['M2'], theta)**2 \
            * np.sin(phi) * np.cos(phi)**2 * chi2['yxy']) \
          + (2 * (1 + fres1p) * (1 - fres1p) * wvec(eps1w['M2'], theta) \
            * np.sin(theta) * np.cos(phi)**2 * chi2['yxz']) \
          + ((1 - fres1p)**2 * wvec(eps1w['M2'], theta)**2 \
            * np.sin(phi)**2 * np.cos(phi) * chi2['yyy']) \
          + (2 * (1 + fres1p) * (1 - fres1p) * wvec(eps1w['M2'], theta) \
            * np.sin(theta) * np.sin(phi) * np.cos(phi) * chi2['yyz']) \
          + ((1 + fres1p)**2 \
            * np.sin(theta)**2 * np.cos(phi) * chi2['yzz'])
    return pre*rps


def rad_sp(energy, eps1w, eps2w, chi2, theta, phi, thick, mref):
    '''
    rsP, see Eq. (60) of PRB 94, 115314 (2016).
    '''
    fres2p = mrc2w(energy,
                   eps2w['M2'],
                   frefp(eps2w['M2'], eps2w['M3'], theta),
                   frefp(eps2w['M1'], eps2w['M2'], theta),
                   theta, thick, mref)
    fres1s = mrc1w(energy,
                   eps1w['M2'],
                   frefs(eps1w['M2'], eps1w['M3'], theta),
                   frefs(eps1w['M1'], eps1w['M2'], theta),
                   theta, thick, mref)
    pre = (ftranp(eps2w['M1'], eps2w['M2'], theta)/np.sqrt(eps2w['M2'])) * \
          (ftrans(eps1w['M1'], eps1w['M2'], theta) * (1 + fres1s))**2
    ### r_{sP}
    rsp = - ((1 - fres2p) * wvec(eps2w['M2'], theta) \
            * np.sin(phi)**2 * np.cos(phi) * chi2['xxx']) \
          + ((1 - fres2p) * wvec(eps2w['M2'], theta) \
            * 2 * np.sin(phi) * np.cos(phi)**2 * chi2['xxy']) \
          - ((1 - fres2p) * wvec(eps2w['M2'], theta) \
            * np.cos(phi)**3 * chi2['xyy']) \
          - ((1 - fres2p) * wvec(eps2w['M2'], theta) \
            * np.sin(phi)**3 * chi2['yxx']) \
          + ((1 - fres2p) * wvec(eps2w['M2'], theta) \
            * 2 * np.sin(phi)**2 * np.cos(phi) * chi2['yxy']) \
          - ((1 - fres2p) * wvec(eps2w['M2'], theta) \
            * np.sin(phi) * np.cos(phi)**2 * chi2['yyy']) \
          + ((1 + fres2p) \
            * np.sin(theta) * np.sin(phi)**2 * chi2['zxx']) \
          - ((1 + fres2p) \
            * np.sin(theta) * 2 * np.sin(phi) * np.cos(phi) * chi2['zxy']) \
          + ((1 + fres2p) \
            * np.sin(theta) * np.cos(phi)**2 * chi2['zyy'])
    return pre*rsp


def rad_ss(energy, eps1w, eps2w, chi2, theta, phi, thick, mref):
    '''
    rsS, see Eq. (65) of PRB 94, 115314 (2016).
    '''
    fres2s = mrc2w(energy,
                   eps2w['M2'],
                   frefs(eps2w['M2'], eps2w['M3'], theta),
                   frefs(eps2w['M1'], eps2w['M2'], theta),
                   theta, thick, mref)
    fres1s = mrc1w(energy,
                   eps1w['M2'],
                   frefs(eps1w['M2'], eps1w['M3'], theta),
                   frefs(eps1w['M1'], eps1w['M2'], theta),
                   theta, thick, mref)
    pre = (ftrans(eps2w['M1'], eps2w['M2'], theta) * (1 + fres2s)) * \
          (ftrans(eps1w['M1'], eps1w['M2'], theta) * (1 + fres1s))**2
    ### r_{sS}
    rss = - (np.sin(phi)**3 * chi2['xxx']) \
          + (2 * np.sin(phi)**2 * np.cos(phi) * chi2['xxy']) \
          - (np.sin(phi) * np.cos(phi)**2 * chi2['xyy']) \
          + (np.sin(phi)**2 * np.cos(phi) * chi2['yxx']) \
          + (np.cos(phi)**3 * chi2['yyy']) \
          - (2 * np.sin(phi) * np.cos(phi)**2 * chi2['yxy'])
    return pre*rss


def shgyield(energy, eps_m1, eps_m2, eps_m3, chi2, theta, phi, thick, gamma=90, sigma_eps=0.0, sigma_chi=0.0, sigma_out=5.0, mode='3-layer', mref=True):
    '''
    Calculates the final SHG yield, see Eq. (38) of PRB 94, 115314 (2016).
    '''

    eps_m1_num = {key: val for key, val in eps_m1.items() if not isinstance(val, dict)}
    eps_m2_num = {key: val for key, val in eps_m2.items() if not isinstance(val, dict)}
    eps_m3_num = {key: val for key, val in eps_m3.items() if not isinstance(val, dict)}
    eps_m1_brd = {key: {'energy': val['energy'], 'data': broadC(val['data'], sigma_eps)} for key, val in eps_m1.items() if isinstance(val, dict)}
    eps_m2_brd = {key: {'energy': val['energy'], 'data': broadC(val['data'], sigma_eps)} for key, val in eps_m2.items() if isinstance(val, dict)}
    eps_m3_brd = {key: {'energy': val['energy'], 'data': broadC(val['data'], sigma_eps)} for key, val in eps_m3.items() if isinstance(val, dict)}
    eps_m1_spl = splineEPS(eps_m1_brd, energy)
    eps_m2_spl = splineEPS(eps_m2_brd, energy)
    eps_m3_spl = splineEPS(eps_m3_brd, energy)
    eps_m1_spl[0].update(eps_m1_num)
    eps_m2_spl[0].update(eps_m2_num)
    eps_m3_spl[0].update(eps_m3_num)
    eps_m1_spl[1].update(eps_m1_num)
    eps_m2_spl[1].update(eps_m2_num)
    eps_m3_spl[1].update(eps_m3_num)

    if mode == '3-layer':
        eps1w = {'M1': avgEPS(eps_m1_spl[0]),
                 'M2': avgEPS(eps_m2_spl[0]),
                 'M3': avgEPS(eps_m3_spl[0])}
        eps2w = {'M1': avgEPS(eps_m1_spl[1]),
                 'M2': avgEPS(eps_m2_spl[1]),
                 'M3': avgEPS(eps_m3_spl[1])}
    elif mode == 'fresnel' or mode == '2-layer':
        ## Fresnel reflection model, see PRB 93, 235304 (2016) and references therein.
        thick = 0
        eps1w = {'M1' : avgEPS(eps_m1_spl[0]),
                 'M2' : avgEPS(eps_m3_spl[0]),
                 'M3' : avgEPS(eps_m3_spl[0])}
        eps2w = {'M1' : avgEPS(eps_m1_spl[1]),
                 'M2' : avgEPS(eps_m1_spl[1]),
                 'M3' : avgEPS(eps_m3_spl[1])}

    chi2_num = {key: val for key, val in chi2.items() if not isinstance(val, dict)}
    chi2_spl = {key: splineC(broadC(val['data'], sigma_chi), val['energy']) for key, val in chi2.items() if isinstance(val, dict)}
    chi2_new = {key: val[0](energy) + 1j*val[1](energy) for key, val in chi2_spl.items()}
    chi2_new.update(chi2_num)
    chi2_rot = rotate(chi2_new, np.radians(gamma))

    m2tocm2 = 1e4 # Convert from m^2 to cm^2
    prefactor = m2tocm2 * ((energy/constants.value("Planck constant over 2 pi in eV s"))**2)/\
                          (2 * constants.epsilon_0 * constants.c**3 * np.cos(np.radians(theta))**2)
    dakkar = {'energy': energy,
              'phi': phi,
              'pp': broad(prefactor * np.absolute(1/np.sqrt(eps1w['M2']) * rad_pp(energy, eps1w, eps2w, chi2_rot, np.radians(theta), np.radians(phi), thick, mref))**2, sigma_out),
              'ps': broad(prefactor * np.absolute(1/np.sqrt(eps1w['M2']) * rad_ps(energy, eps1w, eps2w, chi2_rot, np.radians(theta), np.radians(phi), thick, mref))**2, sigma_out),
              'sp': broad(prefactor * np.absolute(1/np.sqrt(eps1w['M2']) * rad_sp(energy, eps1w, eps2w, chi2_rot, np.radians(theta), np.radians(phi), thick, mref))**2, sigma_out),
              'ss': broad(prefactor * np.absolute(1/np.sqrt(eps1w['M2']) * rad_ss(energy, eps1w, eps2w, chi2_rot, np.radians(theta), np.radians(phi), thick, mref))**2, sigma_out)}
    return dakkar
