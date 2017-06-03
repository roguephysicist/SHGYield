#!/usr/bin/env python
"""
A test script that tests relevant functions and finally compares the numerical
accuracy of the final calculation.
"""

import numpy as np

def read_eps(in_file):
    ''' Reads epsilon file, returns a complex numpy array '''
    data = np.loadtxt(in_file, unpack=True, skiprows=1)
    epsa = data[1] + (1j * data[2]) # complex average
    return epsa

def read_shg(infile):
    ''' Reads chi2 file, returns a complex numpy array of the same length as our energy array '''
    data = np.loadtxt(infile, unpack=True, skiprows=1)
    comp = data[1] + (1j * data[2])
    chi2 = comp[:100]
    return chi2

def test_epsilon():
    ''' Creates and tests the epsilon arrays, then returns dictionary with appropriate values '''
    epsl = read_eps('tests/data/SiH1x1-epsilon') # Epsilon from chi1, layered, normalized
    epsb = read_eps('tests/data/SiBulk-epsilon') # Epsilon from chi1, bulk
    assert epsl.all()
    assert epsb.all()
    eps = {'v1w': 1, 'v2w': 1,
           'b1w': epsb[:100], 'b2w': epsb[1::2][:100],
           'l1w': epsl[:100], 'l2w': epsl[1::2][:100]}
    return eps

def test_chi2():
    ''' Creates and tests the chi2 array, then returns dictionary with appropriate components '''
    xxx = read_shg('tests/data/SiH1x1-chi2-xxx')
    assert xxx.all()
    chi2 = {'xxx': xxx, 'xyy': -1 * xxx, 'yxy': -1 * xxx}
    return chi2

def test_wvecs():
    ''' Creates a dictionary with the wave vectors, tests every value, then returns it '''
    vecs = {'v1w': np.sqrt(EPS['v1w'] - (np.sin(1.1344640137963142)**2)),
            'v2w': np.sqrt(EPS['v2w'] - (np.sin(1.1344640137963142)**2)),
            'b1w': np.sqrt(EPS['b1w'] - (np.sin(1.1344640137963142)**2)),
            'b2w': np.sqrt(EPS['b2w'] - (np.sin(1.1344640137963142)**2)),
            'l1w': np.sqrt(EPS['l1w'] - (np.sin(1.1344640137963142)**2)),
            'l2w': np.sqrt(EPS['l2w'] - (np.sin(1.1344640137963142)**2))}
    for value in vecs.values():
        assert value.all()
    return vecs

def freflc(pol, veci, vecj, epsi, epsj):
    ''' Reflection Fresnel factors '''
    if pol == "p":
        factor = ((veci * epsj) - (vecj * epsi))/((veci * epsj) + (vecj * epsi))
    elif pol == "s":
        factor = (veci - vecj)/(veci + vecj)
    return factor

def ftrans(pol, veci, vecj, epsi, epsj):
    ''' Transmission Fresnel factors '''
    if pol == "p":
        factor = (2 * veci * np.sqrt(epsi * epsj))/(veci * epsj + vecj * epsi)
    elif pol == "s":
        factor = (2 * veci)/(veci + vecj)
    return factor

def test_fresnel():
    ''' Creates a dictionary with the fresnel factors, tests every value, then returns it '''
    fresn = {'tvls': ftrans("s", WVECS['v1w'], WVECS['l1w'], EPS['v1w'], EPS['l1w']),
             'Tvls': ftrans("s", WVECS['v2w'], WVECS['l2w'], EPS['v2w'], EPS['l2w']),
             'rvls': freflc("s", WVECS['v1w'], WVECS['l1w'], EPS['v1w'], EPS['l1w']),
             'rlbs': freflc("s", WVECS['l1w'], WVECS['b1w'], EPS['l1w'], EPS['b1w']),
             'Rvls': freflc("s", WVECS['v2w'], WVECS['l2w'], EPS['v2w'], EPS['l2w']),
             'Rlbs': freflc("s", WVECS['l2w'], WVECS['b2w'], EPS['l2w'], EPS['b2w'])}
    for value in fresn.values():
        assert value.all()
    return fresn

def test_multiref():
    '''
    Calculates the different elements for the multiple reflections framework,
    then creates a dictionary and tests each value, and returns it
    '''
    varphi = 4 * np.pi * ((ONEE * 10 * 1e-9)/(4.135667662e-15 * 299792458.0)) * WVECS['l1w']
    delta = 8 * np.pi * ((ONEE * 10 * 1e-9)/(4.135667662e-15 * 299792458.0)) * WVECS['l2w']
    multi = {'s2w': (1 + ((FRESN['Rlbs'] * np.exp(1j * delta/2))/(1 + (FRESN['Rvls'] * FRESN['Rlbs'] * np.exp(1j * delta))) * np.sin(delta/2)/(delta/2))),
             's1w': (1 + ((FRESN['rlbs'] * np.exp(1j * varphi))/(1 + (FRESN['rvls'] * FRESN['rlbs'] * np.exp(1j * varphi)))))}
    for value in multi.values():
        assert value.all()
    return multi

def test_shgyield():
    '''
    Calculates the reflectance, multiplies everything out, then tests against the
    reference data to verify numerical accuracy
    '''
    multiref = test_multiref()
    chi2 = test_chi2()
    gamma = FRESN['Tvls'] * multiref['s2w'] * (FRESN['tvls'] * multiref['s1w'])**2 # s-in, S-out
    factor = - (np.sin(0.52359877559829882)**3 * chi2['xxx']) \
             - (np.sin(0.52359877559829882) * np.cos(0.52359877559829882)**2 * chi2['xyy']) \
             - (2 * np.sin(0.52359877559829882) * np.cos(0.52359877559829882)**2 * chi2['yxy'])
    prefactor = (2 * 8.854187817620389e-12 * 6.582119514e-16**2 * \
                 299792458.0**3 * np.cos(1.1344640137963142)**2)**-1
    rss = 1e20 * 1e4 * prefactor * (ONEE ** 2) * \
          np.absolute((1/np.sqrt(EPS['l1w'])) * gamma * factor)**2
    np.testing.assert_allclose(rss, REFERENCE, rtol=1e-06, atol=1e-12)


ONEE = np.linspace(0.01, float(100)/100, 100) # 1w energy array, 0.01-10 eV
EPS = test_epsilon() # Creates the dictionary with the epsilons
WVECS = test_wvecs() # Creates the dictionary with the wave vectors
FRESN = test_fresnel() # Creates the dictionary with the fresnel factors
REFERENCE = np.array( # An array with valid reference data, for numerical comparison
    [1.032403505838488195e-11, 4.550832431431746097e-11, 1.325273386471492838e-10,
     3.333508163754209330e-10, 7.554747698066963989e-10, 1.556695234237770624e-09,
     2.941988999332841802e-09, 5.157168742134966043e-09, 8.484673625861601399e-09,
     1.324308967483312284e-08, 1.979021478650005457e-08, 2.852345885277600976e-08,
     3.988427708072623743e-08, 5.435800004817903165e-08, 7.247274260381388879e-08,
     9.480304677776567662e-08, 1.219680461564326177e-07, 1.546347511097495607e-07,
     1.935193991492270539e-07, 2.393815220940141246e-07, 2.930379257226721779e-07,
     3.553525295537086449e-07, 4.272414482187805282e-07, 5.096801600445347622e-07,
     6.036935338622235096e-07, 7.103737319066513703e-07, 8.308602252547417564e-07,
     9.663698498429709728e-07, 1.118171168385683709e-06, 1.287595691327597088e-06,
     1.476063200502429092e-06, 1.685052877003324629e-06, 1.916095668833173497e-06,
     2.170856280534693416e-06, 2.451008203944823028e-06, 2.758378163507136639e-06,
     3.094825349682943973e-06, 3.462323803082928682e-06, 3.862946150516640773e-06,
     4.298846852052466438e-06, 4.772285811578160461e-06, 5.285649321952998811e-06,
     5.841420303006999191e-06, 6.442204926764959899e-06, 7.090755403550167435e-06,
     7.789894907397299111e-06, 8.542631623870533243e-06, 9.352098674582466289e-06,
     1.022155410507798585e-05, 1.115449176750408030e-05, 1.215448941129900609e-05,
     1.322534376754085624e-05, 1.437106807222101180e-05, 1.559573077643159491e-05,
     1.690374158915691161e-05, 1.829969652571109294e-05, 1.978835559189621050e-05,
     2.137487441414569035e-05, 2.306441983598689317e-05, 2.486263870201632124e-05,
     2.677539099976400971e-05, 2.880877329701041957e-05, 3.096933230165414434e-05,
     3.326399829620583587e-05, 3.569962855327095321e-05, 3.828395323839731406e-05,
     4.102498470762885158e-05, 4.393083527062365072e-05, 4.701161749951630371e-05,
     5.027419661173465138e-05, 5.373103924320059218e-05, 5.739117364624867938e-05,
     6.126615955857774140e-05, 6.536727771894309744e-05, 6.970718660051250948e-05,
     7.429824005237595263e-05, 7.915562025273016871e-05, 8.429196764942422141e-05,
     8.972587406413663324e-05, 9.547232530757746244e-05, 1.015486229194248068e-04,
     1.079739761580187610e-04, 1.147685360342658666e-04, 1.219533520956576591e-04,
     1.295520039853262383e-04, 1.375889533546590629e-04, 1.460889787483928806e-04,
     1.550813912972550431e-04, 1.645933039979432340e-04, 1.746610020884219920e-04,
     1.853160275468189509e-04, 1.965975648356654512e-04, 2.085467839250398914e-04,
     2.212026799999379642e-04, 2.346155505875994411e-04, 2.488391475834370015e-04,
     2.639239115208671417e-04, 2.799332689218986974e-04, 2.969349401098667443e-04,
     3.150020689590089578e-04])
