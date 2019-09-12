#!/usr/bin/env python

'''
todo:
* Create polar grid or radial tics
* Labels and axis labels
* Include every symmetry group (see Popov) into menus
* Develop GUI to ingest and pre-process data, provide initial values, etc.
* Allow for saving all data and parameters to NetCDF file
* Develop SHG functions into class for ease of use
* Improve rotation function and avoid running unless changed
* Improve spline function and avoid running unless changed
* Convert to absolute broadening to avoid trouble with polar plots
* Scissors!
'''

import numpy as np

import shgyield.shg as shg
import shgyield.PlotGUI.gui as gui


################################################################################
## Functions to help ingest and pre-process your data
################################################################################

def loadeps(infile, scale):
    ''' loads chi1 from file, converts to epsilon '''
    energy, real, imag = np.loadtxt(infile, unpack=True)
    return {'energy': energy, 'data': 1 + (4 * np.pi * ((real + 1j*imag) * scale))}

def loadshg(infile, scale):
    ''' loads chi2 from file, scales and converts to appropriate units '''
    energy, re1w, im1w, re2w, im2w = np.loadtxt(infile, unpack=True)
    return {'energy': energy, 'data': scale * ((re1w + re2w) + 1j*(im1w + im2w))}


################################################################################
## Data: Epsilon
################################################################################

LSUPER = 180
LSLAB = 142.1885172213904
CHI1NORM = LSUPER/LSLAB # Normalization for layered chi1

EPS_M1 = {'xx' : 1.0,
          'yy' : 1.0,
          'zz' : 1.0}

EPS_M2 = {'xx' : loadeps('example/chi1-linear/SiH1x1-chi1-xx', CHI1NORM),
          'yy' : loadeps('example/chi1-linear/SiH1x1-chi1-yy', CHI1NORM),
          'zz' : loadeps('example/chi1-linear/SiH1x1-chi1-zz', CHI1NORM)}

EPS_M3 = {'xx' : loadeps('example/chi1-linear/SiBulk-chi1-xx', 1),
          'yy' : loadeps('example/chi1-linear/SiBulk-chi1-yy', 1),
          'zz' : loadeps('example/chi1-linear/SiBulk-chi1-zz', 1)}


################################################################################
## Data: Nonlinear susceptibility
################################################################################

## UNITS MUST BE CONVERTED TO m^2/V!!
SCALE = 1e6     # Will multiply your chi2, (your program includes some scaling factor) (DEPRECATE!!)
PM2TOM2 = 1e-24 # Convert from pm^2 to m^2
CHI2NORM = 1    # CHI2 normalization; for bulk: LSUPER * 52.9177210903 for pm^2/V; for surface: 1

CHI2 = {'xxx': loadshg('example/chi2-nonlinear/SiH1x1-chi2-xxx', CHI2NORM*SCALE*PM2TOM2),
        'xxz': loadshg('example/chi2-nonlinear/SiH1x1-chi2-xxz', CHI2NORM*SCALE*PM2TOM2),
        'zxx': loadshg('example/chi2-nonlinear/SiH1x1-chi2-zxx', CHI2NORM*SCALE*PM2TOM2),
        'zzz': loadshg('example/chi2-nonlinear/SiH1x1-chi2-zzz', CHI2NORM*SCALE*PM2TOM2),
        'xzz': 0,
        'xyz': 0,
        'xxy': 0,
        'yxx': 0,
        'yyy': 0,
        'yzz': 0,
        'yxz': 0,
        'zyz': 0,
        'zxz': 0,
        'zxy': 0}

CHI2['xyy'] = {'energy': CHI2['xxx']['energy'], 'data': -CHI2['xxx']['data']}
CHI2['yxy'] = {'energy': CHI2['xxx']['energy'], 'data': -CHI2['xxx']['data']}
CHI2['yyz'] = CHI2['xxz']
CHI2['zyy'] = CHI2['zxx']


################################################################################
## Data: Experiment
################################################################################

exp_e, exp_rpp, exp_rsp, exp_rps = np.loadtxt('example/reference/experiment.dat', unpack=True)

EXP = {'energy': exp_e,
       'pp': exp_rpp,
       'sp': exp_rsp,
       'ps': exp_rps,
       'ss': np.zeros(exp_e.size)}
       # 'ss': np.full(exp_e.size, np.nan)}


################################################################################
## Calculation of the SHG Yield (math time!)
################################################################################

THICKNESS = 10   # thickness "d" of the thin layer, in nm

THETA0 = 65 # angle of incidence in degrees
PHI = 30    # azimuthal angle in degrees
GAMMA = 90  # rotation for CHI2 in degrees; 90 is equal to no rotation.

SIGMA_EPS = 0   # stdev for eps (input) gaussian broadening
SIGMA_CHI = 0   # stdev for chi2 (input) gaussian broadening
SIGMA_OUT = 5   # stdev for output gaussian broadening

if __name__ == '__main__':
    app = gui.QtGui.QApplication([])
    widget = gui.CustomWidget(1.7, (1.25, 2.5), EPS_M1, EPS_M2, EPS_M3, CHI2, THETA0, PHI, GAMMA, SIGMA_EPS, SIGMA_CHI, SIGMA_OUT, EXP)
    widget.show()
    app.exec_()

# ENERGY = np.linspace(1.25, 2.50, 126)   # Establishes energy range over which the
# YIELD = shg.shgyield(energy = ENERGY,
#                      eps_m1 = EPS_M1,
#                      eps_m2 = EPS_M2,
#                      eps_m3 = EPS_M3,
#                      chi2 = CHI2,
#                      theta = THETA0,
#                      phi = PHI,
#                      thick = THICKNESS,
#                      gamma = GAMMA,
#                      sigma_eps = SIGMA_EPS,
#                      sigma_chi = SIGMA_CHI,
#                      sigma_out = SIGMA_OUT)
#
# DATA = np.column_stack((YIELD['energy'],
#                         YIELD['pp'],
#                         YIELD['sp'],
#                         YIELD['ps'],
#                         YIELD['ss']))
# OUTFILE = 'output.dat'
# np.savetxt(OUTFILE, DATA, fmt=('%07.4f', '%.8e', '%.8e', '%.8e', '%.8e'), delimiter='    ',
#            header='RiF (cm^2/W)\n1w(eV)   RpP'+15*" "+'RsP'+15*" "+'RpS'+15*" "+'RsS')
