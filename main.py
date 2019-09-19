#!/usr/bin/env python

'''
todo:
* scissors?
* SHG: SOME Nv=1 INSTANCES ARE HARDCODED, NEED TO GO BACK AND CHANGE
* Allow for saving all data and parameters to NetCDF file, and final data to txt
* Develop SHG functions into class
* Improve rotation function and avoid running unless changed
* Improve spline function and avoid running unless changed
* Convert to absolute broadening to avoid trouble with polar plots

* Develop GUI to ingest and pre-process data, provide initial values, etc.
* Include every symmetry group (see Popov) into menus
'''

import numpy as np
from PyQt5 import QtGui

import shgyield.gui as gui


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
## your material!
################################################################################

ESPECT = np.linspace(1.25, 2.50, 126)
EPOLAR = 1.7

# Epsilon, some normalization stuff
LSUPER = 180
LSLAB = 142.1885172213904
CHI1NORM = LSUPER/LSLAB # Normalization for layered chi1

## Chi2: UNITS MUST BE CONVERTED TO m^2/V!!
SCALE = 1e6     # Will multiply your chi2, (your program includes some scaling factor) (DEPRECATE!!)
PM2TOM2 = 1e-24 # Convert from pm^2 to m^2
CHI2NORM = 1    # CHI2 normalization; for bulk: LSUPER * 52.9177210903 for pm^2/V; for surface: 1

MATERIAL = {
    'Si(111)': {
        'thickness': 10,
        'eps': {
            'm1': {
                'xx': 1.0,
                'yy': 1.0,
                'zz': 1.0
            },
            'm2': {
                'xx' : loadeps('example/chi1-linear/SiH1x1-chi1-xx', CHI1NORM),
                'yy' : loadeps('example/chi1-linear/SiH1x1-chi1-yy', CHI1NORM),
                'zz' : loadeps('example/chi1-linear/SiH1x1-chi1-zz', CHI1NORM)
            },
            'm3': {
                'xx' : loadeps('example/chi1-linear/SiBulk-chi1-xx', 1),
                'yy' : loadeps('example/chi1-linear/SiBulk-chi1-yy', 1),
                'zz' : loadeps('example/chi1-linear/SiBulk-chi1-zz', 1)
            }
        },
        'chi2': {
            'xxx': loadshg('example/chi2-nonlinear/SiH1x1-chi2-xxx', CHI2NORM*SCALE*PM2TOM2),
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
            'zxy': 0
        }
    },
}

MATERIAL['Si(111)']['chi2']['xyy'] = {'energy': MATERIAL['Si(111)']['chi2']['xxx']['energy'], 'data': -MATERIAL['Si(111)']['chi2']['xxx']['data']}
MATERIAL['Si(111)']['chi2']['yxy'] = {'energy': MATERIAL['Si(111)']['chi2']['xxx']['energy'], 'data': -MATERIAL['Si(111)']['chi2']['xxx']['data']}
MATERIAL['Si(111)']['chi2']['yyz'] = MATERIAL['Si(111)']['chi2']['xxz']
MATERIAL['Si(111)']['chi2']['zyy'] = MATERIAL['Si(111)']['chi2']['zxx']

################################################################################
## Data: Experiment
################################################################################

exp_e, exp_rpp, exp_rsp, exp_rps = np.loadtxt('example/reference/experiment.dat', unpack=True)

EXP = {
    'Si(111)': {
        # 'polar': {
        #     'pp': {
        #         'energy': ,
        #         'phi': ,
        #         'data':
        #     },
        #     'sp': {
        #         'energy': ,
        #         'phi': ,
        #         'data':
        #     },
        #     'ps': {
        #         'energy': ,
        #         'phi': ,
        #         'data':
        #     },
        #     'ss': {
        #         'energy': ,
        #         'phi': ,
        #         'data':
        #     }
        # },
        'spect': {
            'pp': {
                'energy': exp_e,
                'phi': 30,
                'data': exp_rpp
            },
            'sp': {
                'energy': exp_e,
                'phi': 30,
                'data': exp_rsp
            },
            'ps': {
                'energy': exp_e,
                'phi': 30,
                'data': exp_rps
            },
            # 'ss': {
            #     'energy': exp_e,
            #     'phi': 30,
            #     'data': exp_rss
            # }
        }
    }
}

if __name__ == '__main__':
    app = QtGui.QApplication([])
    widget = gui.CustomWidget(ESPECT, EPOLAR, EXP, MATERIAL, 1e20)
    widget.show()
    app.exec_()
