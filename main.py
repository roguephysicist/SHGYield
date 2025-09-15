#!/usr/bin/env python

import numpy as np
from PyQt6 import QtGui, QtWidgets

import shgyield.gui as gui
import shgyield.shg as shg


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

MATERIAL = {
    'Si(111)': {
        'prefix': 'silicon',
        'lsuper': 180,
        'lslab': 142.1885172213904,
        'medium 1': {},
        'medium 2': {},
        'medium 3': {}
    }
}

MATERIAL['Si(111)']['medium 1'] = {'eps': {'xx': 1.0, 'yy': 1.0, 'zz': 1.0}}
MATERIAL['Si(111)']['medium 2'] = {
    'eps': {
        'xx' : loadeps('example/chi1-linear/SiH1x1-chi1-xx', MATERIAL['Si(111)']['lsuper']/MATERIAL['Si(111)']['lslab']),
        'yy' : loadeps('example/chi1-linear/SiH1x1-chi1-yy', MATERIAL['Si(111)']['lsuper']/MATERIAL['Si(111)']['lslab']),
        'zz' : loadeps('example/chi1-linear/SiH1x1-chi1-zz', MATERIAL['Si(111)']['lsuper']/MATERIAL['Si(111)']['lslab'])
    },
    'chi2': {
        'xxx': loadshg('example/chi2-nonlinear/SiH1x1-chi2-xxx', 1e6 * 1e-24),
        'xxz': loadshg('example/chi2-nonlinear/SiH1x1-chi2-xxz', 1e6 * 1e-24),
        'zxx': loadshg('example/chi2-nonlinear/SiH1x1-chi2-zxx', 1e6 * 1e-24),
        'zzz': loadshg('example/chi2-nonlinear/SiH1x1-chi2-zzz', 1e6 * 1e-24)
    }
}
MATERIAL['Si(111)']['medium 3'] = {
    'eps': {
        'xx' : loadeps('example/chi1-linear/SiBulk-chi1-xx', 1),
        'yy' : loadeps('example/chi1-linear/SiBulk-chi1-yy', 1),
        'zz' : loadeps('example/chi1-linear/SiBulk-chi1-zz', 1)
    }
}

MATERIAL['Si(111)']['medium 2']['chi2']['xyy'] = {'energy': MATERIAL['Si(111)']['medium 2']['chi2']['xxx']['energy'], 'data': -1*MATERIAL['Si(111)']['medium 2']['chi2']['xxx']['data']}
MATERIAL['Si(111)']['medium 2']['chi2']['yxy'] = {'energy': MATERIAL['Si(111)']['medium 2']['chi2']['xxx']['energy'], 'data': -1*MATERIAL['Si(111)']['medium 2']['chi2']['xxx']['data']}
MATERIAL['Si(111)']['medium 2']['chi2']['yyz'] = MATERIAL['Si(111)']['medium 2']['chi2']['xxz']
MATERIAL['Si(111)']['medium 2']['chi2']['zyy'] = MATERIAL['Si(111)']['medium 2']['chi2']['zxx']
MATERIAL['Si(111)']['medium 2']['chi2']['xzz'] = 0
MATERIAL['Si(111)']['medium 2']['chi2']['xyz'] = 0
MATERIAL['Si(111)']['medium 2']['chi2']['xxy'] = 0
MATERIAL['Si(111)']['medium 2']['chi2']['yxx'] = 0
MATERIAL['Si(111)']['medium 2']['chi2']['yyy'] = 0
MATERIAL['Si(111)']['medium 2']['chi2']['yzz'] = 0
MATERIAL['Si(111)']['medium 2']['chi2']['yxz'] = 0
MATERIAL['Si(111)']['medium 2']['chi2']['zyz'] = 0
MATERIAL['Si(111)']['medium 2']['chi2']['zxz'] = 0
MATERIAL['Si(111)']['medium 2']['chi2']['zxy'] = 0

################################################################################
## Data: Experiment
################################################################################

exp_e, exp_rpp, exp_rsp, exp_rps = np.loadtxt('example/reference/experiment.dat', unpack=True)

EXP = {
    'Si(111)': {
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
            }
        }
    }
}

INIT = {
    'energy': {
        'spect': np.linspace(1.25, 2.50, 126),
        'polar': 1.7
    },
    'theta': 65,
    'phi': 30,
    'gamma': 0,
    'broad': {
        'eps': 0,
        'chi': 0,
        'out': 5
    }
}

if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    widget = gui.CustomWidget(INIT, MATERIAL, EXP, 1)
    widget.show()
    app.exec()

# SPECT = shg.shgyield(energy =    INIT['energy']['spect'],
#                      eps_m1 =    MATERIAL['Si(111)']['medium 1']['eps'],
#                      eps_m2 =    MATERIAL['Si(111)']['medium 2']['eps'],
#                      eps_m3 =    MATERIAL['Si(111)']['medium 3']['eps'],
#                      chi2 =      MATERIAL['Si(111)']['medium 2']['chi2'],
#                      theta =     INIT['theta'],
#                      phi =       INIT['phi'],
#                      gamma =     INIT['gamma'],
#                      thick =     MATERIAL['Si(111)']['lslab']*5.2918E-2,
#                      sigma_eps = INIT['broad']['eps'],
#                      sigma_chi = INIT['broad']['chi'],
#                      sigma_out = INIT['broad']['out'])

# np.savetxt('spect.dat',
#            np.column_stack((SPECT['energy'], SPECT['pp'], SPECT['sp'], SPECT['ps'], SPECT['ss'])),
#            fmt = '%07.4f  %.8e  %.8e  %.8e  %.8e',
#            header = 'w(eV)  RpP             RsP             RpS             RsS')
