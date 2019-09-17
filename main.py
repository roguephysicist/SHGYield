#!/usr/bin/env python

'''
todo:
* scissors?

* Allow for saving all data and parameters to NetCDF file, and final data to txt
* Develop SHG functions into class
* Improve rotation function and avoid running unless changed
* Improve spline function and avoid running unless changed
* Convert to absolute broadening to avoid trouble with polar plots

* Develop GUI to ingest and pre-process data, provide initial values, etc.
* Include every symmetry group (see Popov) into menus
'''

import numpy as np
import pyqtgraph as pg
from PyQt5 import QtGui, QtCore

import shgyield.shg as shg
from shgyield.QtLayout import Ui_CustomWidget

pg.setConfigOption('background', '#EAEAEA')
pg.setConfigOption('foreground', 'k')
pg.setConfigOption('antialias', True)
pg.setConfigOption('useWeave', True)


class CustomWidget(QtGui.QWidget):

    def __init__(self, experiment, parent=None):
        super(CustomWidget, self).__init__(parent=parent)

        # set up the form class as a `ui` attribute
        self.ui = Ui_CustomWidget()
        self.ui.setupUi(self)

        self.scale = 1e20

        self.colors = ['k', '#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']

        self.marker = {
            'pen': None,
            'symbolSize': 12,
            'symbol': 'o',
            'symbolPen': pg.mkPen(None),
            'symbolBrush': pg.mkBrush(204,0,0,150)
        }

        self.widgets = {
            'tab_polar': {
                'pp': self.ui.plot_all_polar_rpp.getPlotItem(),
                'sp': self.ui.plot_all_polar_rsp.getPlotItem(),
                'ps': self.ui.plot_all_polar_rps.getPlotItem(),
                'ss': self.ui.plot_all_polar_rss.getPlotItem()
            },
            'tab_spect': {
                'pp': self.ui.plot_all_spect_rpp.getPlotItem(),
                'sp': self.ui.plot_all_spect_rsp.getPlotItem(),
                'ps': self.ui.plot_all_spect_rps.getPlotItem(),
                'ss': self.ui.plot_all_spect_rss.getPlotItem()
            },
            'tab_pp': {
                'exp_polar': self.ui.plot_rpp_polar_exp.getPlotItem(),
                'exp_spect': self.ui.plot_rpp_spect_exp.getPlotItem(),
                'thr_polar': self.ui.plot_rpp_polar.getPlotItem(),
                'thr_spect': self.ui.plot_rpp_spect.getPlotItem()
            },
            'tab_sp': {
                'exp_polar': self.ui.plot_rsp_polar_exp.getPlotItem(),
                'exp_spect': self.ui.plot_rsp_spect_exp.getPlotItem(),
                'thr_polar': self.ui.plot_rsp_polar.getPlotItem(),
                'thr_spect': self.ui.plot_rsp_spect.getPlotItem()
            },
            'tab_ps': {
                'exp_polar': self.ui.plot_rps_polar_exp.getPlotItem(),
                'exp_spect': self.ui.plot_rps_spect_exp.getPlotItem(),
                'thr_polar': self.ui.plot_rps_polar.getPlotItem(),
                'thr_spect': self.ui.plot_rps_spect.getPlotItem()
            },
            'tab_ss': {
                'exp_polar': self.ui.plot_rss_polar_exp.getPlotItem(),
                'exp_spect': self.ui.plot_rss_spect_exp.getPlotItem(),
                'thr_polar': self.ui.plot_rss_polar.getPlotItem(),
                'thr_spect': self.ui.plot_rss_spect.getPlotItem()
            }
        }

        # initializes legends, polar lines
        for pol in ['pp', 'sp', 'ps', 'ss']:
            self.widgets['tab_polar'][pol].addLegend()
            self.widgets['tab_spect'][pol].addLegend()
            self.widgets['tab_' + pol]['thr_polar'].addLegend()
            self.widgets['tab_' + pol]['thr_spect'].addLegend()
            self.widgets['tab_' + pol]['exp_polar'].addLegend()
            self.widgets['tab_' + pol]['exp_spect'].addLegend()
            for ang in np.arange(0, 361, 30):
                self.widgets['tab_polar'][pol].addItem(
                    pg.InfiniteLine(angle = ang, pen = pg.mkPen(color=0.8, width=0.5, style=QtCore.Qt.DashLine)))

        # calculation time!
        # init: energy
        self.einc = 0.01
        self.epolar = 1.7
        self.erng = np.linspace(1.25, 2.50, 126)
        # self.erng = np.arange(self.ui.box_energy_spect_min.value(),
        #                       self.ui.box_energy_spect_max.value()+self.einc,
        #                       self.einc)
        self.ui.box_energy_polar.setValue(self.epolar)
        self.ui.box_energy_polar.setMinimum(self.erng.min())
        self.ui.box_energy_polar.setMaximum(self.erng.max())
        self.ui.box_energy_polar.setSingleStep(self.einc)

        # plots for experimental data, that do not change or need calculations
        deco = 0
        for case, measurement in experiment.items():
            for kind, values in measurement.items():
                for polarization, data in values.items():
                    if kind == 'polar':
                        self.widgets['tab_' + polarization]['exp_' + kind].plot(
                            transform_polar(data['energy'], data['data']*self.scale),
                            pen = None,
                            symbol = deco,
                            symbolSize = 5,
                            symbolBrush = pg.mkBrush(color=self.colors[deco]),
                            symbolPen = pg.mkPen('k', width=0.5),
                            name = case)
                    elif kind == 'spect':
                        if 'stdev' in data.keys():
                            self.widgets['tab_' + polarization]['exp_' + kind].addItem(pg.ErrorBarItem(
                                x = data['energy'],
                                y = data['data']*self.scale,
                                height = data['stdev'],
                                pen = pg.mkPen(color=self.colors[deco], width=1)))
                        self.widgets['tab_' + polarization]['exp_' + kind].plot(
                            x = data['energy'],
                            y = data['data']*self.scale,
                            pen = None,
                            symbol = deco,
                            symbolSize = 5,
                            symbolBrush = pg.mkBrush(color=self.colors[deco]),
                            symbolPen = pg.mkPen('k', width=0.5),
                            name = case)
            deco += 1

        deco = 0
        self.test = {}
        for case in MATERIAL.keys():
            self.polar = shg.shgyield(energy =    self.ui.box_energy_polar.value(),
                                      eps_m1 =    MATERIAL[case]['eps']['m1'],
                                      eps_m2 =    MATERIAL[case]['eps']['m2'],
                                      eps_m3 =    MATERIAL[case]['eps']['m3'],
                                      chi2 =      MATERIAL[case]['chi2'],
                                      theta =     self.ui.sld_angle_theta.value(),
                                      phi =       np.arange(0, 361),
                                      gamma =     self.ui.sld_angle_gamma.value(),
                                      thick =     MATERIAL[case]['thickness'],
                                      sigma_eps = self.ui.sld_broad_eps.value(),
                                      sigma_chi = self.ui.sld_broad_chi.value(),
                                      sigma_out = 0)

            self.spect = shg.shgyield(energy =    self.erng,
                                      eps_m1 =    MATERIAL[case]['eps']['m1'],
                                      eps_m2 =    MATERIAL[case]['eps']['m2'],
                                      eps_m3 =    MATERIAL[case]['eps']['m3'],
                                      chi2 =      MATERIAL[case]['chi2'],
                                      theta =     self.ui.sld_angle_theta.value(),
                                      phi =       self.ui.sld_angle_phi.value(),
                                      gamma =     self.ui.sld_angle_gamma.value(),
                                      thick =     MATERIAL[case]['thickness'],
                                      sigma_eps = self.ui.sld_broad_eps.value(),
                                      sigma_chi = self.ui.sld_broad_chi.value(),
                                      sigma_out = self.ui.sld_broad_out.value())

            self.pidx = np.where(np.isclose(self.spect['phi'], self.polar['phi']))
            self.eidx = np.where(np.isclose(self.spect['energy'], self.polar['energy']))

            self.test[case] = {
                pol: {
                    'test1': self.widgets['tab_polar'][pol].plot(self.trans_polar(self.polar['phi'], self.polar[pol]*self.scale), pen=pg.mkPen(color=self.colors[deco], width=1.5), name=case),
                    'test2': self.widgets['tab_polar'][pol].plot(self.trans_polar2(self.spect['phi'], self.polar[pol][self.pidx][0]*self.scale), **self.marker),
                    'test3': self.widgets['tab_spect'][pol].plot(x = self.spect['energy'], y = self.spect[pol]*self.scale, pen=pg.mkPen(color=self.colors[deco], width=1.5), name=case),
                    'test4': self.widgets['tab_spect'][pol].plot(x = [self.polar['energy']], y = [self.spect[pol][self.eidx][0]*self.scale], **self.marker),
                    'test5': self.widgets['tab_' + pol]['thr_polar'].plot(self.trans_polar(self.polar['phi'], self.polar[pol]*self.scale), pen=pg.mkPen(color=self.colors[deco], width=1.5), name=case),
                    'test6': self.widgets['tab_' + pol]['thr_polar'].plot(self.trans_polar2(self.spect['phi'], self.polar[pol][self.pidx][0]*self.scale), **self.marker),
                    'test7': self.widgets['tab_' + pol]['thr_spect'].plot(x = self.spect['energy'], y = self.spect[pol]*self.scale, pen=pg.mkPen(color=self.colors[deco], width=1.5), name=case),
                    'test8': self.widgets['tab_' + pol]['thr_spect'].plot(x = [self.polar['energy']], y = [self.spect[pol][self.eidx][0]*self.scale], **self.marker)
                } for pol in ['pp', 'sp', 'ps', 'ss']
            } 
            
            deco += 1

        # decorations and behavior for plots        
        for pol in ['pp', 'sp', 'ps', 'ss']:
            self.widgets['tab_polar'][pol].setAspectLocked(True)
            self.widgets['tab_polar'][pol].disableAutoRange()
            self.widgets['tab_polar'][pol].setLabels(bottom='R (10<sup>-20</sup> cm<sup>2</sup>/W)', left='R (10<sup>-20</sup> cm<sup>2</sup>/W)')
            self.widgets['tab_spect'][pol].setLabels(bottom='Photon Energy (eV)', left='R (10<sup>-20</sup> cm<sup>2</sup>/W)')
            self.widgets['tab_spect'][pol].setXLink(self.widgets['tab_spect']['pp'])
            self.widgets['tab_' + pol]['thr_polar'].setLabels(bottom='R (10<sup>-20</sup> cm<sup>2</sup>/W)', left='R (10<sup>-20</sup> cm<sup>2</sup>/W)')
            self.widgets['tab_' + pol]['thr_polar'].setAspectLocked(True)
            self.widgets['tab_' + pol]['thr_polar'].disableAutoRange()
            self.widgets['tab_' + pol]['thr_spect'].setLabels(bottom='Photon Energy (eV)', left='R (10<sup>-20</sup> cm<sup>2</sup>/W)')
            self.widgets['tab_' + pol]['exp_spect'].setXLink(self.widgets['tab_' + pol]['thr_spect'])
            self.widgets['tab_' + pol]['exp_spect'].setLabels(bottom='Photon Energy (eV)')

        # updating
        self.ui.box_energy_polar.valueChanged.connect(self.update_plot)
        # self.ui.box_energy_spect_min.valueChanged.connect(self.update_plot)
        # self.ui.box_energy_spect_max.valueChanged.connect(self.update_plot)
        
        self.ui.sld_angle_theta.valueChanged.connect(self.update_plot)
        self.ui.sld_angle_phi.valueChanged.connect(self.update_plot)
        self.ui.sld_angle_gamma.valueChanged.connect(self.update_plot)
        
        self.ui.sld_broad_eps.valueChanged.connect(self.update_plot)
        self.ui.sld_broad_chi.valueChanged.connect(self.update_plot)
        self.ui.sld_broad_out.valueChanged.connect(self.update_plot)

        # simple demonstration of pure Qt widgets interacting with pyqtgraph
        # self.ui.checkBox.stateChanged.connect(self.toggleMouse)

    # def toggleMouse(self, state):
    #     if state == QtCore.Qt.Checked:
    #         enabled = True
    #     else:
    #         enabled = False

    #     self.ui.plotWidget.setMouseEnabled(x=enabled, y=enabled)
    def trans_polar(self, angle, radius):
        return {'x': radius*np.cos(np.radians(angle)), 'y': radius*np.sin(np.radians(angle))}

    def trans_polar2(self, angle, radius):
        return {'x': [radius*np.cos(np.radians(angle))], 'y': [radius*np.sin(np.radians(angle))]}

    def update_plot(self):
        for case in MATERIAL.keys():
            polar = shg.shgyield(energy =    self.ui.box_energy_polar.value(),
                                 eps_m1 =    MATERIAL[case]['eps']['m1'],
                                 eps_m2 =    MATERIAL[case]['eps']['m2'],
                                 eps_m3 =    MATERIAL[case]['eps']['m3'],
                                 chi2 =      MATERIAL[case]['chi2'],
                                 theta =     self.ui.sld_angle_theta.value(),
                                 phi =       np.arange(0, 361),
                                 gamma =     self.ui.sld_angle_gamma.value(),
                                 thick =     MATERIAL[case]['thickness'],
                                 sigma_eps = self.ui.sld_broad_eps.value(),
                                 sigma_chi = self.ui.sld_broad_chi.value(),
                                 sigma_out = 0)

            spect = shg.shgyield(energy =    self.erng,
                                 eps_m1 =    MATERIAL[case]['eps']['m1'],
                                 eps_m2 =    MATERIAL[case]['eps']['m2'],
                                 eps_m3 =    MATERIAL[case]['eps']['m3'],
                                 chi2 =      MATERIAL[case]['chi2'],
                                 theta =     self.ui.sld_angle_theta.value(),
                                 phi =       self.ui.sld_angle_phi.value(),
                                 gamma =     self.ui.sld_angle_gamma.value(),
                                 thick =     MATERIAL[case]['thickness'],
                                 sigma_eps = self.ui.sld_broad_eps.value(),
                                 sigma_chi = self.ui.sld_broad_chi.value(),
                                 sigma_out = self.ui.sld_broad_out.value())

            pidx = np.where(np.isclose(spect['phi'], polar['phi']))
            eidx = np.where(np.isclose(spect['energy'], polar['energy']))

            for pol in ['pp', 'sp', 'ps', 'ss']:
                self.test[case][pol]['test1'].setData(self.trans_polar(polar['phi'], polar[pol]*self.scale)),
                self.test[case][pol]['test2'].setData(self.trans_polar2(spect['phi'], polar[pol][pidx][0]*self.scale)),
                self.test[case][pol]['test3'].setData(x = spect['energy'], y = spect[pol]*self.scale),
                self.test[case][pol]['test4'].setData(x = [polar['energy']], y = [spect[pol][eidx][0]*self.scale]),
                self.test[case][pol]['test5'].setData(self.trans_polar(polar['phi'], polar[pol]*self.scale)),
                self.test[case][pol]['test6'].setData(self.trans_polar2(spect['phi'], polar[pol][pidx][0]*self.scale)),
                self.test[case][pol]['test7'].setData(x = spect['energy'], y = spect[pol]*self.scale),
                self.test[case][pol]['test8'].setData(x = [polar['energy']], y = [spect[pol][eidx][0]*self.scale])

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

# Epsilon, some normalization stuff
LSUPER = 180
LSLAB = 142.1885172213904
CHI1NORM = LSUPER/LSLAB # Normalization for layered chi1

## Chi2: UNITS MUST BE CONVERTED TO m^2/V!!
SCALE = 1e6     # Will multiply your chi2, (your program includes some scaling factor) (DEPRECATE!!)
PM2TOM2 = 1e-24 # Convert from pm^2 to m^2
CHI2NORM = 1    # CHI2 normalization; for bulk: LSUPER * 52.9177210903 for pm^2/V; for surface: 1

MATERIAL = {
    '1 QL': {
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

MATERIAL['1 QL']['chi2']['xyy'] = {'energy': MATERIAL['1 QL']['chi2']['xxx']['energy'], 'data': -MATERIAL['1 QL']['chi2']['xxx']['data']}
MATERIAL['1 QL']['chi2']['yxy'] = {'energy': MATERIAL['1 QL']['chi2']['xxx']['energy'], 'data': -MATERIAL['1 QL']['chi2']['xxx']['data']}
MATERIAL['1 QL']['chi2']['yyz'] = MATERIAL['1 QL']['chi2']['xxz']
MATERIAL['1 QL']['chi2']['zyy'] = MATERIAL['1 QL']['chi2']['zxx']

################################################################################
## Data: Experiment
################################################################################

exp_e, exp_rpp, exp_rsp, exp_rps = np.loadtxt('example/reference/experiment.dat', unpack=True)

EXP = {
    '1 QL': {
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
    widget = CustomWidget(EXP)
    widget.show()
    app.exec_()
