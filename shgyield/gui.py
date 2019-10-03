import numpy as np
import pyqtgraph as pg

from PyQt5 import QtCore, QtGui, QtWidgets

import shgyield.shg as shg
from shgyield.QtLayout import Ui_CustomWidget


pg.setConfigOption('background', '#EAEAEA')
pg.setConfigOption('foreground', 'k')
pg.setConfigOption('antialias', True)
pg.setConfigOption('useWeave', True)


class CustomWidget(QtGui.QWidget):

    def __init__(self, params, espect, epolar, experiment, material, scale, parent=None):
        super(CustomWidget, self).__init__(parent=parent)

        # set up the form class as a `ui` attribute
        self.ui = Ui_CustomWidget()
        self.ui.setupUi(self)

        self.params = params
        self.scale = scale
        self.experiment = experiment
        self.material = material

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
                self.widgets['tab_' + pol]['thr_polar'].addItem(
                    pg.InfiniteLine(angle = ang, pen = pg.mkPen(color=0.8, width=0.5, style=QtCore.Qt.DashLine)))
                self.widgets['tab_' + pol]['exp_polar'].addItem(
                    pg.InfiniteLine(angle = ang, pen = pg.mkPen(color=0.8, width=0.5, style=QtCore.Qt.DashLine)))

        # calculation time!

################################################################################
################################################################################
        
        # loc = 0
        # self.gapgui = {}
        # self.group_gap = QtWidgets.QGroupBox(self.ui.widget)
        # self.group_gap.setObjectName("group_gap")
        # self.grid_gap = QtWidgets.QGridLayout(self.group_gap)
        # self.grid_gap.setObjectName("grid_gap")
        # self.ui.verticalLayout.addWidget(self.group_gap)
        # self.group_gap.setTitle(QtCore.QCoreApplication.translate("CustomWidget", "Gap"))
        # for case in self.material.keys():
        #     self.gapgui['label_gap_' + case] = QtWidgets.QLabel(self.group_gap)
        #     self.gapgui['label_gap_' + case].setObjectName('label_gap_' + case)
        #     self.gapgui['label_gap_' + case].setText(case)

        #     self.gapgui['box_gap_' + case] = QtWidgets.QDoubleSpinBox(self.group_gap)
        #     self.gapgui['box_gap_' + case].setObjectName('box_gap_' + case)
        #     self.gapgui['box_gap_' + case].setSingleStep(0.01)
        #     self.gapgui['box_gap_' + case].setSuffix(" eV")
        #     self.gapgui['box_gap_' + case].setMaximum(self.material[case]['gap']['rng'].max())
        #     self.gapgui['box_gap_' + case].setMinimum(self.material[case]['gap']['rng'].min())

        #     self.gapgui['box_gap_' + case].setValue(self.material[case]['gap']['gw'])

        #     self.grid_gap.addWidget(self.gapgui['label_gap_' + case], loc, 0, 1, 1)
        #     self.grid_gap.addWidget(self.gapgui['box_gap_' + case], loc, 1, 1, 1)

        #     self.gapgui['box_gap_' + case].valueChanged.connect(self.update_plot)
        #     self.gapgui['box_gap_' + case].valueChanged.connect(self.update_plot)
        #     loc += 1

################################################################################
################################################################################

        # init: energy
        self.einc = 0.01
        self.epolar = epolar
        self.erng = espect
        # self.erng = np.arange(self.ui.box_energy_spect_min.value(),
        #                       self.ui.box_energy_spect_max.value()+self.einc,
        #                       self.einc)
        self.ui.box_energy_polar.setValue(self.epolar)
        self.ui.box_energy_polar.setMinimum(self.erng.min())
        self.ui.box_energy_polar.setMaximum(self.erng.max())
        self.ui.box_energy_polar.setSingleStep(self.einc)

        self.ui.sld_angle_theta.setValue(self.params['theta'])
        self.ui.sld_angle_phi.setValue(self.params['phi'])
        self.ui.sld_broad_eps.setValue(self.params['broad']['eps'])
        self.ui.sld_broad_chi.setValue(self.params['broad']['chi'])
        self.ui.sld_broad_out.setValue(self.params['broad']['out'])

        # plots for experimental data, that do not change or need calculations
        deco = 0
        for case, measurement in self.experiment.items():
            for kind, values in measurement.items():
                for polarization, data in values.items():
                    if kind == 'polar':
                        self.widgets['tab_' + polarization]['exp_' + kind].plot(
                            self.trans_polar(data['phi'], data['data']),
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
                                y = data['data'],
                                height = data['stdev'],
                                pen = pg.mkPen(color=self.colors[deco], width=1)))
                        self.widgets['tab_' + polarization]['exp_' + kind].plot(
                            x = data['energy'],
                            y = data['data'],
                            pen = None,
                            symbol = deco,
                            symbolSize = 5,
                            symbolBrush = pg.mkBrush(color=self.colors[deco]),
                            symbolPen = pg.mkPen('k', width=0.5),
                            name = case)
            deco += 1

################################################################################
################################################################################

        deco = 0
        self.test = {}
        for case in self.material.keys():
            self.polar = shg.shgyield(energy =    self.ui.box_energy_polar.value(),
                                      eps_m1 =    self.material[case]['medium 1']['eps'],
                                      # eps_m2 =    self.material[case]['medium 2']['{:4.2f}'.format(self.gapgui['box_gap_' + case].value())]['eps'],
                                      eps_m2 =    self.material[case]['medium 2']['eps'],
                                      eps_m3 =    self.material[case]['medium 3']['eps'],
                                      # chi2 =      self.material[case]['medium 2']['{:4.2f}'.format(self.gapgui['box_gap_' + case].value())]['chi2'],
                                      chi2 =      self.material[case]['medium 2']['chi2'],
                                      theta =     self.ui.sld_angle_theta.value(),
                                      phi =       np.arange(0, 361),
                                      gamma =     self.ui.sld_angle_gamma.value(),
                                      thick =     self.material[case]['lslab'],
                                      sigma_eps = self.ui.sld_broad_eps.value(),
                                      sigma_chi = self.ui.sld_broad_chi.value(),
                                      sigma_out = 0)

            self.spect = shg.shgyield(energy =    self.erng,
                                      eps_m1 =    self.material[case]['medium 1']['eps'],
                                      # eps_m2 =    self.material[case]['medium 2']['{:4.2f}'.format(self.gapgui['box_gap_' + case].value())]['eps'],
                                      eps_m2 =    self.material[case]['medium 2']['eps'],
                                      eps_m3 =    self.material[case]['medium 3']['eps'],
                                      # chi2 =      self.material[case]['medium 2']['{:4.2f}'.format(self.gapgui['box_gap_' + case].value())]['chi2'],
                                      chi2 =      self.material[case]['medium 2']['chi2'],
                                      theta =     self.ui.sld_angle_theta.value(),
                                      phi =       self.ui.sld_angle_phi.value(),
                                      gamma =     self.ui.sld_angle_gamma.value(),
                                      thick =     self.material[case]['lslab'],
                                      sigma_eps = self.ui.sld_broad_eps.value(),
                                      sigma_chi = self.ui.sld_broad_chi.value(),
                                      sigma_out = self.ui.sld_broad_out.value())

            self.pidx = np.where(np.isclose(self.spect['phi'], self.polar['phi']))
            self.eidx = np.where(np.isclose(self.spect['energy'], self.polar['energy']))

            self.test[case] = {
                pol: {
                    'test1': self.widgets['tab_polar'][pol].plot(self.trans_polar(self.polar['phi'], self.polar[pol]*self.scale), pen=pg.mkPen(color=self.colors[deco], width=1.5), name=case),
                    'test2': self.widgets['tab_polar'][pol].plot(self.trans_polar2(self.spect['phi'], self.polar[pol][self.pidx][0]*self.scale), **self.marker),
                    'test3': self.widgets['tab_spect'][pol].plot(x=self.spect['energy'], y=self.spect[pol]*self.scale, pen=pg.mkPen(color=self.colors[deco], width=1.5), name=case),
                    'test4': self.widgets['tab_spect'][pol].addLine(x=self.polar['energy'], pen=pg.mkPen(color=(204,0,0,150), width=2)),
                    'test5': self.widgets['tab_' + pol]['thr_polar'].plot(self.trans_polar(self.polar['phi'], self.polar[pol]*self.scale), pen=pg.mkPen(color=self.colors[deco], width=1.5), name=case),
                    'test6': self.widgets['tab_' + pol]['thr_polar'].plot(self.trans_polar2(self.spect['phi'], self.polar[pol][self.pidx][0]*self.scale), **self.marker),
                    'test7': self.widgets['tab_' + pol]['thr_spect'].plot(x = self.spect['energy'], y = self.spect[pol]*self.scale, pen=pg.mkPen(color=self.colors[deco], width=1.5), name=case),
                    'test8': self.widgets['tab_' + pol]['thr_spect'].addLine(x = self.polar['energy'], pen = pg.mkPen(color=(204,0,0,150), width=2))
                } for pol in ['pp', 'sp', 'ps', 'ss']
            } 
            
            deco += 1

################################################################################
################################################################################

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
            self.widgets['tab_' + pol]['exp_polar'].setAspectLocked(True)
            self.widgets['tab_' + pol]['exp_polar'].disableAutoRange()

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
        for case in self.material.keys():
            polar = shg.shgyield(energy =    self.ui.box_energy_polar.value(),
                                 eps_m1 =    self.material[case]['medium 1']['eps'],
                                 # eps_m2 =    self.material[case]['medium 2']['{:4.2f}'.format(self.gapgui['box_gap_' + case].value())]['eps'],
                                 eps_m2 =    self.material[case]['medium 2']['eps'],
                                 eps_m3 =    self.material[case]['medium 3']['eps'],
                                 # chi2 =      self.material[case]['medium 2']['{:4.2f}'.format(self.gapgui['box_gap_' + case].value())]['chi2'],
                                 chi2 =      self.material[case]['medium 2']['chi2'],
                                 theta =     self.ui.sld_angle_theta.value(),
                                 phi =       np.arange(0, 361),
                                 gamma =     self.ui.sld_angle_gamma.value(),
                                 thick =     self.material[case]['lslab'],
                                 sigma_eps = self.ui.sld_broad_eps.value(),
                                 sigma_chi = self.ui.sld_broad_chi.value(),
                                 sigma_out = 0)

            spect = shg.shgyield(energy =    self.erng,
                                 eps_m1 =    self.material[case]['medium 1']['eps'],
                                 # eps_m2 =    self.material[case]['medium 2']['{:4.2f}'.format(self.gapgui['box_gap_' + case].value())]['eps'],
                                 eps_m2 =    self.material[case]['medium 2']['eps'],
                                 eps_m3 =    self.material[case]['medium 3']['eps'],
                                 # chi2 =      self.material[case]['medium 2']['{:4.2f}'.format(self.gapgui['box_gap_' + case].value())]['chi2'],
                                 chi2 =      self.material[case]['medium 2']['chi2'],
                                 theta =     self.ui.sld_angle_theta.value(),
                                 phi =       self.ui.sld_angle_phi.value(),
                                 gamma =     self.ui.sld_angle_gamma.value(),
                                 thick =     self.material[case]['lslab'],
                                 sigma_eps = self.ui.sld_broad_eps.value(),
                                 sigma_chi = self.ui.sld_broad_chi.value(),
                                 sigma_out = self.ui.sld_broad_out.value())

            pidx = np.where(np.isclose(spect['phi'], polar['phi']))
            eidx = np.where(np.isclose(spect['energy'], polar['energy']))

            for pol in ['pp', 'sp', 'ps', 'ss']:

                self.test[case][pol]['test1'].setData(self.trans_polar(polar['phi'], polar[pol]*self.scale)),
                self.test[case][pol]['test2'].setData(self.trans_polar2(spect['phi'], polar[pol][pidx][0]*self.scale)),
                self.test[case][pol]['test3'].setData(x = spect['energy'], y = spect[pol]*self.scale),
                self.test[case][pol]['test4'].setValue(polar['energy']),
                self.test[case][pol]['test5'].setData(self.trans_polar(polar['phi'], polar[pol]*self.scale)),
                self.test[case][pol]['test6'].setData(self.trans_polar2(spect['phi'], polar[pol][pidx][0]*self.scale)),
                self.test[case][pol]['test7'].setData(x = spect['energy'], y = spect[pol]*self.scale),
                self.test[case][pol]['test8'].setValue(polar['energy'])
