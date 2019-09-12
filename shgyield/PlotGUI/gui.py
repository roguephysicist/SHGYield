import numpy as np

from PyQt5 import QtGui, QtCore

import pyqtgraph as pg

import shgyield.shg as shg
from shgyield.PlotGUI.QtLayout import Ui_CustomWidget

COLORBG = '#EAEAEA'
COLOR00 = '#4C72B0'
COLOR01 = '#C44E52'
COLOR02 = '#55A868'
COLOR03 = '#64B5CD'
COLOR04 = '#DD8452'
COLOR05 = '#8172B3'
COLOR06 = '#937860'
COLOR07 = '#DA8BC3'
COLOR08 = '#CCB974'

pg.setConfigOption('background', COLORBG)
pg.setConfigOption('foreground', 'k')
pg.setConfigOption('antialias', True)
pg.setConfigOption('useWeave', True)


class CustomWidget(QtGui.QWidget):

    def __init__(self, epolar, espect, eps1, eps2, eps3, chi2, theta, phi, gamma, beps, bchi, bout, exp, parent=None):
        super(CustomWidget, self).__init__(parent=parent)

        # set up the form class as a `ui` attribute
        self.ui = Ui_CustomWidget()
        self.ui.setupUi(self)

        self.ui.box_energy_polar.setValue(epolar)
        self.ui.box_energy_spect_min.setValue(espect[0])
        self.ui.box_energy_spect_max.setValue(espect[1])
        self.ui.sld_angle_theta.setValue(theta)
        self.ui.sld_angle_phi.setValue(phi)
        self.ui.sld_angle_gamma.setValue(gamma)
        self.ui.sld_broad_eps.setValue(beps)
        self.ui.sld_broad_chi.setValue(bchi)
        self.ui.sld_broad_out.setValue(bout)

        self.increment = 0.01
        self.energy = np.arange(self.ui.box_energy_spect_min.value(), self.ui.box_energy_spect_max.value()+self.increment, self.increment)
        self.eps_m1 = eps1
        self.eps_m2 = eps2
        self.eps_m3 = eps3
        self.chi2r = chi2

        polar = shg.shgyield(energy =    self.ui.box_energy_polar.value(),
                             eps_m1 =    self.eps_m1,
                             eps_m2 =    self.eps_m2,
                             eps_m3 =    self.eps_m3,
                             chi2 =      self.chi2r,
                             theta =     self.ui.sld_angle_theta.value(),
                             phi =       np.arange(0, 361),
                             gamma =     self.ui.sld_angle_gamma.value(),
                             thick =     10,
                             sigma_eps = self.ui.sld_broad_eps.value(),
                             sigma_chi = self.ui.sld_broad_chi.value(),
                             sigma_out = self.ui.sld_broad_out.value())

        spect = shg.shgyield(energy =    self.energy,
                             eps_m1 =    eps1,
                             eps_m2 =    eps2,
                             eps_m3 =    eps3,
                             chi2 =      chi2,
                             theta =     self.ui.sld_angle_theta.value(),
                             phi =       self.ui.sld_angle_phi.value(),
                             gamma =     self.ui.sld_angle_gamma.value(),
                             thick =     10,
                             sigma_eps = self.ui.sld_broad_eps.value(),
                             sigma_chi = self.ui.sld_broad_chi.value(),
                             sigma_out = self.ui.sld_broad_out.value())

        # Combined tabs, polar
        self.plot_all_polar_rpp = self.ui.plot_all_polar_rpp.plot(x=polar['pp']*1e20*np.cos(np.radians(polar['phi'])), y=polar['pp']*1e20*np.sin(np.radians(polar['phi'])), pen=pg.mkPen(COLOR00, width=2))
        self.plot_all_polar_rsp = self.ui.plot_all_polar_rsp.plot(x=polar['sp']*1e20*np.cos(np.radians(polar['phi'])), y=polar['sp']*1e20*np.sin(np.radians(polar['phi'])), pen=pg.mkPen(COLOR01, width=2))
        self.plot_all_polar_rps = self.ui.plot_all_polar_rps.plot(x=polar['ps']*1e20*np.cos(np.radians(polar['phi'])), y=polar['ps']*1e20*np.sin(np.radians(polar['phi'])), pen=pg.mkPen(COLOR02, width=2))
        self.plot_all_polar_rss = self.ui.plot_all_polar_rss.plot(x=polar['ss']*1e20*np.cos(np.radians(polar['phi'])), y=polar['ss']*1e20*np.sin(np.radians(polar['phi'])), pen=pg.mkPen(COLOR03, width=2))
        self.plot_all_polar_rpp.getViewBox().setAspectLocked(True)
        self.plot_all_polar_rsp.getViewBox().setAspectLocked(True)
        self.plot_all_polar_rps.getViewBox().setAspectLocked(True)
        self.plot_all_polar_rss.getViewBox().setAspectLocked(True)
        self.plot_all_polar_rpp.getViewBox().disableAutoRange()
        self.plot_all_polar_rsp.getViewBox().disableAutoRange()
        self.plot_all_polar_rps.getViewBox().disableAutoRange()
        self.plot_all_polar_rss.getViewBox().disableAutoRange()
        self.plot_all_polar_rpp_marker = self.ui.plot_all_polar_rpp.plot(x=[polar['pp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['pp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_all_polar_rsp_marker = self.ui.plot_all_polar_rsp.plot(x=[polar['sp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['sp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_all_polar_rps_marker = self.ui.plot_all_polar_rps.plot(x=[polar['ps'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['ps'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_all_polar_rss_marker = self.ui.plot_all_polar_rss.plot(x=[polar['ss'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['ss'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))

        # Combined tabs, spectra
        self.plot_all_spect_rpp = self.ui.plot_all_spect_rpp.plot(x=spect['energy'], y=spect['pp']*1e20, pen=pg.mkPen(COLOR00, width=2))
        self.plot_all_spect_rsp = self.ui.plot_all_spect_rsp.plot(x=spect['energy'], y=spect['sp']*1e20, pen=pg.mkPen(COLOR01, width=2))
        self.plot_all_spect_rps = self.ui.plot_all_spect_rps.plot(x=spect['energy'], y=spect['ps']*1e20, pen=pg.mkPen(COLOR02, width=2))
        self.plot_all_spect_rss = self.ui.plot_all_spect_rss.plot(x=spect['energy'], y=spect['ss']*1e20, pen=pg.mkPen(COLOR03, width=2))
        self.plot_all_spect_rpp.getViewBox().setXLink(self.plot_all_spect_rsp.getViewBox())
        self.plot_all_spect_rsp.getViewBox().setXLink(self.plot_all_spect_rps.getViewBox())
        self.plot_all_spect_rps.getViewBox().setXLink(self.plot_all_spect_rss.getViewBox())
        self.plot_all_spect_rss.getViewBox().setXLink(self.plot_all_spect_rpp.getViewBox())
        self.plot_all_spect_rpp_marker = self.ui.plot_all_spect_rpp.plot(x=[polar['energy']], y=[spect['pp'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_all_spect_rsp_marker = self.ui.plot_all_spect_rsp.plot(x=[polar['energy']], y=[spect['sp'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_all_spect_rps_marker = self.ui.plot_all_spect_rps.plot(x=[polar['energy']], y=[spect['ps'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_all_spect_rss_marker = self.ui.plot_all_spect_rss.plot(x=[polar['energy']], y=[spect['ss'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))

        # RpP, theory vs. experiment
        self.plot_rpp_polar_exp = self.ui.plot_rpp_polar_exp.plot()
        self.plot_rpp_polar_exp.getViewBox().setAspectLocked(True)
        self.plot_rpp_spect_exp = self.ui.plot_rpp_spect_exp.plot(x=exp['energy'], y=exp['pp']*1e20, pen=None, symbol='o', symbolSize=5, symbolBrush=pg.mkBrush(COLOR07), symbolPen=pg.mkPen('k', width=0.5))
        self.plot_rpp_polar = self.ui.plot_rpp_polar.plot(x=polar['pp']*1e20*np.cos(np.radians(polar['phi'])), y=polar['pp']*1e20*np.sin(np.radians(polar['phi'])), pen=pg.mkPen(COLOR00, width=2))
        self.plot_rpp_polar.getViewBox().setAspectLocked(True)
        self.plot_rpp_spect = self.ui.plot_rpp_spect.plot(x=spect['energy'], y=spect['pp']*1e20, pen=pg.mkPen(COLOR00, width=2))
        self.plot_rpp_spect.getViewBox().setXLink(self.plot_rpp_spect_exp.getViewBox())
        self.plot_rpp_polar_marker = self.ui.plot_rpp_polar.plot(x=[polar['pp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['pp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_rpp_spect_marker = self.ui.plot_rpp_spect.plot(x=[polar['energy']], y=[spect['pp'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))

        # RsP, theory vs. experiment
        self.plot_rsp_polar_exp = self.ui.plot_rsp_polar_exp.plot()
        self.plot_rsp_polar_exp.getViewBox().setAspectLocked(True)
        self.plot_rsp_spect_exp = self.ui.plot_rsp_spect_exp.plot(x=exp['energy'], y=exp['sp']*1e20, pen=None, symbol='o', symbolSize=5, symbolBrush=pg.mkBrush(COLOR07), symbolPen=pg.mkPen('k', width=0.5))
        self.plot_rsp_polar = self.ui.plot_rsp_polar.plot(x=polar['sp']*1e20*np.cos(np.radians(polar['phi'])), y=polar['sp']*1e20*np.sin(np.radians(polar['phi'])), pen=pg.mkPen(COLOR01, width=2))
        self.plot_rsp_polar.getViewBox().setAspectLocked(True)
        self.plot_rsp_spect = self.ui.plot_rsp_spect.plot(x=spect['energy'], y=spect['sp']*1e20, pen=pg.mkPen(COLOR01, width=2))
        self.plot_rsp_spect.getViewBox().setXLink(self.plot_rsp_spect_exp.getViewBox())
        self.plot_rsp_polar_marker = self.ui.plot_rsp_polar.plot(x=[polar['sp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['sp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_rsp_spect_marker = self.ui.plot_rsp_spect.plot(x=[polar['energy']], y=[spect['sp'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        
        # RpS, theory vs. experiment
        self.plot_rps_polar_exp = self.ui.plot_rps_polar_exp.plot()
        self.plot_rps_polar_exp.getViewBox().setAspectLocked(True)
        self.plot_rps_spect_exp = self.ui.plot_rps_spect_exp.plot(x=exp['energy'], y=exp['ps']*1e20, pen=None, symbol='o', symbolSize=5, symbolBrush=pg.mkBrush(COLOR07), symbolPen=pg.mkPen('k', width=0.5))
        self.plot_rps_polar = self.ui.plot_rps_polar.plot(x=polar['ps']*1e20*np.cos(np.radians(polar['phi'])), y=polar['ps']*1e20*np.sin(np.radians(polar['phi'])), pen=pg.mkPen(COLOR02, width=2))
        self.plot_rps_polar.getViewBox().setAspectLocked(True)
        self.plot_rps_spect = self.ui.plot_rps_spect.plot(x=spect['energy'], y=spect['ps']*1e20, pen=pg.mkPen(COLOR02, width=2))
        self.plot_rps_spect.getViewBox().setXLink(self.plot_rps_spect_exp.getViewBox())
        self.plot_rps_polar_marker = self.ui.plot_rps_polar.plot(x=[polar['ps'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['ps'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_rps_spect_marker = self.ui.plot_rps_spect.plot(x=[polar['energy']], y=[spect['ps'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        
        # RsS, theory vs. experiment
        self.plot_rss_polar_exp = self.ui.plot_rss_polar_exp.plot()
        self.plot_rss_polar_exp.getViewBox().setAspectLocked(True)
        self.plot_rss_spect_exp = self.ui.plot_rss_spect_exp.plot(x=exp['energy'], y=exp['ss']*1e20, pen=None, symbol='o', symbolSize=5, symbolBrush=pg.mkBrush(COLOR07), symbolPen=pg.mkPen('k', width=0.5))
        self.plot_rss_polar = self.ui.plot_rss_polar.plot(x=polar['ss']*1e20*np.cos(np.radians(polar['phi'])), y=polar['ss']*1e20*np.sin(np.radians(polar['phi'])), pen=pg.mkPen(COLOR03, width=2))
        self.plot_rss_polar.getViewBox().setAspectLocked(True)
        self.plot_rss_spect = self.ui.plot_rss_spect.plot(x=spect['energy'], y=spect['ss']*1e20, pen=pg.mkPen(COLOR03, width=2))
        self.plot_rss_spect.getViewBox().setXLink(self.plot_rss_spect_exp.getViewBox())
        self.plot_rss_polar_marker = self.ui.plot_rss_polar.plot(x=[polar['ss'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['ss'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_rss_spect_marker = self.ui.plot_rss_spect.plot(x=[polar['energy']], y=[spect['ss'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))

        # Updating
        self.ui.box_energy_polar.valueChanged.connect(self.update_plot)
        self.ui.box_energy_spect_min.valueChanged.connect(self.update_plot)
        self.ui.box_energy_spect_max.valueChanged.connect(self.update_plot)
        
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


    def update_plot(self):
        self.energy = np.arange(self.ui.box_energy_spect_min.value(), self.ui.box_energy_spect_max.value()+self.increment, self.increment)
        polar = shg.shgyield(energy =    self.ui.box_energy_polar.value(),
                             eps_m1 =    self.eps_m1,
                             eps_m2 =    self.eps_m2,
                             eps_m3 =    self.eps_m3,
                             chi2 =      self.chi2r,
                             theta =     self.ui.sld_angle_theta.value(),
                             phi =       np.arange(0, 361),
                             gamma =     self.ui.sld_angle_gamma.value(),
                             thick =     10,
                             sigma_eps = self.ui.sld_broad_eps.value(),
                             sigma_chi = self.ui.sld_broad_chi.value(),
                             sigma_out = self.ui.sld_broad_out.value())

        spect = shg.shgyield(energy =    self.energy,
                             eps_m1 =    self.eps_m1,
                             eps_m2 =    self.eps_m2,
                             eps_m3 =    self.eps_m3,
                             chi2 =      self.chi2r,
                             theta =     self.ui.sld_angle_theta.value(),
                             phi =       self.ui.sld_angle_phi.value(),
                             gamma =     self.ui.sld_angle_gamma.value(),
                             thick =     10,
                             sigma_eps = self.ui.sld_broad_eps.value(),
                             sigma_chi = self.ui.sld_broad_chi.value(),
                             sigma_out = self.ui.sld_broad_out.value())

        # Combined tab, polar
        self.plot_all_polar_rpp.setData(polar['pp']*1e20*np.cos(np.radians(polar['phi'])), polar['pp']*1e20*np.sin(np.radians(polar['phi'])))
        self.plot_all_polar_rsp.setData(polar['sp']*1e20*np.cos(np.radians(polar['phi'])), polar['sp']*1e20*np.sin(np.radians(polar['phi'])))
        self.plot_all_polar_rps.setData(polar['ps']*1e20*np.cos(np.radians(polar['phi'])), polar['ps']*1e20*np.sin(np.radians(polar['phi'])))
        self.plot_all_polar_rss.setData(polar['ss']*1e20*np.cos(np.radians(polar['phi'])), polar['ss']*1e20*np.sin(np.radians(polar['phi'])))
        self.plot_all_polar_rpp_marker.setData(x=[polar['pp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['pp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_all_polar_rsp_marker.setData(x=[polar['sp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['sp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_all_polar_rps_marker.setData(x=[polar['ps'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['ps'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_all_polar_rss_marker.setData(x=[polar['ss'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['ss'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))

        # Combined tab, polar
        self.plot_all_spect_rpp.setData(spect['energy'], spect['pp']*1e20)
        self.plot_all_spect_rsp.setData(spect['energy'], spect['sp']*1e20)
        self.plot_all_spect_rps.setData(spect['energy'], spect['ps']*1e20)
        self.plot_all_spect_rss.setData(spect['energy'], spect['ss']*1e20)
        self.plot_all_spect_rpp_marker.setData(x=[polar['energy']], y=[spect['pp'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_all_spect_rsp_marker.setData(x=[polar['energy']], y=[spect['sp'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_all_spect_rps_marker.setData(x=[polar['energy']], y=[spect['ps'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))
        self.plot_all_spect_rss_marker.setData(x=[polar['energy']], y=[spect['ss'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20], pen=None, symbolSize=12, symbol='o', symbolPen=pg.mkPen(None), symbolBrush=pg.mkBrush(204,0,0,150))

        # RpP, theory vs. experiment
        self.plot_rpp_polar.setData(polar['pp']*1e20*np.cos(np.radians(polar['phi'])), polar['pp']*1e20*np.sin(np.radians(polar['phi'])))
        self.plot_rpp_spect.setData(spect['energy'], spect['pp']*1e20)
        self.plot_rpp_spect_marker.setData(x=[polar['energy']], y=[spect['pp'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20])
        self.plot_rpp_polar_marker.setData(x=[polar['pp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['pp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))])
        
        # RsP, theory vs. experiment
        self.plot_rsp_polar.setData(polar['sp']*1e20*np.cos(np.radians(polar['phi'])), polar['sp']*1e20*np.sin(np.radians(polar['phi'])))
        self.plot_rsp_spect.setData(spect['energy'], spect['sp']*1e20)
        self.plot_rsp_spect_marker.setData(x=[polar['energy']], y=[spect['sp'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20])
        self.plot_rsp_polar_marker.setData(x=[polar['sp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['sp'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))])

        # RpS, theory vs. experiment
        self.plot_rps_polar.setData(polar['ps']*1e20*np.cos(np.radians(polar['phi'])), polar['ps']*1e20*np.sin(np.radians(polar['phi'])))
        self.plot_rps_spect.setData(spect['energy'], spect['ps']*1e20)
        self.plot_rps_spect_marker.setData(x=[polar['energy']], y=[spect['ps'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20])
        self.plot_rps_polar_marker.setData(x=[polar['ps'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['ps'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))])

        # RsS, theory vs. experiment
        self.plot_rss_polar.setData(polar['ss']*1e20*np.cos(np.radians(polar['phi'])), polar['ss']*1e20*np.sin(np.radians(polar['phi'])))
        self.plot_rss_spect.setData(spect['energy'], spect['ss']*1e20)
        self.plot_rss_spect_marker.setData(x=[polar['energy']], y=[spect['ss'][np.where(np.isclose(spect['energy'], polar['energy']))][0]*1e20])
        self.plot_rss_polar_marker.setData(x=[polar['ss'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.cos(np.radians(spect['phi']))], y=[polar['ss'][np.where(np.isclose(spect['phi'], polar['phi']))][0]*1e20*np.sin(np.radians(spect['phi']))])
