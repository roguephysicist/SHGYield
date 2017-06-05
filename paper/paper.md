---
title: 'SHGYield'
tags:
  - python
  - shg
  - second-harmonic
  - spectroscopy
  - semiconductor
  - surfaces
authors:
 - name: Sean M. Anderson
   orcid: 0000-0002-4809-921X
   affiliation: none
 - name: Bernardo S. Mendoza
   orcid: 0000-0002-8546-0262
   affiliation: 1
affiliations:
 - name: Centro de Investigaciones en Óptica, A. C.
   index: 1
date: 07 February 2017
bibliography: paper.bib
---

# Summary

`SHGYield.py` is a python program for calculating the surface second-harmonic generation (SSHG) yield (in reflectance) for semiconductor surfaces. 

SSHG is an effective, nondestructive, and noninvasive probe for studying surface and interface properties, and even for characterizing buried interfaces and nanostructures. The high surface sensitivity of SSHG spectroscopy is due to the fact that within the dipole approximation, the bulk SHG response in centrosymmetric materials is identically zero. The SHG process can occur only at the surface where the inversion symmetry is broken. 

This program has several potential applications and uses:
* determining and analyzing the physical origin of SSHG spectra
* predicting and characterizing the radiated SSHG for interesting new materials
* characterizing thin films based on measured SH spectra
* allowing the experimenter to calculate and analyze the SSHG yield to optimize experiments

For example, the figure below is an overview of the angular dependence of the reflected SHG Yield from the Si(111)(1x1)H surface. Experimentalists will find this very useful, as they can plan the experiment accordingly in order to optimize the output signal strength and polarization.

![An overview of the angular dependence of the SHG Yield for the Si(111)(1x1)H surface](../example/figures/3D-Si1x1.png)


# References

The complete theory is derived step-by-step in [Phys. Rev. B 94, 115314 (2016)](https://doi.org/10.1103/PhysRevB.94.115314). This software has been developed and used in the following publications:

* [Front. Mater. 4:12 (2017)](https://doi.org/10.3389/fmats.2017.00012)
* [Phys. Rev. B 94, 115314 (2016)](https://doi.org/10.1103/PhysRevB.94.115314)
* [Phys. Rev. B 93, 235304 (2016)](https://doi.org/10.1103/PhysRevB.93.235304)
* [Phys. Rev. B 91, 075302 (2015)](https://doi.org/10.1103/PhysRevB.91.075302)
* [arXiv:1604.07722 (2016)](https://arxiv.org/abs/1604.07722)
* [Theoretical Optical Second-Harmonic Calculations for Surfaces (2016)](https://doi.org/10.13140/RG.2.2.35619.66082)
