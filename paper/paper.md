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
   affiliation: 2
 - name: Bernardo S. Mendoza
   orcid: 0000-0002-8546-0262
   affiliation: 1
affiliations:
 - name: Centro de Investigaciones en Óptica, A. C.
   index: 1
 - name: None
   index: 2

date: 07 February 2017
bibliography: paper.bib
nocite: |
  @andersonFMATS17, @andersonPRB16b, @andersonPRB16a, @andersonPRB15, @andersonARXIV16, @andersonthesis
---

# Summary

`SHGYield.py` is a python program for calculating the surface second-harmonic generation (SSHG) yield (in reflectance) for semiconductor surfaces.

SSHG is an effective, nondestructive, and noninvasive probe for studying surface and interface properties, and even for characterizing buried interfaces and nanostructures. The high surface sensitivity of SSHG spectroscopy is due to the fact that within the dipole approximation, the bulk SHG response in centrosymmetric materials is identically zero. The SHG process can occur only at the surface where the inversion symmetry is broken.

This program has several potential applications and uses:

- determining and analyzing the physical origin of SSHG spectra  
- predicting and characterizing the radiated SSHG for interesting new materials  
- characterizing thin films based on measured SH spectra  
- allowing the experimenter to calculate and analyze the SSHG yield to optimize experiments  

For example, the figure below is an overview of the angular dependence of the reflected SHG Yield from the Si(111)(1x1)H surface. Experimentalists will find this very useful, as they can plan the experiment accordingly in order to optimize the output signal strength and polarization.

![An overview of the angular dependence of the SHG Yield for the Si(111)(1x1)H surface](../example/figures/3D-Si1x1.png)


# References

The complete theory is derived step-by-step in Ref. [@andersonPRB16b]. This software has been developed and used in the following publications:
