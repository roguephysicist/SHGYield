SHG Yield for Semiconductor Surfaces
====================================

[![DOI](https://zenodo.org/badge/11697217.svg)](https://zenodo.org/badge/latestdoi/11697217)

SHGYield is a python program for calculating the second-harmonic generation (SHG) yield for semiconductor surfaces. This is useful for conducting theoretical second-harmonic spectroscopy of any crystalline semiconductor or nanostructure. The physical theory behind this phenomenon is presented in the references below. We have applied this software to several silicon test cases and have achieved very accurate results that compare well with experiment.

These calculations can predict effects in novel materials that are still under study, or save valuable time and resources for the experimentalist. For example, the figure below is an overview of the angular dependence of the reflected SHG Yield from the Si(111)(1x1)H surface. Experimentalists will find this very useful, as they can plan the experiment accordingly in order to optimize the output signal strength and polarization.

![An overview of the angular dependence of the SHG Yield for the Si(111)(1x1)H surface](paper/3D-Si1x1.png)

You must first calculate the different components of the nonlinear susceptibility tensor. The theory surrounding this problem is still being developed, and there are many ways to go about it. We leave it to the reader to find the best method for their particular problem. As an example, we use [ABINIT](http://www.abinit.org) to calculate the electron density and matrix elements and then [TINIBA](https://github.com/bemese/tiniba) to calculate the tensor components. 


References
------------------------------------

The complete theory is derived step-by-step in my [Ph.D. thesis](https://github.com/roguephysicist/thesis-phd). This software has been used in the following publications:

* [Front. Mater. 4:12 (2017)](http://journal.frontiersin.org/article/10.3389/fmats.2017.00012/abstract)
* [Phys. Rev. B 94, 115314 (2016)](https://doi.org/10.1103/PhysRevB.94.115314)
* [Phys. Rev. B 93, 235304 (2016)](https://doi.org/10.1103/PhysRevB.93.235304)
* [Phys. Rev. B 91, 075302 (2015)](https://doi.org/10.1103/PhysRevB.91.075302)
* [arXiv:1604.07722 (2016)](https://arxiv.org/abs/1604.07722)
* [Theoretical Optical Second-Harmonic Calculations for Surfaces (2016)](https://doi.org/10.13140/RG.2.2.35619.66082)


Installation
------------------------------------

SHGYield has been tested with Python 2 and 3, and Anaconda Python 4+ on Mac OS X and Linux. It should work on any system with the required Python packages installed. 

Python requirements:
`sys`, `math`, `numpy`, `scipy`

Usage:
`python shgyield.py <sample.in>`

A sample input file `sample.in` and data set `sample-data/` are included for your convenience.


License
------------------------------------

Copyright 2017 Sean M. Anderson and Bernardo S. Mendoza.

SHGYield is free software made available under the BSD-3-Clause License. For details please see the LICENSE file.
