SHG Yield for Semiconductor Surfaces
====================================

shgyield.py is a python script designed to calculate the nonlinear reflection
coefficient for semiconductor surfaces. It works in conjunction with the
matrix elements calculated with ABINIT, an open source ab initio software, and
TINIBA, our in-house optical calculation software.

The work coded in this software can be found in the following publications:

* [Phys. Rev. B 94, 115314](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.115314)
* [Phys. Rev. B 93, 235304](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.235304)
* [arXiv:1604.07722](https://arxiv.org/abs/1604.07722)

All the theory is derived step-by-step in my Ph.D. thesis, which you can find in
the appropriately named repository right here on Github.

This script has been tested with Python 2.7.11 and Anaconda 4.0.0 on Mac OS X
and Linux. It should work on any system with the required Python packages
installed.

The script reads an input file that specifies all the necessary filenames,
angles, broadening, and other variables. It allows the user to include the
effects of multiple reflections if desired, and select several parameters for
calculating these effects. It will automatically try to read all 18 independent
components of the nonlinear susceptibility (for SHG). If you do not wish to
calculate all the components or know which of them are zero due to symmetry
relations, you can just comment out the lines in the input file.

requirements:
sys, math, numpy, scipy

usage:
`python shgyield.py <sample.in>`
