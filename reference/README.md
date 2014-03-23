BMS Reference Program
=====================

This is BMS' original program for calculating the nonlinear reflection
coefficients for silicon slabs. His program is written in FORTRAN and has 
the following input/outputs:

* fort.14 (input)contains the linear data,
* fort.25 (input)contains the chi2 values,
* fort.36 (output) contains Rpp, Rsp, Rps, and Rss,
* fort.38 (output) contains the xxz, zxx, and zzz components of Rpp,
* fort.91 and fort.93 (output) are for debugging only.

My program in python can replicate the results of the FORTRAN program exactly.
It uses the appropriately named files in the 'res/' directory and outputs to
the 'nrc/' directory.
