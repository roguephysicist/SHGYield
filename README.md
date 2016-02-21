SHG Yield for Semiconductor Surfaces
====================================

shgyield.py is a python program designed to calculate the nonlinear reflection
coefficient for semiconductor surfaces. It works in conjunction with the
matrix elements calculated with ABINIT, an open source ab initio software, and
TINIBA, our in-house optical calculation software.

The work coded in this software can be found in an upcoming publication and is
explicitly derived in 'phd-thesis'.

Tested with Anaconda Python.

requirements:
sys, math, numpy, scipy

usage:
python shgyield.py <sample.in>
