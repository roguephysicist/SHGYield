"""
A quick and dirty script to convert from an xyz file from bohrs to angstroms
"""

import sys
import os
from numpy import loadtxt, savetxt

IN_FILE = sys.argv[1]
OUT_FILE = os.path.splitext(IN_FILE)[0] + "_xangst.xyz"
BOHR_COORDS = loadtxt(IN_FILE)
ANGS_COORDS = BOHR_COORDS * 0.529177249
savetxt(OUT_FILE, ANGS_COORDS, fmt=('%5.14f'))
