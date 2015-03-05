""" A quick and dirty script to convert an xyz file to and from coordinate """

import sys
import os
import numpy as np
from scipy import constants

def conversion():
    """ Does the actual conversion """
    in_file = sys.argv[1]
    target = sys.argv[2]
    if target == "bohr": # converts to bohrs
        factor = constants.angstrom / constants.value("Bohr radius")
    elif target == "angstrom": # converts to angstroms
        factor = constants.value("Bohr radius") / constants.angstrom
    else:
        print "Bad choice!"
        exit(0)
    source = np.loadtxt(in_file)
    new = source * factor
    save_coords(in_file, new)

def save_coords(output, data):
    """ Saves the new coordinates to a file """
    out_file = os.path.splitext(output)[0] + "_angstrom.xyz"
    np.savetxt(out_file, data, fmt=('% 2.15e'), delimiter='        ')

conversion()
