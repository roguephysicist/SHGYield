""" A quick and dirty script to convert an xyz file to and from coordinate """

import sys
import os
from scipy import constants
from numpy import loadtxt, savetxt

def conversion():
    """ Does the actual conversion """
    target = sys.argv[1]
    in_file = sys.argv[2]
    if target == "b": # converts to bohrs
        factor = constants.angstrom / constants.value("Bohr radius")
    elif target == "a": # converts to angstroms
        factor = constants.value("Bohr radius") / constants.angstrom
    else:
        print "Fuck you! You done fucked up!"
        exit(0)
    source = loadtxt(in_file)
    new = source * factor
    save_coords(in_file, new)

def save_coords(output, data):
    """ Saves the new coordinates to a file """
    out_file = os.path.splitext(output)[0] + "_mod.xyz"
    savetxt(out_file, data, fmt=('%5.14f'), delimiter='\t')

conversion()
