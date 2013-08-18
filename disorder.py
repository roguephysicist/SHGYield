"""
This python script disorders an input XYX file (in bohrs) and lets you choose
which atoms (rows) you want disordered and by how much.
"""

from random import random
from math import sin, cos
from scipy import constants
from numpy import loadtxt, savetxt

DISORDER_PATH = "./sample_data/disorder/si_h_14.xyz"
BOND_LENGTH = 2.351
ATOMS_DISORDERED = [2, 3, 4]
DISORDER_AMOUNT = [0.1, 0.3, 0.9]
REPEAT = 1

def disorder():
    """ extracts selected rows from input, converts to matrix,
    disorders, and returns matrix which is written to file """
    count = 1
    chosen = [x - 1 for x in ATOMS_DISORDERED]
    blbohr = BOND_LENGTH * constants.angstrom / constants.value("Bohr radius")
    xyz = loadtxt(DISORDER_PATH)
    data = xyz[chosen]
    new_xyz = xyz.copy()
    final = zip(data, DISORDER_AMOUNT)
    while count <= REPEAT:
        for atom in range(0, len(final)):
            polar = random() * constants.pi
            azimuthal = 2 * polar
            new_xyz[chosen[atom], 0] = final[atom][0][0] + (final[atom][1] *
                    (blbohr / 2) * sin(polar) * cos(azimuthal))
            new_xyz[chosen[atom], 1] = final[atom][0][1] + (final[atom][1] *
                    (blbohr / 2) * sin(polar) * sin(azimuthal))
            new_xyz[chosen[atom], 2] = final[atom][0][2] + (final[atom][1] *
                    (blbohr / 2) * cos(polar))
            out = "sample_data/disorder/si_h_14_" + str(count).zfill(3) + ".xyz"
            savetxt(out, new_xyz, fmt=('%5.14e'), delimiter='\t')
        count += 1
    return polar

disorder()
