#!/Users/sma/anaconda/bin/python
"""
Allows you to add stretch or shrink bonds between atoms in a xyz file.
* Make lists that contain what changes in order of change, like 
[2,3,+4]
[3,4,-4]
"""

import sys
import math

DELTA = 10

def read_atoms():
    """ reads atoms from file and converts strings to floats """
    atoms = []
    infile = sys.argv[1]
    strings = [line.strip().split() for line in open(infile)]
    for atom in strings:
        temp = [float(item) for item in atom]
        atoms.append(temp)
    return atoms

def atom_lengths(atom1, atom2):
    """ distances between each cartesian coordinate for two atoms """
    distx = atom2[0] - atom1[0]
    disty = atom2[1] - atom1[1]
    distz = atom2[2] - atom1[2]
    lengths = [distx, disty, distz]
    return lengths

def spherical_coords(atom1, atom2):
    """ converts cartesian to spherical coordinates """
    dist = atom_lengths(atom1, atom2)
    radius = (((atom2[0] - atom1[0]) ** 2) +
              ((atom2[1] - atom1[1]) ** 2) +
              ((atom2[2] - atom1[2]) ** 2)) ** 0.5
    theta = math.atan2(dist[1], dist[0])
    phi = math.acos(dist[2]/radius)
    coords = [radius, theta, phi]
    return coords

def cartesian_coords(radius, theta, phi):
    """ convertes spherical to cartesian coordinates """
    xcart = radius * math.cos(theta) * math.sin(phi)
    ycart = radius * math.sin(theta) * math.sin(phi)
    zcart = radius * math.cos(phi)
    deltas = [xcart, ycart, zcart]
    return deltas

def newcoords(origin, deltas):
    """ adds calculated deltas to original coordinate """
    xnew = deltas[0] + origin[0]
    ynew = deltas[1] + origin[1]
    znew = deltas[2] + origin[2]
    coords = [xnew, ynew, znew]
    return coords

def control(origin, target):
    """ obtains the new coordinate. put inside for loop """
    porcent = 1 + DELTA/100.0
    sphere = spherical_coords(origin, target)
    new_radius = sphere[0] * porcent
    diff = cartesian_coords(new_radius, sphere[1], sphere[2])
    new = newcoords(origin, diff)
    return new

ATOMS = read_atoms()
print control(ATOMS[2], ATOMS[1])



# def print_table():
#     """ prints table with atoms, symbol names, and bond lengths """
#     os.system('clear')
#     print "Available atoms and bond lengths (Bohrs):"
#     for idx, atom in enumerate(ATOMS):
#         print "{0} {1:18.14f}   {2:18.14f}    {3:18.14f}"\
#             .format(SYMBOL[idx].ljust(2),
#                     float(atom[0]),
#                     float(atom[1]),
#                     float(atom[2]))
#         if idx != len(ATOMS) - 1:
#             print '|' + 21 * 3 * '=' + "| {0}: {1:16.14f}"\
#                 .format(str(idx+1).zfill(2), BONDS[idx])
#    print

# print_table()
# SELECTION = input("Which bond do you want to modify? ")
# DELTA = input("By how much do you wish to modify it? (%) ")
# print "You have selected bond {0}".format(SELECTION)
# print "Old bond length = {0:18.14f}".format(SELECTION)
# percent = 1 + DELTA/100.0
