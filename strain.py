#!/Users/sma/anaconda/bin/python
"""
Allows you to add stretch or shrink bonds between atoms in a xyz file.
"""

import sys
import math

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

def add_deltas(origin, deltas):
    """ adds calculated deltas to original coordinate """
    xnew = deltas[0] + origin[0]
    ynew = deltas[1] + origin[1]
    znew = deltas[2] + origin[2]
    coords = [xnew, ynew, znew]
    return coords

def control(origin, target, porcent):
    """ obtains the new coordinate. put inside for loop """
    delta = 1 + porcent/100.0
    sphere = spherical_coords(origin, target)
    new_radius = sphere[0] * delta
    diff = cartesian_coords(new_radius, sphere[1], sphere[2])
    new = add_deltas(origin, diff)
    return new

#STRAIN = [[3, 4], [2, 4]]
#OUTFILE = 'plus01.xyz'
#STRAIN = [[5, 4], [4, 4], [3, 4], [2, 4]]
#OUTFILE = 'plus02.xyz'
#STRAIN = [[7, 4], [6, 4], [5, 4], [4, 4], [3, 4], [2, 4]]
#OUTFILE = 'plus03.xyz'
#STRAIN = [[9, 4], [8, 4], [7, 4], [6, 4], [5, 4], [4, 4], [3, 4], [2, 4]]
#OUTFILE = 'plus04.xyz'
STRAIN = [[11, 4], [10, 4], [9, 4], [8, 4], [7, 4], [6, 4],  [5, 4], [4, 4], [2, 4], [3, 4]]
OUTFILE = 'plus05.xyz'
ATOMS = read_atoms()
for trans in STRAIN:
    NEWATOM = control(ATOMS[trans[0]], ATOMS[trans[0] - 1], trans[1])
    DELTAS = atom_lengths(ATOMS[trans[0] - 1], NEWATOM)
    for num in xrange(trans[0]):
        newatom = add_deltas(ATOMS[num], DELTAS)
        ATOMS[num] = newatom
with open(OUTFILE, 'w') as file:
    for item in ATOMS:
        file.write("{0:17.14f}   {1:17.14f}    {2:18.14f}".format(item[0], item[1], item[2]) + '\n')
