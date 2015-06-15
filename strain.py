#!/Users/sma/anaconda/bin/python
"""
Functionality:
* Build an array (dict) with element symbol, coordinates, and bond length
* Function that refreshes lengths and prints table
* Allow user to select desired bond, then how much to grow or shrink
* Saves new coordinates
* Saves config file with all changes saved for easy repeatability
"""

import sys, os
import math

INFILE = sys.argv[1]
ELEMENTS = sys.argv[2]

def distance(atom1, atom2):
    """ calculates distance between two points """
    dist = (((float(atom2[0]) - float(atom1[0])) ** 2) +
            ((float(atom2[1]) - float(atom1[1])) ** 2) +
            ((float(atom2[2]) - float(atom1[2])) ** 2)) ** 0.5
    return dist

def parse_elements():
    """ parses element symbol names from command line """
    elem = []
    for name in ELEMENTS.split():
        symbol = name.split('*')
        for x in xrange(int(symbol[0])): # don't need this motherfucker
            elem.append(symbol[1])
    return elem

def bond_lengths():
    """ calculates all bond lengths between atoms """
    bonds = []
    for idx, atom in enumerate(ATOMS): # why do we need atom here? or enumerate atoms for that matter
        if idx != 0:
            length = distance(ATOMS[idx-1], ATOMS[idx])
            bonds.append(length)
    return bonds

def print_table():
    """ prints table with atoms, symbol names, and bond lengths """
    os.system('clear')
    print "Available atoms and bond lengths (Bohrs):"
    for idx, atom in enumerate(ATOMS):
        print "{0} {1:18.14f}   {2:18.14f}    {3:18.14f}"\
            .format(SYMBOL[idx].ljust(2),
                    float(atom[0]),
                    float(atom[1]),
                    float(atom[2]))
        if idx != len(ATOMS) - 1:
            print '|' + 21 * 3 * '=' + "| {0}: {1:15.14f}"\
                .format(str(idx+1).zfill(2), BONDS[idx])
    print

ATOMS = [line.strip().split() for line in open(INFILE)]
SYMBOL = parse_elements()
BONDS = bond_lengths()
NUMBER = len(ATOMS)

print_table()
#SEL_BOND = input("Which bond do you want to modify? ")

#percent = 1 + DELTA/100.0
#hnew = distance(Si1, Si2) * percent
#zdist = Si1[2] - Si2[2]
#xdist = Si1[0] - Si2[0]
#ang = math.atan(zdist/xdist) + math.pi
#znew = hnew * math.sin(ang) + Si2[2]
#xnew = hnew * math.cos(ang) + Si2[0]
