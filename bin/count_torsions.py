#!/usr/bin/env python3

# count number of torsion angles in each line of a docked_Torsion_Strain.csv file

import re, sys

input = open(sys.argv[1])

for line in input.readlines():
    name=line.split('|')[0]
    # each torsion angle identifies four consecutive atoms in the molecule: [i, j, k, l]
    print("%s\t%d" % (name, len(re.findall('\[[0-9]+, [0-9]+, [0-9]+, [0-9]+\]', line))))
