#!/usr/bin/env python3

# compute geometric center of a .pqrs file
# both ligands and proteins are supported
# proteins: only (xyz, charge, symbol) lines
# ligands: same but slightly different header line
#          and additional info after coordinates

import sys

fn = sys.argv[1]

i = 0
num_atoms = 0
sum_x = 0.0
sum_y = 0.0
sum_z = 0.0

def parse_num_atoms(line):
    return int(line.split(':')[0])

to_read = open(fn, 'r')
lines = list(to_read.readlines())
to_read.close()

header = lines[0].strip()
num_atoms = parse_num_atoms(header)
print(header)
# skip header
i += 1

while i <= num_atoms:
    split = lines[i].strip().split(' ')
    x = float(split[0])
    y = float(split[1])
    z = float(split[2])
    sum_x += x
    sum_y += y
    sum_z += z
    i += 1

print('%g %g %g' % (sum_x / num_atoms, \
                    sum_y / num_atoms, \
                    sum_z / num_atoms))
