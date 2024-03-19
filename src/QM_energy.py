#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# Tsuda laboratory, The University of Tokyo,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# compute ligand internal energy using a DNN approximation of a QM FF

import sys, torch, torchani

# elements supported by ANI2x
supported_atoms = set(["C","H","N","O","S","F","Cl"])

# elements supported by our UFF implementation
elt2anum = {}
elt2anum["C"]  =  6
elt2anum["H"]  =  1
elt2anum["N"]  =  7
elt2anum["O"]  =  8
elt2anum["P"]  = 15
elt2anum["S"]  = 16
elt2anum["F"]  =  9
elt2anum["Cl"] = 17
elt2anum["Br"] = 35
elt2anum["I"]  = 53

class UnsupportedAtom(Exception):
    pass

class EOF(Exception):
    pass

def readline(f):
    l = f.readline()
    if l == '':
        raise EOF
    return l

# retrieve a whole molecule .xyz block from a read-only opened file
def read_xyz_block(f):
    lines = []
    line0 = readline(f).strip()
    lines.append(line0)
    num_atoms = int(line0)
    comment = readline(f).strip()
    lines.append(comment)
    name = comment.split('\t')[1]
    print("I: num_atoms: %d" % num_atoms, file=sys.stderr)
    for i in range(num_atoms):
        lines.append(readline(f).strip())
    return lines

def parse_xyz_block(lines):
    coords = [] # atom 3D coords
    species = [] # chemical elements
    n = len(lines)
    num_atoms = int(lines[0])
    comment = lines[1]
    name = comment.split('\t')[1]
    for i in range(2, n):
        elt, x, y, z = lines[i].split(' ')
        if elt not in supported_atoms:
            print("E: %s: unsupported atom: %s" % (name, elt),
                  file=sys.stderr)
            print('%s\tnan' % name)
            raise UnsupportedAtom # /!\ EXCEPTION /!\
        xyz = [float(x), float(y), float(z)]
        anum = elt2anum[elt]
        coords.append(xyz)
        species.append(anum)
    return (name, [coords], [species])

cpu_or_gpu = torch.device('cpu') # FBR: faster on my computer
# cpu_or_gpu = None
# if torch.cuda.is_available():
#     print("I: running on GPU", file=sys.stderr)
#     cpu_or_gpu = torch.device('cuda')
# else:
#     print("W: running on CPU", file=sys.stderr)
#     cpu_or_gpu = torch.device('cpu')
model = torchani.models.ANI2x(periodic_table_index=True).to(cpu_or_gpu)
hartree_to_kcal_per_mol = 627.5094740631

input_fn = sys.argv[1]
output_fn = sys.argv[2]

with open(input_fn, 'r') as input:
    with open(output_fn, 'w') as output:
        while True:
            try:
              xyz_block = read_xyz_block(input)
              mol_name, xyz, anums = parse_xyz_block(xyz_block)
              # coordinates(Angstrom) and chemical species(anum); energies in Hartree
              coordinates = torch.tensor(xyz, requires_grad=False, device=cpu_or_gpu)
              species = torch.tensor(anums, device=cpu_or_gpu)
              energy = hartree_to_kcal_per_mol * model((species, coordinates)).energies.item()
              # print result
              print('%s\t%f' % (mol_name, energy), file=output)
            except UnsupportedAtom:
                pass
            except EOF:
                break
