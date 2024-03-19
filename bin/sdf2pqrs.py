#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# Tsuda laboratory, The University of Tokyo,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# transform _ONE_ molecule in an sdf file into a pseudo PQR file format
# named .pqrs (s for chemical element Symbol)
# Output line format:
# ---
# ^NUM_ATOMS:MOL_NAME$
# ^X Y Z Q R Chemical_symbol$
# ...
# ---

import argparse, math, os, sys, time
from rdkit import Chem
from rdkit.Chem import AllChem, rdForceFieldHelpers

ptable = Chem.GetPeriodicTable()

# better than default readline() which never throws an exception
def read_line_EOF(input):
    line = input.readline()
    if line == "":
        raise EOFError
    else:
        return line

# list all molecule names
def names_of_sdf_file(input_fn):
    try:
        with open(input_fn, 'r') as input:
            fst_name = read_line_EOF(input).strip()
            yield fst_name
            while True:
                line = read_line_EOF(input).strip()
                while line != "$$$$":
                    line = read_line_EOF(input).strip()
                next_name = read_line_EOF(input).strip()
                yield next_name
    except EOFError:
        pass

# output in a "pseudo PQR" file format
def mol2pqr(out, mol, name):
    conf0 = mol.GetConformer(0) # WARNING: all other conformers are ignored
    n = conf0.GetNumAtoms()
    print("%d:%s" % (n, name), file=out)
    props = AllChem.MMFFGetMoleculeProperties(mol)
    for a_i in mol.GetAtoms():
        i = a_i.GetIdx()        
        xyz = conf0.GetAtomPosition(i)
        # !!! WE SHOULD BE USING QEq PARTIAL CHARGES WITH UFF !!!
        q = props.GetMMFFPartialCharge(i)
        anum = a_i.GetAtomicNum()        
        r = ptable.GetRvdw(anum)
        # #UFF vdW parameters
        # x_i, D_i = rdForceFieldHelpers.GetUFFVdWParams(mol, i, i)
        elt = a_i.GetSymbol()
        print("%g %g %g %g %g %s" %
              (xyz.x, xyz.y, xyz.z, q, r, elt), file=out)

# FBR:  - output the torsion tree
#       - detect rotatable bonds
#       - for each RotBond: - create group for all atoms on the left
#                           - group for all atoms on the right
#                           - make smaller group the right one

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "sdf to pqr conversion, with MMFF94 partial charges")
    parser.add_argument("-i", metavar = "input.sdf", dest = "input_fn",
                        help = "one 3D conformer input file")
    parser.add_argument("-o", metavar = "output.pqrs", dest = "output_fn",
                        help = "output file")
    # parse CLI ---------------------------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output_fn = args.output_fn
    mol_names = names_of_sdf_file(input_fn)
    mol_supplier = Chem.SDMolSupplier(input_fn, sanitize=True, removeHs=False)
    # parse CLI end -----------------------------------------------------------
    count = 0
    errors = 0
    with open(output_fn, 'w') as out:
        for mol, name in zip(mol_supplier, mol_names):
            # print("%d atoms" % mol.GetNumHeavyAtoms(), file=sys.stderr)
            if mol == None:
                errors += 1
            else:
                mol2pqr(out, mol, name)
            count += 1
    after = time.time()
    dt = after - before
    print("%d molecules @ %.2fHz; %d errors" % (count, count / dt, errors),
          file=sys.stderr)
