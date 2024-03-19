#!/usr/bin/env python3

import argparse, mol2, re, sys, time
from rdkit import Chem
from rdkit.Chem.rdchem import RWMol, BondType
from rdkit.Chem.rdmolops import (AddHs,
                                 RemoveHs,
                                 AssignAtomChiralTagsFromStructure,
                                 SanitizeMol)

# aliases
sp2 = Chem.rdchem.HybridizationType.SP2
sp3 = Chem.rdchem.HybridizationType.SP3
single_bond = Chem.rdchem.BondType.SINGLE

# single bonds out of rings; each end being (sp2 or sp3)
def list_rotatable_bonds(mol):
    res = []
    for b in mol.GetBonds():
        if b.GetBondType() == single_bond and not b.IsInRing():
            j_a = b.GetBeginAtom()
            k_a = b.GetEndAtom()
            j_h = j_a.GetHybridization()
            k_h = k_a.GetHybridization()
            if (j_h == sp2 or j_h == sp3) and \
               (k_h == sp2 or k_h == sp3):
                jk = (j_a.GetIdx(), k_a.GetIdx())
                res.append(jk)
    return res

# atoms in smallest group will be rotatable; others are fixed for this RotBond
def smallest_group(mol):
    d = {}
    for a in mol.GetAtoms():
        g = a.GetIntProp('group')
        if g in d:
            d[g] += 1
        else:
            d[g] = 1
    # one RotBond makes two groups
    assert(len(d.keys()) == 2)
    smallest_group = -1
    smallest_size = sys.maxsize # ~= integer infinity
    for g in d.keys():
        group_size = d[g]
        if group_size < smallest_size:
            smallest_size = group_size
            smallest_group = g
    return smallest_group

#create atom groups; each end of a RotBond belongs to a different group
def atom_rotatable_groups(mol, rot_bonds):
    # remove RotBond from the molecular graph, so that we don't walk over it
    for j, k in rot_bonds:
        rw = Chem.RWMol(mol)
        rw.RemoveBond(j, k)
        print('%d-%d=' % (j, k), end='')
        # initially, each atom in a different group
        for a in rw.GetAtoms():
            i = a.GetIdx()
            a.SetIntProp('group', i)
        # loop until groups are stable
        changed = True
        while changed:
            changed = False
            for b in rw.GetBonds():
                a_i = b.GetBeginAtom()
                a_j = b.GetEndAtom()
                g_i = a_i.GetIntProp('group')
                g_j = a_j.GetIntProp('group')
                if g_i < g_j:
                    a_j.SetIntProp('group', g_i)
                    changed = True
                elif g_j < g_i:
                    a_i.SetIntProp('group', g_j)
                    changed = True
        movable = smallest_group(rw)
        # print them
        for a in rw.GetAtoms():
            if a.GetIntProp('group') == movable:
                print(' 1', end='')
            else:
                print(' 0', end='')
        print()

ptable = Chem.GetPeriodicTable()

def print_dist_matrix(mol):
    mtx = Chem.GetDistanceMatrix(mol)
    num_atoms = mol.GetNumAtoms() # checked below
    assert(num_atoms == len(mtx)) # be extra cautious
    for i in range(num_atoms):
        for j in range(num_atoms):
            if j == 0:
                print('%d' % mtx[i][j], end='')
            else:
                print(' %d' % mtx[i][j], end='')
        print()

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "read in a mol2 file _WITH_ partial charges")
    parser.add_argument("-i", metavar = "input.mol2", dest = "input_fn",
                        help = "one 3D conformer input file")
    # parse CLI ---------------------------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    mol_supplier = mol2.Mol2MolSupplier(input_fn, sanitize=True, removeHs=False)
    # parse CLI end -----------------------------------------------------------
    count = 0
    errors = 0
    for mol in mol_supplier:
        name = ""
        if mol == None:
            errors += 1
            continue
        else:
            count += 1
            name = mol.GetProp("name")
        if mol2.all_null_charges(mol):
            print("no partial charges in %s" % name, file=sys.stderr)
            continue
        # print stuffs about this molecule
        conf0 = mol.GetConformer(0)
        num_atoms = conf0.GetNumAtoms()
        rot_bonds = list_rotatable_bonds(mol)
        num_rot_bonds = len(rot_bonds)
        print('%d:%d:%s' % (num_atoms, num_rot_bonds, name))
        # pqrs
        for a_i in mol.GetAtoms():
            i = a_i.GetIdx()
            xyz = conf0.GetAtomPosition(i)
            q = a_i.GetDoubleProp('q')
            anum = a_i.GetAtomicNum()
            r = ptable.GetRvdw(anum)
            elt = a_i.GetSymbol()
            print("%g %g %g %g %g %s" % (xyz.x, xyz.y, xyz.z, q, r, elt))
        # # rotBonds
        # for j, k in rot_bonds:
        #     print('%d %d' % (j, k))
        # atom rot. groups
        atom_rotatable_groups(mol, rot_bonds)
        # needed to find 1,4 interactions
        # FBR: should be all atoms but will probably skip hydrogens
        print_dist_matrix(mol)
    after = time.time()
    dt = after - before
    print("%d molecules @ %.2fHz; %d errors" % (count, count / dt, errors),
          file=sys.stderr)
