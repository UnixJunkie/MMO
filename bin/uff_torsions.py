#!/usr/bin/env python3

# list all rotatable bonds

# for each single bond not in a ring, classify bond as
# - sp3-sp2sp2 or sp2sp2-sp3
# - sp3-sp3
# - sp3-sp2 or sp2-sp3

import rdkit, scanf, sys
from rdkit import Chem
from rdkit.Chem import RWMol
from enum import Enum

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
                res.append(b)
    return res

# get all i-j-k-l quads for this bond
def list_quads(bond):
    j_a = b.GetBeginAtom()
    k_a = b.GetEndAtom()
    j = j_a.GetIdx()
    k = k_a.GetIdx()
    iset = set()
    lset = set()
    for i_a in j_a.GetNeighbors():
        iset.add(i_a.GetIdx())
    for l_a in k_a.GetNeighbors():
        lset.add(l_a.GetIdx())
    iset.remove(k)
    lset.remove(j)
    # ignore 3-membered rings
    inter = iset.intersection(lset)
    icands = iset.difference(inter)
    lcands = lset.difference(inter)
    res = []
    # enumerate quads
    for i in icands:
        for l in lcands:
            res.append((i,j,k,l))
    return res

class UFF_RotBond(Enum):
    SP3_SP3    = 0
    SP3_SP2SP2 = 1
    SP3_SP2    = 2

def classify_quads(mol, quads):
    res = []
    for (i, j, k, l) in quads:
        a_j = mol.GetAtomWithIdx(j)
        a_k = mol.GetAtomWithIdx(k)
        j_h = a_j.GetHybridization()
        k_h = a_k.GetHybridization()
        if j_h == sp3 and k_h == sp3:
            res.append(UFF_RotBond.SP3_SP3)
        else:
            a_i = mol.GetAtomWithIdx(i)
            a_l = mol.GetAtomWithIdx(l)
            i_h = a_i.GetHybridization()
            l_h = a_l.GetHybridization()
            if (i_h == sp2 and j_h == sp2 and k_h == sp3) or \
               (j_h == sp3 and k_h == sp2 and l_h == sp2):
                res.append(UFF_RotBond.SP3_SP2SP2)
            elif (j_h == sp2 and k_h == sp3) or \
                 (j_h == sp3 and k_h == sp2):
                res.append(UFF_RotBond.SP3_SP2)
            else:
                assert(false)
    return res

input_fn = sys.argv[1]
mol_supplier = Chem.SDMolSupplier(input_fn, sanitize=True, removeHs=False)

for mol in mol_supplier:
    bonds = list_rotatable_bonds(mol)
    print("bonds: %d" % len(bonds))
    for b in bonds:
        quads = list_quads(b)
        for q in quads:
            print(q)
        for c in classify_quads(mol, quads):
            print(c)
