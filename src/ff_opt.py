#!/usr/bin/env python3

# generate low energy conformers for given FF

import argparse, mol2, rdkit, sys, time
from rdkit import Chem
from rdkit.Chem import AllChem

# 3D conf. gen
# m = Chem.MolFromSmiles('C1CCC1OC')
# m2=Chem.AddHs(m)
# AllChem.EmbedMolecule(m2)
# res = AllChem.MMFFOptimizeMoleculeConfs(m2)

# 3D confs gen
# m = Chem.MolFromSmiles('C1CCC1OC')
# m2=Chem.AddHs(m)
# # run ETKDG
# cids = AllChem.EmbedMultipleConfs(m2, numConfs=50)
# rmslist = []
# AllChem.AlignMolConformers(m2, RMSlist=rmslist)
# #rms = AllChem.GetConformerRMS(m2, 1, 9, prealigned=True)
# res = AllChem.MMFFOptimizeMoleculeConfs(m2)

# # get FF enregy w/o any minimization
# # UFF
# ffu = AllChem.UFFGetMoleculeForceField(mol)
# energy_value_U = ffu.CalcEnergy()
# print(energy_value_U)
# # MMFF
# mol = rdmolops.AddHs(mol, addCoords = True)
# mp = AllChem.MMFFGetMoleculeProperties(mol)
# ffm = AllChem.MMFFGetMoleculeForceField(mol, mp)
# energy_value_M = ffm.CalcEnergy()
# print(energy_value_M)

# # single molecule minimization, without check of convergence
# ff = AllChem.UFFGetMoleculeForceField(mol, confId = 0)
# # print("E before: %f" % ff.CalcEnergy())
# ff.Minimize()
# energy = ff.CalcEnergy()
# # print("E after: %f" % energy)

def minimize_conformer(ff, mol):
    # ligand is supposed to be already properly protonated
    # for given pH and in 3D
    assert(mol.GetNumConformers() == 1)
    conv_ene = []
    if ff == "uff":
        conv_ene = AllChem.UFFOptimizeMoleculeConfs(mol)
    elif ff == "mmff":
        conv_ene = AllChem.MMFFOptimizeMoleculeConfs(mol)
    else:
        print("minimize_conformer: unsupported FF: %s" % ff,
              file=sys.stderr)
        assert(False)
    not_converged, ene = conv_ene[0]
    assert(not_converged == 0)
    return (mol, ene)

# dump conf0 in .xyz format
def xyz_dump(mol, ene):
    name = mol.GetProp("name")
    conf0 = mol.GetConformer(0)
    n = conf0.GetNumAtoms() # ask num atoms to conformer, NOT to molecule
    print(n)
    print('%s:%g' % (name, ene))
    for i in range(n):
        a_i = mol.GetAtomWithIdx(i)
        xyz = conf0.GetAtomPosition(i)
        elt = a_i.GetSymbol()
        print('%s %g %g %g' % (elt, xyz.x, xyz.y, xyz.z))

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "read in a mol2 file _WITH_ partial charges; compute E_intra")
    parser.add_argument("-i", metavar = "input.mol2", dest = "input_fn",
                        help = "one 3D conformer input file")
    parser.add_argument("-uff", dest = "uff", action='store_true',
                        help = "use UFF force field")
    parser.add_argument("-mmff", dest = "mmff", action='store_true',
                        help = "use MMFF94 force field")
    # parse CLI ---------------------------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    use_uff = args.uff
    use_mmff = args.mmff
    assert(use_uff != use_mmff)
    mol_supplier = mol2.Mol2MolSupplier(input_fn, sanitize=True, removeHs=False)
    # parse CLI end -----------------------------------------------------------
    ff = "uff"
    if use_mmff:
        ff = "mmff"
    count = 0
    errors = 0
    for mol in mol_supplier:
        if mol2.all_null_charges(mol):
            print("no partial charges in %s" % mol.GetProp("name"),
                  file=sys.stderr)
            mol = None
        if mol == None:
            errors += 1
        else:
            count += 1
        mol_w_min_conf0, ene = minimize_conformer(ff, mol)
        xyz_dump(mol_w_min_conf0, ene)
    after = time.time()
    dt = after - before
    print("%d molecules @ %.2fHz; %d errors" % (count, count / dt, errors),
          file=sys.stderr)
