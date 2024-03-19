#!/usr/bin/env python3

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
    count = 0
    errors = 0
    for i, mol in enumerate(mol_supplier):
        if mol == None:
            errors += 1
        else:
            if mol2.all_null_charges(mol):
                print("no partial charges in mol %d" % i, file=sys.stderr)
                errors += 1
            else:
                name = mol.GetProp("name")
                count += 1
                if use_mmff:
                    if AllChem.MMFFHasAllMoleculeParams(mol):
                        props = AllChem.MMFFGetMoleculeProperties(mol)
                        mmff = AllChem.MMFFGetMoleculeForceField(mol, props)
                        mmff_ene = mmff.CalcEnergy()
                        print('%s: %f MMFF' % (name, mmff_ene))
                    else:
                        errors += 1
                elif use_uff:
                    if AllChem.UFFHasAllMoleculeParams(mol):
                        uff = AllChem.UFFGetMoleculeForceField(mol)
                        uff_ene = uff.CalcEnergy()
                        print('%s: %f UFF' % (name, uff_ene))
                    else:
                        errors += 1
    after = time.time()
    dt = after - before
    print("%d molecules @ %.2fHz; %d errors" % (count, count / dt, errors),
          file=sys.stderr)
