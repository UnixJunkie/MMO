import rdkit, typing
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

def mol2_read_one_conf(input_fn, lig_i):
    # Mol2MolSupplier comes from mol2.py; which will be prepended to the current file
    mol_supplier = Mol2MolSupplier(input_fn, sanitize=True, removeHs=False)
    for i, mol in enumerate(mol_supplier):
        if i == lig_i:
            return mol

class Rdkit:
    # this is needed because the OCaml side needs to know how
    # to get an object of type t
    def __init__(self, mol2: str, i: int):
        self.mol = mol2_read_one_conf(mol2, i)

    # could rdkit read the molecule from a mol2 file (very often: NO ! )
    # and is this molecule OK w/ UFF
    # FBR: maybe I need a much more robust rdkit mol2 parser
    def is_valid_UFF(self):
        return (self.mol != None and
                AllChem.UFFHasAllMoleculeParams(self.mol))

    def is_valid_MMFF(self):
        return (self.mol != None and
                AllChem.MMFFHasAllMoleculeParams(self.mol))

    # allow conf0's coordinates to be updated
    def update_coords(self, xs, ys, zs):
        conf0 = self.mol.GetConformer(0)
        n = conf0.GetNumAtoms()
        assert(n == len(xs) == len(ys) == len(zs)) # FBR: maybe comment out later
        for i in range(n):
            conf0.SetAtomPosition(i, Point3D(xs[i], ys[i], zs[i]))

    # score conf0 using UFF
    def score_UFF(self) -> float:
        uff = AllChem.UFFGetMoleculeForceField(self.mol)
        return uff.CalcEnergy()

    def score_MMFF(self) -> float:
        mp = AllChem.MMFFGetMoleculeProperties(self.mol)
        mmff = AllChem.MMFFGetMoleculeForceField(self.mol, mp)
        return mmff.CalcEnergy()
