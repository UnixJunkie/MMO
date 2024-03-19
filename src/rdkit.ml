module Rdkit : sig
  type t

  val of_pyobject : Pytypes.pyobject -> t
  val to_pyobject : t -> Pytypes.pyobject
  val __init__ : mol2:string -> i:int -> unit -> t
  val is_valid_UFF : t -> unit -> bool
  val is_valid_MMFF : t -> unit -> bool

  val update_coords :
    t -> xs:float array -> ys:float array -> zs:float array -> unit -> unit

  val score_UFF : t -> unit -> float
  val score_MMFF : t -> unit -> float
end = struct
  let filter_opt l = List.filter_map Fun.id l

  let py_module =
    lazy
      (let source =
         {pyml_bindgen_string_literal|# Copyright (C) 2001-2008 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ANY CORRECTION IN THIS FILE NEEDS TO BE REPLICATED INTO rdkit_uff.py
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import re, sys
from rdkit import Chem
from rdkit.Chem.rdchem import RWMol, Conformer, Atom, BondType
from rdkit.Chem.rdmolfiles import MolFromMol2Block
from rdkit.Chem.rdmolops import (AddHs,
                                 RemoveHs,
                                 AssignAtomChiralTagsFromStructure,
                                 SanitizeMol)
from rdkit.Geometry import Point3D

def all_null_charges(mol):
    res = True
    for a in mol.GetAtoms():
        if a.GetDoubleProp('q') != 0.0:
            res = False
    return res

def symbol_of_mol2_type(mol2_type):
    try:
        i = mol2_type.index('.')
        return mol2_type[0:i]
    except ValueError:
        # only chemical element string; no hybridization suffix
        return mol2_type

def MolFromCommonMol2Block(block, sanitize=True, removeHs=False):
    """Patch MolFromMol2Block to be more flexible and parse mol2 like SD files,
    thus alowing to read/write mol2 files without infering the non-canonical
    atom types used by other software.
    """
    mol = RWMol()
    expected_num_atoms = 0
    expected_num_bonds = 0
    actual_num_atoms = 0
    actual_num_bonds = 0
    coords = []
    mode = 0 # 0 - meta, 1 - atoms, 2 - bonds, 3 - exit
    for n, line in enumerate(block.split('\n')):
        stripped = line.strip()
        if stripped == '':
            continue

        if line[:1] == '@':
            rline = line.rstrip()
            if rline == '@<TRIPOS>MOLECULE':
                mode = 0
            elif rline == '@<TRIPOS>ATOM':
                mode = 1
            elif rline == '@<TRIPOS>BOND':
                mode = 2
            else:
                mode = 3
                #break # ???
            continue

        # 0. Get molecule meta-data
        if mode == 0:
            if n == 1:
                mol.SetProp('name', stripped)
            elif n == 2:
                nums = stripped.split()
                expected_num_atoms = int(nums[0])
                expected_num_bonds = int(nums[1])
            # n = 3: SMALL/PROTEIN
            # n = 4: GASTEIGER/USER_CHARGES
            elif n > 5:
                raise ValueError('Too many lines in @<TRIPOS>MOLECULE block')
        # 1. Add atoms
        elif mode == 1:
            data = re.split('\s+', stripped)
            idx = int(data[0]) - 1
            #atom_name = data[1]
            x = float(data[2])
            y = float(data[3])
            z = float(data[4])
            coords.append(Point3D(x, y, z))
            symbol = symbol_of_mol2_type(data[5])
            atom = Atom(symbol)
            charge = float(data[8])
            new_idx = mol.AddAtom(atom)
            assert(new_idx == idx)
            # atom properties must be set on the atom
            # AFTER it has been added to the molecule !!!
            atom_i = mol.GetAtomWithIdx(new_idx)
            atom_i.SetDoubleProp('q', charge)
            if data[5] == "N.4":
                atom_i.SetFormalCharge(1)
            actual_num_atoms += 1
        # 2. Add bonds
        elif mode == 2:
            data = re.split('\s+', stripped)
            idx = int(data[0]) - 1
            begin_atom = int(data[1]) - 1
            end_atom = int(data[2]) - 1
            if data[3] == 'ar':
                order = BondType.AROMATIC
            elif data[3] == 'am':
                order = BondType.SINGLE
            else:
                order = BondType.values[int(data[3])]
            mol.AddBond(begin_atom, end_atom, order)
            actual_num_bonds += 1

    # check we got a complete MOL2 block
    if actual_num_atoms != expected_num_atoms:
        # MOL2 is an all atom format; hydrogens _included_
        print('mol2.MolFromCommonMol2Block: incorrect number of atoms for %s' %
              mol.GetProp('name'), file=sys.stderr)
        return None
    if actual_num_bonds != expected_num_bonds:
        print('mol2.MolFromCommonMol2Block: incorrect number of bonds for %s' %
              mol.GetProp('name'), file=sys.stderr)
        return None

    # 3. Remove Hs, sanitize

    # convert to ROMol
    mol = mol.GetMol()

    try:
        for atom in mol.GetAtoms():
            atom.UpdatePropertyCache()
    except:
        # FBR: molecules with valence(N) == 4 could become N+ and be corrected !
        return None
    # There is no such function, just marking it TODO
    #Chem.DetectAtomStereoChemistry(mol, conf)
    AssignAtomChiralTagsFromStructure(mol)

    # create the 3D conformer
    conf = Conformer(len(coords))
    for i, xyz in enumerate(coords):
        conf.SetAtomPosition(i, xyz)

    cid = mol.AddConformer(conf)
    assert(cid == 0)

    if sanitize:
        try:
            if removeHs:
                mol = RemoveHs(mol)
            else:
                SanitizeMol(mol)
        except:
            return None
    # Chem.DetectBondStereoChemistry(mol, conf)

    return mol

class Mol2MolSupplier:
    def __init__(self, filename, *args, **kwargs):
        """Reads a multi-mol Mol2 file
          ARGUMENTS:

            - filename: the file to read  or file-like object
            - args, kwargs: arbitrary arguments to pass to internal MolFromMol2Block

          RETURNS:

            None
        """
        self.f = filename
        self._args = args
        self._kwargs = kwargs

    def __iter__(self):
        """ Iterates over molecules in file """
        block = ''
        data = ''
        n = 0
        if hasattr(self.f, 'read') and hasattr(self.f, 'close'):
            f = self.f
        else:
            f = open(self.f)
        for line in f:
            if line[:1] == '#':
                data += line
            elif line[:17] == '@<TRIPOS>MOLECULE':
                if n > 0:  #skip `zero` molecule (any preciding comments and spaces)
                    yield MolFromCommonMol2Block(block, *self._args, **self._kwargs)
                n += 1
                block = data
                data = ''
            block += line
        # open last molecule
        if block:
            yield MolFromCommonMol2Block(block, *self._args, **self._kwargs)
        f.close()
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
|pyml_bindgen_string_literal}
       in
       let filename =
         {pyml_bindgen_string_literal|rdkit_uff.py|pyml_bindgen_string_literal}
       in
       let bytecode = Py.compile ~filename ~source `Exec in
       Py.Import.exec_code_module
         {pyml_bindgen_string_literal|rdkit_uff|pyml_bindgen_string_literal}
         bytecode)

  let import_module () = Lazy.force py_module

  type t = Pytypes.pyobject

  let of_pyobject pyo = pyo
  let to_pyobject x = x

  let __init__ ~mol2 ~i () =
    let callable = Py.Module.get (import_module ()) "Rdkit" in
    let kwargs =
      filter_opt
        [ Some ("mol2", Py.String.of_string mol2); Some ("i", Py.Int.of_int i) ]
    in
    of_pyobject @@ Py.Callable.to_function_with_keywords callable [||] kwargs

  let is_valid_UFF t () =
    let callable = Py.Object.find_attr_string t "is_valid_UFF" in
    let kwargs = filter_opt [] in
    Py.Bool.to_bool
    @@ Py.Callable.to_function_with_keywords callable [||] kwargs

  let is_valid_MMFF t () =
    let callable = Py.Object.find_attr_string t "is_valid_MMFF" in
    let kwargs = filter_opt [] in
    Py.Bool.to_bool
    @@ Py.Callable.to_function_with_keywords callable [||] kwargs

  let update_coords t ~xs ~ys ~zs () =
    let callable = Py.Object.find_attr_string t "update_coords" in
    let kwargs =
      filter_opt
        [
          Some ("xs", Py.List.of_array_map Py.Float.of_float xs);
          Some ("ys", Py.List.of_array_map Py.Float.of_float ys);
          Some ("zs", Py.List.of_array_map Py.Float.of_float zs);
        ]
    in
    ignore @@ Py.Callable.to_function_with_keywords callable [||] kwargs

  let score_UFF t () =
    let callable = Py.Object.find_attr_string t "score_UFF" in
    let kwargs = filter_opt [] in
    Py.Float.to_float
    @@ Py.Callable.to_function_with_keywords callable [||] kwargs

  let score_MMFF t () =
    let callable = Py.Object.find_attr_string t "score_MMFF" in
    let kwargs = filter_opt [] in
    Py.Float.to_float
    @@ Py.Callable.to_function_with_keywords callable [||] kwargs
end
