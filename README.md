# MMO: Molecular Mechanics in OCaml

Some code for flexible ligand / rigid protein calculations.

Warning: this code is prototypical in nature.
Maybe, part of it will be released as a proper standalone library in the long-term
(the Mol and Mol2 modules plus their dependencies).

# Dataset

The dataset that was produced using this software can be downloaded from:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10682034.svg)](https://doi.org/10.5281/zenodo.10682034)

## To compile the software

Tested on Ubuntu Linux 24.04 w/ ocaml-4.14.1.

On Linux or a UNIX-like:
```bash
mkdir -p ~/src
cd src
git clone https://github.com/UnixJunkie/MMO.git
cd MMO
sudo apt install opam
opam init
eval `opam env --shell=bash`
opam pin add mmo .
make
```
## How to score a ligand conformer w/ psi4 (true QM calculation)

This requires mayachemtools and psi4 to be installed.

```bash
~/src/MMO/bin/psi4_ene.sh ligand.sdf
```

## How to score a ligand conformer w/ ANI-2 (approximated QM)

This requires the Python package torchani to be installed:

```bash
~/src/MMO/bin/QM_score_mol2.sh ligand.mol2
```

## How to run a MC simulation in dihedral space to minimize a docked ligand using ANI-2?

If you have compiled the software and you have torchani installed, you can run:
```bash
./lds -f --no-interp --ligand-only --intra-QM -lig docked_ligand.mol2 \
  -rec prot_rec.pqrs -roi prot_ROI.bild -steps 10k --ene 50 --no-compress --out-dir MC10kQM
```
