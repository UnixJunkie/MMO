# MMO
Molecular Mechanics in OCaml

A bunch of code for flexible ligand / rigid protein calculations.

Warning: this code is prototypical in nature.
Maybe, part of it will be released as a proper standalone library in the mid-term future
(the Mol and Mol2 modules plus their dependencies).

## to compile the software

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
