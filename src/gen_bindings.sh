#!/bin/bash

# embed all python files into one ---------------------------------------------
cat mol2.py uff.py > rdkit_uff.py

pyml_bindgen rdkit_uff_specs.txt rdkit_uff Rdkit \
             --embed-python-source rdkit_uff.py \
             --caml-module=Rdkit --of-pyo-ret-type=no_check > rdkit.ml

# format generated code
ocamlformat --inplace --enable-outside-detected-project rdkit.ml

# torchani QM ene wrapper -----------------------------------------------------
pyml_bindgen QM_ene_specs.txt QM_ene Torchani \
             --embed-python-source QM_ene.py \
             --caml-module=QM_ene --of-pyo-ret-type=no_check > QM_ene.ml

ocamlformat --inplace --enable-outside-detected-project QM_ene.ml
