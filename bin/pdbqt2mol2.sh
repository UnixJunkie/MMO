#!/bin/bash

# Assign MMFF94 partial charges and hydrogens (protonate) at given pH (7.4) usinb openbabel
obabel $1 -O $2 --partialcharge mmff94 -p 7.4
