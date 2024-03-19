#!/bin/bash
#
# QM score conformer from a .mol2 file using torchani

set -u

# where MMCLD is installed
MMCLD_DIR=~/src/MMCLD

XYZ=`mktemp --suffix ".xyz"`
QM_ENE=${XYZ}.qm_ene

# mol2 -> xyz
${MMCLD_DIR}/mol2pqrs -xyz -lig $1 -o ${XYZ}

# QM score it
${MMCLD_DIR}/src/QM_energy.py ${XYZ} ${QM_ENE} 2>/dev/null 1>/dev/null

# display result
cat ${QM_ENE}

# cleanup
rm -f ${XYZ} ${QM_ENE}
