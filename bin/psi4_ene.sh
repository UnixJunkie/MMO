#!/bin/bash

set -u
#set -x

# make everything happen in this temp. dir
TMP=`mktemp -d`
IN=${TMP}"/in.sdf"
OUT=${TMP}"/out.sdf"

cp $1 ${IN}
cd ${TMP}
~/usr/mayachemtools/bin/Psi4CalculateEnergy.py --quiet yes -i ${IN} -o ${OUT}
name=`head -1 ${IN}`
ene=`grep -A1 Psi4_Energy ${OUT} | tail -1`

echo $1"	"${name}"	"${ene}

rm -rf ${TMP}
