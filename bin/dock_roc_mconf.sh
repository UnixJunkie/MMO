#!/bin/bash
#
# only keep lowest scoring conformer from a multi-conformer
# name-scores file

set -u
#set -x # DEBUG

IN=$1
OUT=`mktemp`

# sort by incr. docking score; remove conf_id suffix to mol_name
sort -n -k2 $IN | sed 's/_1.*\t/\t/g' > $OUT
# keep only lowest scoring conformer to stdout
molenc_uniq -i $OUT -f 1

# cleanup
rm -f $OUT
