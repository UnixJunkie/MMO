#!/bin/bash
#
# ROC curve from a MMO dock docking_scores.tsv file;
# actives are supposed to have a molecule name prefixed with "active"
# so that they can be labeled

#set -x # DEBUG
set -u

INPUT_FN=$1
SCORE_LABELS_FN=${INPUT_FN}.score_labels
ROC_FN=${INPUT_FN}.roc
AUC_FN=${INPUT_FN}.auc

echo ${SCORE_LABELS_FN}

# warn about strange values in input file
echo "INF values:" `grep -c -w inf ${INPUT_FN}`
echo "NaN values:" `grep -c -w nan ${INPUT_FN}`

# ignore inf and nan values from input file
grep -v -w inf ${INPUT_FN} | grep -v -w nan | \
    awk '/^active/{printf("%f 1\n", -$2)}!/^active/{printf("%f 0\n", -$2)}' | sort -n -r -k1 > ${SCORE_LABELS_FN}

# ~/src/UNIX-fun/bin/roc_auc -i ${SCORE_LABELS_FN}

croc-curve < ${SCORE_LABELS_FN} > ${ROC_FN} 2> ${AUC_FN}
AUC=`awk '{printf("%.2f\n",$5)}' ${AUC_FN}`
TITLE=`echo $INPUT_FN AUC=${AUC} | sed 's/_/\\\_/g'`

gnuplot --persist <<EOF
set tics out nomirror
set size square
set key out center top
set xlabel 'FPR'
set ylabel 'TPR'
f(x) = x
plot f(x) not, '${ROC_FN}' u 1:2 w l t '${TITLE}'
EOF
