#!/bin/bash
set -xe
TXT_CONFIGS="IDX_NN IDX_NN_QUANTILE IDX_NN_QUANTILE_LG IDX_NN_QUANTILE_LG_16"
INT_CONFIGS="IDX_NN_INT IDX_NN_DOCID_SMART_INT IDX_NN_K3_DAAT_INT"

scripts/build_config.sh $TXT_CONFIGS #$INT_CONFIGS
scripts/build.sh gen_patterns

coll=collections/ENWIKISML
if [[ $# > 0 ]]; then
    coll="$1"
fi

scripts/compare.py -c "$coll" $TXT_CONFIGS \
    -b build/release -q 3000 -r 1 --no_multi_occ
#scripts/compare.py -c "$coll" $INT_CONFIGS \
    #-b build/release -q 3000 -r 1
