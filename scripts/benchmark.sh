#!/bin/bash
set -xe
TXT_CONFIGS="IDX_NN IDX_NN_DOCID_SMART IDX_NN_K3_DAAT"
INT_CONFIGS="IDX_NN_INT IDX_NN_DOCID_SMART_INT IDX_NN_K3_DAAT_INT"

scripts/build_config.sh $TXT_CONFIGS $INT_CONFIGS
scripts/build.sh gen_patterns

coll=collections/TEST_TXT
if [[ $# > 0 ]]; then
    coll="$1"
fi

scripts/compare.py -c "$coll" $TXT_CONFIGS \
    -b build/release -q 1000 -r 1 $@
scripts/compare.py -c "$coll" $INT_CONFIGS \
    -b build/release -q 1000 -r 1
