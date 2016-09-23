#!/bin/bash
set -xe
TXT_CONFIGS="BRUTE_TXT IDX_NN IDX_NN_DOCID_SMART IDX_NN_K3_DAAT"
INT_CONFIGS="BRUTE_INT IDX_NN_INT IDX_NN_DOCID_SMART_INT IDX_NN_K3_DAAT_INT"
INTERSECT_TXT_CONFIGS="BRUTE_TXT IDX_NN_K3_DAAT"
INTERSECT_INT_CONFIGS="BRUTE_INT"

test_txt() {
    coll="$1"
    scripts/build_config.sh -d $TXT_CONFIGS $INTERSECT_TXT_CONFIGS
    scripts/compare.py -c "$coll" $TXT_CONFIGS -b build/debug
    scripts/compare.py -c "$coll" -i 2 $INTERSECT_TXT_CONFIGS -b build/debug
    scripts/compare.py -c "$coll" -i 3 $INTERSECT_TXT_CONFIGS -b build/debug
}

test_int() {
    coll="$1"
    scripts/build_config.sh -d $INT_CONFIGS $INTERSECT_INT_CONFIGS
    scripts/compare.py -c "$coll" $INT_CONFIGS -b build/debug
    scripts/compare.py -c "$coll" -i 2 $INTERSECT_INT_CONFIGS -b build/debug
    scripts/compare.py -c "$coll" -i 3 $INTERSECT_INT_CONFIGS -b build/debug
}

scripts/build.sh -d gen_patterns

if [[ $# > 0 ]]; then
    coll="$1"
    if [[ -e "$coll/text_SURF.sdsl" ]]; then
        echo "Testing character-based collection"
        test_txt "$coll"
    else
        echo "Testing integer collection"
        test_int "$coll"
    fi
else
    test_txt collections/TEST_TXT
    test_int collections/TEST_INT
fi
