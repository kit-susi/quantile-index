#!/bin/bash
set -xe
TXT_CONFIGS="BRUTE_TXT IDX_NN IDX_NN_DOCID_SMART IDX_NN_K3_DAAT"
INT_CONFIGS="BRUTE_INT IDX_NN_INT IDX_NN_DOCID_SMART_INT IDX_NN_K3_DAAT_INT"

scripts/build_config.sh -d $TXT_CONFIGS $INT_CONFIGS
scripts/build.sh -d gen_patterns

scripts/compare.py -c collections/TEST_TXT $TXT_CONFIGS -b build/debug
scripts/compare.py -c collections/TEST_INT $INT_CONFIGS -b build/debug
