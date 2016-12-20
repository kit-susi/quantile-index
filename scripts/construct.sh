#!/bin/bash
set -xe
TXT_CONFIGS="IDX_NN IDX_NN_QUANTILE IDX_NN_QUANTILE_LG IDX_NN_QUANTILE_LG_16"

scripts/build_config.sh $TXT_CONFIGS

coll=collections/TEST_TXT
if [[ $# > 0 ]]; then
    coll="$1"
fi

for CONFIG in $TXT_CONFIGS
do
    build/release/surf_index-$CONFIG -c "$coll"
done

set +x
# Output sizes
for CONFIG in $TXT_CONFIGS
do 
    echo $CONFIG
    build/release/surf_index-$CONFIG -c "$coll" -m m | tail -n1
done
