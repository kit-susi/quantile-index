#!/bin/bash
#
# Usage:
# scripts/build.sh [-d] [config [config ...]]
#
# Flags:
# -d: Make debug build with -O1 -g2 and assertions
#
# Example:
# scripts/build.sh IDX_NN BRUTE_TXT
set -e

build_type=Release
if [[ "$1" == "-d" ]]; then
  build_type=Debug
  shift 1
fi

mkdir -p build
cd build
set -x
cmake -DCMAKE_BUILD_TYPE=$build_type ..
set +x

targets=
for config in $@; do
  targets="$targets surf_index-$config surf_query-$config"
done

set -x
eval make -j $targets
