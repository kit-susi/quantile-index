#!/bin/bash
#
# Usage:
# scripts/build_config.sh [-d] [config [config ...]]
#
# Flags:
# -d: Make debug build with -O1 -g2 and assertions
#
# Example:
# scripts/build_config.sh IDX_NN BRUTE_TXT
set -e

build_args=
if [[ "$1" == "-d" ]]; then
  build_args="$build_args -d"
  shift 1
fi

for config in $@; do
  build_args="$build_args surf_index-$config surf_query-$config"
done

scripts/build.sh $build_args
