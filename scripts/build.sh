#!/bin/bash
#
# Usage:
# scripts/build.sh [-d] [target [target ...]]
#
# Flags:
# -d: Make debug build with -O1 -g2 and assertions
#
# Example:
# scripts/build.sh surf_query-IDX_NN gen_patterns
set -e

build_type=Release
build_dir=build/release
if [[ "$1" == "-d" ]]; then
  build_type=Debug
  build_dir=build/debug
  shift 1
fi

mkdir -p $build_dir
cd $build_dir
set -x
cmake -DCMAKE_BUILD_TYPE=$build_type ../..
set +x

echo "Building $@"
for target in $@; do
  set -x
  make -j16 $target &
  set +x
done

for job in `jobs -p`; do
  wait $job
done
