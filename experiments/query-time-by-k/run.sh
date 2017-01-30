#!/bin/bash
set -e

configs="IDX_NN_16 IDX_NN_LG_16 IDX_NN_QUANTILE_LG_16_32 IDX_NN_QUANTILE_LG_16_64"
colls="TEST_TXT"
patlen=5
ks="1 2 4 8 16 32 64 128 256 512 1024 2048"
queries=10000
seed=1

dir=`dirname "$0"`
cd "$dir"
dir="`pwd`"
root="$dir/../.."
cd "$root"

result_file="$dir/results.csv"

scripts/build_config.sh $configs > /dev/null

tmpdir=`mktemp -d`
echo "Tmpdir is $tmpdir"
function cleanup() {
  [ -z "$tmpdir" ] && return
  echo "Cleaning up $tmpdir"
  rm -rf "$tmpdir"
}
trap cleanup EXIT

echo "instance;algo;k;time;result_count" > "$result_file"

for coll in $colls; do
  query_file="$tmpdir/queries_$coll"
  echo "Building queries for $coll into $query_file"
  coll_dir=collections/$coll
  scripts/gen_queries.py --seed $seed -c "$coll_dir" -q $queries -n "$patlen" "$query_file"
  echo Done
  for config in $configs; do
    for k in $ks; do
      echo "Running $config with k=$k"
      build/release/surf_index-$config -c "$coll_dir" > /dev/null
      build/release/surf_query-$config -k $k -c "$coll_dir" -q "$query_file" -t \
          | grep TIME | cut -d';' -f2,3 > "$tmpdir/timing_${coll}_${config}_${k}"
      yes "$coll;$config;$k" | head -n $queries > "$tmpdir/common_cols"
      paste -d';' "$tmpdir/common_cols" "$tmpdir/timing_${coll}_${config}_$k" \
        >> "$result_file"
    done
  done
done

cd "$dir"
pwd
Rscript plot.R
