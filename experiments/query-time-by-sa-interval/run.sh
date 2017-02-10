#!/bin/bash
set -e

configs="IDX_NN_16 IDX_NN_LG_16 IDX_NN_QUANTILE_16_64 IDX_NN_QUANTILE_LG_16_64"
colls="ENWIKIBIG REVISIONS SOURCES ENWIKISML"
config_intervalsz=IDX_NN_QUANTILE_LG_16_64
patlen="3-10"
k=10
queries=500000
seed=1
numactl=~/numactl

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

echo "instance;algo;intervalsz;time;result_count" > "$result_file"

for coll in $colls; do
  query_file="$tmpdir/queries_$coll"
  echo "Building queries for $coll into $query_file"
  coll_dir=collections/$coll
  scripts/gen_queries.py --seed $seed -c "$coll_dir" -q $queries -n "$patlen" "$query_file"
  for config in $configs; do
    echo "Running $config"
    build/release/surf_index-$config -c "$coll_dir" > /dev/null
    extra_args=()
    if [[ $config == $config_intervalsz ]]; then
      extra_args+=(-d "$tmpdir/intervalsz_raw_$coll")
    fi
    $numactl --physcpubind=0 --membind=0 \
      build/release/surf_query-$config -k $k -c "$coll_dir" -q "$query_file" -t "${extra_args[@]}" \
        | grep TIME | cut -d';' -f2,3 > "$tmpdir/timing_${coll}_${config}"
    yes "$coll;$config" | head -n $queries > "$tmpdir/common_cols_${coll}_${config}"
  done
  cut -d';' -f2 "$tmpdir/intervalsz_raw_$coll" > "$tmpdir/intervalsz_$coll"
  for config in $configs; do
    paste -d';' "$tmpdir/common_cols_${coll}_$config" \
                "$tmpdir/intervalsz_$coll" \
                "$tmpdir/timing_${coll}_$config" \
      >> "$result_file"
  done
done

cd "$dir"
pwd
Rscript plot.R
