# SuSI - SUccint Search Index

## Directory structure

* `build/release` + `build/debug`: Release/Debug build
* `collections`: contains collection data; each collection in its own subdirectory.
* `external`: Contains external libraries
* `external/sdsl-lite`: Mirror of git@github.com:kit-susi/sdsl-lite master
* `include/surf`: Contains headers.
* `src`: Contains surf sources.
  - `surf_index.cpp` - Build an index
  - `surf_query.cpp` - Query an index
* `scripts`:
  - `build.sh`/`build_config.sh`: Build a binary / index config
  - `smoke_test.sh`: Test all important index implementations for correctness
  - `benchmark.sh`: Quick performance overview
  - `format.sh`: Quickly format a C++ file (header or source)
  - `sdsl.sh`: pull or push from/to susi-sdsl-lite

## Installation

    $ git clone git@github.com:kit-susi/susi.git susi
    $ cd susi
    $ git submodule init
    $ git submodule update --recursive

## Building

To do a release build:

    $ scripts/build.sh

Or just a selected subset of configs:

    $ scripts/build.sh IDX_NN IDX_NN_K3_DAAT

For debug builds, add the `-d` flag to `build.sh` or `build_config.sh`. The two
builds will live in `build/release` and `build/debug`, respectively.

## Keeping `external/sdsl-lite/` up to date

Use the shell script `scripts/sdsl.sh`. We are mirroring
git@github.com:kit-susi/sdsl-lite.git (branch master) in
subtree `external/sdsl-lite/`. See
[https://git-scm.com/book/en/v1/Git-Tools-Subtree-Merging](subtree merging)
for an overview over the subtree merging workflow.

