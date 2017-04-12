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

## Index construction

To construct a specific index with the configuration 'IDX' and a collection
directory 'COLDIR'.

    $ ./build/release/surf_build_IDX -c COLDIR

Make sure that the collection directory contains a file
named 'text_SURF.sdsl' with a sdsl::int_vector. The file should the
concatenation of all documents separated by '\1'.

## Reproducing experiments
To the experiments of the work
'The Quantile Index - Succinct Self-Index for Top-k Document retrieval' by
Niklas Baumstark, Simon Gog, Tobias Heuer and Julian Labeit use the following
instructions.

First build all binaries as descibed above.
Then download the collections by running

    $ ./scripts/download-collections.sh

After completion the directory ./collections should contain the collections
'ENWIKISML', ENWIKIBIG', 'REVISIONS' and 'SOURCES'.

Finnaly run the experiments by executing the scripts in the directory
'./experiments'.  The directory with all the indexes have to be specified
through the command line parameters of the 'run' scripts. Each script should
output a *.csv file with the results.

In case of any questions please contact us via the github issue or e-mail.
