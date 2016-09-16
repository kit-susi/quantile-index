# SuSIE - SUccint Search & Information retrieval Engine?

## Directory structure

    * build:
    * collections: contains collection data;
        each collection in its own subdirectory.
    * external: Contains external libraries
    * include/surf: Contains headers.
    * src: Contains surf sources.
      - surf_index.cpp
      - surf_search.cpp
      - surf_profile.cpp
    * tools: Contains source of tools. E.g.
      - qry2intqry
      - indri2surf converter

## Installation

    * git clone git@github.com:niklasb/susi.git susie
    * cd susie
    * git submodule init
    * git submodule update --recursive
    * cd build
    * cmake ..

## Keeping `external/sdsl-lite/` up to date

Use the shell script `scripts/sdsl.sh`. We are mirroring
git@github.com:niklasb/susi-sdsl-lite.git (branch master) in
subtree `external/sdsl-lite/`. See
[https://git-scm.com/book/en/v1/Git-Tools-Subtree-Merging](subtree merging)
for an overview over the subtree merging workflow.
