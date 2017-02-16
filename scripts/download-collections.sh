#!/bin/bash
set -e
url=http://quantile.s3.amazonaws.com/collections
colls=(ENWIKISML SOURCES REVISIONS ENWIKIBIG)
sha1sums=(
  fbfd9dca28595e86435622d7a0987bc6478781a5
  893eb57ebe951fbcf3079e4b09f6a150c7421805
  fc901b6116c8c540fa83b6e369120398013acb31
  7043dfeb9440921c38b289e35d916ea0668eea20
)
cd `dirname "$0"`/..

for (( i=0; i<${#colls[@]}; i++ )); do
    coll=${colls[$i]}
    file=collections/$coll/text_SURF.sdsl
    sha1=${sha1sums[$i]}
    if [ -e $file ]; then
        echo "$coll already existing, checking hash"
        has_sha1=(`sha1sum "$file"`)
        if [[ "$has_sha1" != "$sha1" ]]; then
            echo 2>&1 "File $file has wrong SHA-1 checksum (should be $sha1)"
        fi
    else
        mkdir -p "collections/$coll"
        wget "$URL/$coll/text_SURF.sdsl" -O "$file"
    fi
done
