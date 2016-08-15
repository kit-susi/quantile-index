#!/bin/bash

for x in `seq 1 3`; do
    echo "${x}"
    rm results/experiment2.*.byte.txt
    rm results/experiment2.txt
    make experiment2
    cat results/experiment2.txt >> results/all2.txt
done
