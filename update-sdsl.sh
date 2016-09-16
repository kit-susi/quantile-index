#!/bin/bash
cd external/sdsl-lite && \
git checkout k2_algos && \
git pull && \
cd ../.. && \
git add external/sdsl-lite && \
git commit -m "forwarded sdsl-lite to current master"
