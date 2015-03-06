#!/bin/sh

SASVIEW=$(ls -d ../sasview/build/lib.*)
PYTHONPATH=../bumps:../periodictable:$SASVIEW
export PYTHONPATH

echo PYTHONPATH=$PYTHONPATH
set -x

for dim in 1d100 2d32; do
    ./compare_many.py all 100 $dim 0 > ${dim}_poly_0.csv
    ./compare_many.py all 100 $dim 1e-5 > ${dim}_poly_1e-5.csv
    ./compare_many.py all 100 $dim 1e-3 > ${dim}_poly_1e-3.csv
    ./compare_many.py all 10000 $dim mono > ${dim}_mono.csv
done
