#!/bin/sh

SASVIEW=$(ls -d ../sasview/build/lib.*)
PYTHONPATH=../bumps:../periodictable:$SASVIEW
export PYTHONPATH

#echo PYTHONPATH=$PYTHONPATH

./compare_many.py $*

#for dim in 1d100 2d32; do
#    ./compare_many.py all 100 $dim 0 > ${dim}_poly_0.csv
#    ./compare_many.py all 100 $dim 1e-5 > ${dim}_poly_1e-5.csv
#    ./compare_many.py all 100 $dim 1e-3 > ${dim}_poly_1e-3.csv
#    ./compare_many.py all 100 $dim mono > ${dim}_mono.csv
#done

