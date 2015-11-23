#!/bin/sh

SASVIEW=$(ls -d ../sasview/build/lib.*)
PYTHONPATH=../bumps:../periodictable:$SASVIEW
export PYTHONPATH

#echo PYTHONPATH=$PYTHONPATH

python -m sasmodels.compare_many $*

#for dim in 1d100 2d32; do
#    python -m sasmodels.compare_many all 100 $dim 0 > ${dim}_poly_0.csv
#    python -m sasmodels.compare_many all 100 $dim 1e-5 > ${dim}_poly_1e-5.csv
#    python -m sasmodels.compare_many all 100 $dim 1e-3 > ${dim}_poly_1e-3.csv
#    python -m sasmodels.compare_many all 100 $dim mono > ${dim}_mono.csv
#done

