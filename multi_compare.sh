#!/bin/sh

sasview=( ../sasview/build/lib.* )
sep=$(python -c "import os;print(os.pathsep)")
PYTHONPATH=../bumps${sep}../periodictable${sep}$sasview
export PYTHONPATH
#echo PYTHONPATH=$PYTHONPATH

#PYOPENCL_COMPILER_OUTPUT=1; export PYOPENCL_COMPILER_OUTPUT
#PYOPENCL_CTX=${CTX:-1}; export PYOPENCL_CTX

#SAS_OPENMP=1; export SAS_OPENMP

python -m sasmodels.compare_many $*

#for dim in 1d100 2d32; do
#    python -m sasmodels.compare_many all 100 $dim 0 > ${dim}_poly_0.csv
#    python -m sasmodels.compare_many all 100 $dim 1e-5 > ${dim}_poly_1e-5.csv
#    python -m sasmodels.compare_many all 100 $dim 1e-3 > ${dim}_poly_1e-3.csv
#    python -m sasmodels.compare_many all 100 $dim mono > ${dim}_mono.csv
#done

