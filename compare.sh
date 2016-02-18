#!/bin/bash

sasview=( ../sasview/build/lib.* )
sep=$(python -c "import os;print(os.pathsep)")
PYTHONPATH=../bumps${sep}../periodictable${sep}$sasview
export PYTHONPATH
#echo PYTHONPATH=$PYTHONPATH

#PYOPENCL_COMPILER_OUTPUT=1; export PYOPENCL_COMPILER_OUTPUT
#PYOPENCL_CTX=${CTX:-1}; export PYOPENCL_CTX

#SAS_OPENMP=1; export SAS_OPENMP

${PYTHON:-python} -m sasmodels.compare $*
