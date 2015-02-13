#!/bin/sh

#SASVIEW=$(ls -d ../sasview/build/lib.*)
SASVIEW=( ../sasview/build/lib.* )
#PYOPENCL_CTX=${CTX:-1}
PYTHONPATH=../bumps:../periodictable:$SASVIEW
export PYOPENCL_CTX PYTHONPATH

echo PYTHONPATH=$PYTHONPATH
set -x

./compare.py $*
