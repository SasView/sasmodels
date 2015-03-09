#!/bin/bash

#SASVIEW=$(ls -d ../sasview/build/lib.*)
SASVIEW=( ../sasview/build/lib.* )
#PYOPENCL_CTX=${CTX:-1}
PYTHONPATH=../bumps:../periodictable:$SASVIEW
export PYTHONPATH

set -x

./compare.py $*
