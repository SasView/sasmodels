#!/bin/bash

SASVIEW=$PWD/../../sasview/src
PYTHONPATH=$PWD/../:$PWD/../../bumps:$PWD/../../periodictable:$SASVIEW
export PYOPENCL_CTX PYTHONPATH

echo PYTHONPATH=$PYTHONPATH
set -x

python -m bumps.cli $*
