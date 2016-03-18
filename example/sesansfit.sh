#!/bin/bash

# Need to fix the paths to sasmodels and sasview if no eggs present
echo $PWD
ONEUP="$(dirname "$PWD")"
PROJECTS="$(dirname "$ONEUP")"
CCOLON="C:/"
CSLASH="/c/"
SASMODELSBASE=$PROJECTS/sasmodels/
SASMODELS="${SASMODELSBASE/$CSLASH/$CCOLON}"
SASVIEWBASE=$PROJECTS/sasview/src/
SASVIEW="${SASVIEWBASE/$CSLASH/$CCOLON}"
PYTHONPATH="$SASVIEW;$SASMODELS"
export PYOPENCL_CTX PYTHONPATH

echo PYTHONPATH=$PYTHONPATH
set -x

python -m bumps.cli $*
