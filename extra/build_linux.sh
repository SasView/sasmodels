#!/bin/bash

EASY_INSTALL=${EASY_INSTALL:-`which easy_install`}
cd $WORKSPACE
if [ ! -d "utils" ]; then
    mkdir utils
fi
export PYTHONPATH=$WORKSPACE/utils:$PYTHONPATH

"$EASY_INSTALL" -d "$WORKSPACE/utils" unittest-xml-reporting pylint

cd $WORKSPACE
python -m sasmodels.model_test opencl_and_dll all || exit 1

python extra/run-pylint.py > pylint_violations.txt || exit 0


