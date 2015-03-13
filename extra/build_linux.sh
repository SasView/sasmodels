#!/bin/bash

test -z $WORKSPACE && WORKSPACE=.

EASY_INSTALL=${EASY_INSTALL:-`which easy_install`}
cd $WORKSPACE
if [ ! -d "utils" ]; then
    mkdir utils
fi
export PYTHONPATH=$WORKSPACE/utils:$PYTHONPATH

"$EASY_INSTALL" -d "$WORKSPACE/utils" unittest-xml-reporting pylint

# Run tests
STATUS=0
python -m sasmodels.model_test opencl_and_dll all || STATUS=$?
python -m sasmodels.resolution || STATUS=$?

# Run lint
python extra/run-pylint.py > pylint_violations.txt

exit $STATUS


