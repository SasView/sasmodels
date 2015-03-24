#!/usr/bin/env python

import os
import sys
from os.path import dirname, join as joinpath, abspath

def main():
    envpath = os.environ.get('PYTHONPATH',None)
    path = [envpath] if envpath else []
    path.append(abspath(dirname(__file__))) # so we can find the plugins
    os.environ['PYTHONPATH'] = ':'.join(path)
    root = abspath(joinpath(dirname(__file__), '..'))
    os.chdir(root)
    cmd = "pylint --rcfile extra/pylint.rc -f parseable sasmodels"
    status = os.system(cmd)
    sys.exit(status)

if __name__ == "__main__":
    main()
