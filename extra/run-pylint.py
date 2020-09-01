#!/usr/bin/env python

import os
import sys
from os.path import dirname, join as joinpath, abspath

def main():
    extra = abspath(dirname(__file__))
    root = abspath(joinpath(extra, '..'))

    envpath = os.environ.get('PYTHONPATH', None)
    path = [envpath] if envpath else []
    path.append(extra)

    #bumps = abspath(joinpath(root, "..", "bumps"))
    #periodictable = abspath(joinpath(root, "..", "periodictable"))
    #sasview = abspath(joinpath(root, "..", "sasview", "src"))
    #path.extend((bumps, periodictable, sasview))

    os.environ['PYTHONPATH'] = ':'.join(path)

    # Run the lint command
    cmd = "pylint --rcfile extra/pylint.rc -f parseable sasmodels"
    os.chdir(root)
    status = os.system(cmd)
    sys.exit(status)

if __name__ == "__main__":
    main()
