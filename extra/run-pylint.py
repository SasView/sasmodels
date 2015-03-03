import os
import sys
from os.path import dirname, join as joinpath, abspath

def main():
    root = abspath(joinpath(dirname(__file__), '..'))
    os.chdir(root)
    cmd = "pylint --rcfile extra/pylint.rc -f parseable sasmodels > pylint_violations.txt"
    status = os.system(cmd)
    sys.exit(status)

if __name__ == "__main__":
    main()
