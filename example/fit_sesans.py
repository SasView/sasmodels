#!/usr/bin/env python
import os
import sys

# Need to fix the paths to sasmodels and sasview if no eggs present
ONEUP=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECTS=os.path.dirname(ONEUP)
SASMODELS=os.path.join(PROJECTS, 'sasmodels')
SASVIEW=os.path.join(PROJECTS, 'sasview', 'src')
BUMPS=os.path.join(PROJECTS, 'bumps')

sys.path.insert(0, BUMPS)
sys.path.insert(0, SASVIEW)
sys.path.insert(0, SASMODELS)

if sys.argv[-1].startswith('-'):
    from bumps.cli import main as run_bumps
else:
    from bumps.gui.gui_app import main as run_bumps
run_bumps()
