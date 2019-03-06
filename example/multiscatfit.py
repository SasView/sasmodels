#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fit model using multiple scattering.

As of this writing, multiscattering isn't integrated into sasmodels, and a
number of hacks are required to get it running as a fit.

The appropriate items need to be on the python path.  These include
sasview (for reading the data), sasmodels and bumps.  The multiscat module
(currently in the sasmodels/explore directory) is also needed, either beside
this example fit file, or by putting sasmdoels/explore on the python path.

On Unix/Mac running as developer I do::

    # Show the model without fitting
    PYTHONPATH=..:../../bumps:../../sasview/src python multiscatfit.py

    # Run the fit
    PYTHONPATH=..:../../bumps:../../sasview/src ../../bumps/run.py \
    multiscatfit.py --store=/tmp/t1

You may be able to run multiscatfit.py against the distributed sasview
application (if it is new enough, and if you have multiscat.py in the
same directory).  You probably need a command such as::

    sasview.exe bumps.cli multiscatfit.py --store=t1
"""

import sys
from bumps.names import *
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasmodels.data import load_data, set_beam_stop, set_top

from sasmodels.multiscat import MultipleScattering

## Load the data
#data = load_data('DEC07267.DAT')
#set_beam_stop(data, 0.003, outer=0.025)
data = load_data('latex_smeared.xml', index=0)

## Define the model
kernel = load_model("ellipsoid")

model = Model(
    kernel,
    scale=0.005, background=0.05,
    radius_polar=2200, radius_equatorial=2200,
    sld=.291, sld_solvent=7.105,
    #theta=90, theta_pd=0, theta_pd_n=0, theta_pd_nsigma=3,
    #phi=90, phi_pd=0, phi_pd_n=20, phi_pd_nsigma=3,
    radius_polar_pd=0.222296, radius_polar_pd_n=1, radius_polar_pd_nsigma=0,
    radius_equatorial_pd=.000128, radius_equatorial_pd_n=1, radius_equatorial_pd_nsigma=0,
    )

# Tie the model to the data
M = Experiment(data=data, model=model)

# Stack mulitple scattering on top of the existing resolution function.
M.resolution = MultipleScattering(resolution=M.resolution, probability=0.)

# SET THE FITTING PARAMETERS
model.radius_polar.range(15, 3000)
model.radius_equatorial.range(15, 3000)
#model.theta.range(0, 90)
#model.theta_pd.range(0,10)
#model.phi_pd.range(0,20)
#model.phi.range(0, 180)
model.background.range(0,1000)
model.scale.range(0, 0.1)

# The multiple scattering probability parameter is in the resolution function
# instead of the scattering function, so access it through M.resolution
M.scattering_probability.range(0.0, 0.9)

# Let bumps know that we are fitting this experiment
problem = FitProblem(M)

if __name__ == "__main__":
    #M.theory()
    M.plot()
    import pylab; pylab.show()
