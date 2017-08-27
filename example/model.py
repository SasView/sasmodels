#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from bumps.names import *
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasmodels.data import load_data, set_beam_stop, set_top

""" IMPORT THE DATA USED """
radial_data = load_data('DEC07267.DAT')
set_beam_stop(radial_data, 0.00669, outer=0.025)
set_top(radial_data, -.0185)

kernel = load_model("ellipsoid")

model = Model(kernel,
    scale=0.08,
    radius_polar=15, radius_equatorial=800,
    sld=.291, sld_solvent=7.105,
    background=0,
    theta=90, phi=0,
    theta_pd=15, theta_pd_n=40, theta_pd_nsigma=3,
    radius_polar_pd=0.222296, radius_polar_pd_n=1, radius_polar_pd_nsigma=0,
    radius_equatorial_pd=.000128, radius_equatorial_pd_n=1, radius_equatorial_pd_nsigma=0,
    phi_pd=0, phi_pd_n=20, phi_pd_nsigma=3,
    )

# SET THE FITTING PARAMETERS
model.radius_polar.range(15, 1000)
model.radius_equatorial.range(15, 1000)
model.theta_pd.range(0, 360)
model.background.range(0,1000)
model.scale.range(0, 10)

#cutoff = 0     # no cutoff on polydisperisity loops
#cutoff = 1e-5  # default cutoff
cutoff = 1e-3  # low precision cutoff
M = Experiment(data=radial_data, model=model, cutoff=cutoff)
problem = FitProblem(M)