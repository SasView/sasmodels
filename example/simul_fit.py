#!/usr/bin/env python
# -*- coding: utf-8 -*-

# To Sasview/documents/scripts

from bumps.names import *
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasmodels.data import load_data, plot_data


""" IMPORT THE DATA USED """
datafiles = ['latex_smeared_out_0.txt', 'latex_smeared_out_1.txt']
datasets = [load_data(el) for el in datafiles]

for data in datasets:
    data.qmin = 0.0
    data.qmax = 10.0

#sphere model
kernel = load_model('sphere', dtype="single")
pars = dict(scale=0.01, background=0.0, sld=1.0, sld_solvent=6.0, radius=1500.)
model = Model(kernel, **pars)
model.radius.range(0, inf)
#model.background.range(-inf, inf)
#model.scale.range(0, inf)
model.sld.range(-inf, inf)
model.sld_solvent.range(-inf, inf)

free = FreeVariables(
    names=[data.filename for data in datasets],
    background=model.background,
    scale=model.scale,
    )
free.background.range(-inf, inf)
free.scale.range(0, inf)

M = [Experiment(data=data, model=model) for data in datasets]

problem = FitProblem(M, freevars=free)

print(problem._parameters)
