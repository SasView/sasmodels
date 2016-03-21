"""
This is a data file  used to load in sesans data and fit it using the bumps engine
"""
from bumps.names import *

import sesansfit

# Enter the model name to use
model_name = "core_shell_sphere*hardsphere"

# DO NOT MODIFY THIS LINE
model = sesansfit.get_bumps_model(model_name)

# Enter any custom parameters
phi = Parameter(0.45, name='phi')
pen = Parameter(0.95, name='solvent penetration')
custom_params = {"phi" : phi, "pen" : pen}

# SESANS data file
sesans_file = "core_shell.ses"

# Initial parameter values (if other than defaults)
initial_vals = {
    "sld_core" : 1.05,
    "sld_shell" : 2.88*pen-0.05*(1-pen),
    "sld_solvent" : 2.88,
    "radius" : 730,
    "thickness" : 20,
    "volfraction" : phi,
    "scale" : (1-phi)
}

# Ranges for parameters if other than default
param_range = {
    "phi" : [0.2, 0.5],
    "pen" : [0,1],
    "radius" : [500, 3000],
    "thickness" : [0,200]
}

# Constraints
# model.param_name = f(other params)
# EXAMPLE: model.scale = model.radius*model.radius*(1 - phi) - where radius and scale are model functions and phi is
# a custom parameter
model.scale = phi*(1-phi)
model.volfraction = phi
model.shell_sld = pen*2.88

# Send to the fitting engine
problem = sesansfit.sesans_fit(sesans_file, model_name, initial_vals, custom_params, param_range)
