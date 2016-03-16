"""
This is a data file  used to load in sesans data and fit it using the bumps engine
"""
from bumps.names import *

import sesansfit

# Enter the model name to use
model_name = "core_shell_sphere*hardsphere"

# Enter any custom parameters
phi = Parameter(0.45, name='phi')
pen = Parameter(0.95, name='solvent penetration')
custom_params = {"phi" : phi, "pen" : pen}

# SESANS data file
sesans_file = "core_shell.ses"

# Initial parameter values (if other than defaults)
initial_vals = {
    "scale" : 0.09,
    "core_sld" : 1.0592,
    "solvent_sld" : 2.88,
    "shell_sld" : 2.88,
    "radius" : 890,
    "thickness" : 130,
    "volfraction" : 0.45
}

# Ranges for parameters if other than default
param_range = {
    "phi" : [0.2, 0.5],
    "pen" : [0,1],
    "radius" : [500, 3000],
    "thickness" : [0,200]
}

# Send to the fitting engine
problem = sesansfit.sesans_fit(sesans_file, model_name, initial_vals, custom_params, param_range)