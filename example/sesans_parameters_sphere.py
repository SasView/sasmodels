"""
This is a data file  used to load in sesans data and fit it using the bumps engine
"""
from bumps.names import *

import sesansfit

# Enter the model name to use
model_name = "sphere"

# Enter any custom parameters
phi = Parameter(0.10, name='phi')
custom_params = {"phi" : phi}

# SESANS data file
sesans_file = "sphere.ses"

# Initial parameter values (if other than defaults)
initial_vals = {
    "scale" : phi*(1 - phi),
    "sld" : 7.0,
    "solvent_sld" : 1.0,
    "radius" : 1000,
}

# Ranges for parameters if other than default
param_range = {
    "phi" : [0.001, 0.5],
    "radius" : [1, 10000]
}

# Send to the fitting engine
problem = sesansfit.sesans_fit(sesans_file, model_name, initial_vals, custom_params, param_range)