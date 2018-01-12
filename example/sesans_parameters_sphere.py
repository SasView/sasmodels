"""
This is a data file  used to load in sesans data and fit it using the bumps engine
"""
from bumps.names import *

import sesansfit

# Enter the model name to use
model_name = "sphere"

# DO NOT MODIFY THIS LINE
model = sesansfit.get_bumps_model(model_name)
model.radius.range(1,10000)

# Enter any custom parameters
# name = Parameter(initial_value, name='name')
phi = Parameter(0.10, name='phi')
# Add the parameters to this list that should be displayed in the fitting window
custom_params = {"phi" : phi}

# SESANS data file name
sesans_file = "sphere.ses"

# Initial parameter values (if other than defaults)
# "model_parameter_name" : value
initial_vals = {
    "sld" : 7.0,
    "radius" : 1000,
    "sld_solvent" : 1.0,
}

# Ranges for parameters if other than default
# "model_parameter_name" : [min, max]
param_range = {
    "phi" : [0.001, 0.5],
    "radius" : [1, 10000]
}

# Constraints
# model.param_name = f(other params)
# EXAMPLE: model.scale = model.radius*model.radius*(1 - phi) - where radius and scale are model functions and phi is
# a custom parameter
model.scale = phi*(1-phi)

# Send to the fitting engine
# DO NOT MODIFY THIS LINE
problem = sesansfit.sesans_fit(sesans_file, model, initial_vals, custom_params, param_range)
