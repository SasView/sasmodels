"""
This is a data file used to load in sesans data and fit it using the bumps engine
"""
from bumps.names import *

import sesansfit

# Enter the model name to use
model_name = "sphere"

# DO NOT MODIFY THIS LINE
model = sesansfit.get_bumps_model(model_name)

# Enter any custom parameters
# name = Parameter(initial_value, name='name')
phi = Parameter(0.0855, name='phi')
# Add the parameters to this list that should be displayed in the fitting window
custom_params = {"phi" : phi}

# SESANS data file name
sesans_file = "spheres2micron.ses"

# Initial parameter values (if other than defaults)
# "model_parameter_name" : value
initial_vals = {
    "sld" : 1.41,
    "radius" : 10000,
    "sld_solvent" : 2.70,
}

# Ranges for parameters if other than default
# "model_parameter_name" : [min, max]
param_range = {
    "phi" : [0.001, 0.5],
    "radius" : [100, 100000]
}

# Constraints
# model.param_name = f(other params)
# EXAMPLE: model.scale = model.radius*model.radius*(1 - phi) - where radius
# and scale are model functions and phi is a custom parameter
model.scale = phi*(1-phi)

# Send to the fitting engine
# DO NOT MODIFY THIS LINE
problem = sesansfit.sesans_fit(sesans_file, model, initial_vals, custom_params, param_range)
