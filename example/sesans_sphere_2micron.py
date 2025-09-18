"""
This is a data file used to load in sesans data and fit it using the bumps engine

Usage:

    bumps sesans_sphere_2micron.py
"""
from bumps.names import FitProblem, Parameter

from sasdata import data_path

from sasmodels.bumps_model import Experiment, Model
from sasmodels.core import load_model
from sasmodels.data import load_data

# Enter the model name and the datafile path
model_name = "sphere"
data_file = data_path / "sesans_data" / "sphere2micron.ses"

# Custom parameters for use in expressions
# name = Parameter(initial_value, name='name')
phi = Parameter(0.0855, name='phi')  # scale = phi*(1-phi)

# Initial parameter values and expressions (if other than defaults)
# "model_parameter_name" : value or expression
pars = {
    "scale": phi*(1-phi),
    "sld" : 1.41,
    "radius" : 10000,
    "sld_solvent" : 2.70,
}

# DO NOT MODIFY THIS LINE
model = Model(load_model(model_name), **pars)

# Bounds constraints
# model.param_name = f(other params)
model.radius.range(100, 100000)
phi.range(0.001, 0.5)


# Send to the fitting engine
# DO NOT MODIFY THESE LINES
data = load_data(str(data_file))
M = Experiment(data=data, model=model)
problem = FitProblem([M])

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    print(f"==== {model_name} model for {data_file.name} has χ² = {problem.chisq_str()} ====")
    print(problem.summarize())
    problem.plot()
    plt.show()
