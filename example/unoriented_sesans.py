from bumps.names import *

from sasmodels.data import load_data
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment

# Spherical particle data, not ellipsoids
data = load_data('../../sasview/sasview/test/sesans_data/sphere2micron.ses')
kernel = load_model("sphere*hardsphere")

model = Model(
    kernel,
    scale=0.08, background=0,
    sld=.291, sld_solvent=7.105,
    radius=10000, radius_pd=0.222296, radius_pd_n=0,
    volfraction=0.2,
    )

# SET THE FITTING PARAMETERS
model.radius.range(1000, 100000)
model.volfraction.range(0, 1)
model.background.range(0, 1000)
model.scale.range(0, 10)

M = Experiment(data=data, model=model)
problem = FitProblem(M)
