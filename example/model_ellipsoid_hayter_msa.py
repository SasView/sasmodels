import sys
import numpy as np
#sys.path.append('path_to_sasmodels')


from bumps.fitProblem import FitProblem
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasmodels.data import load_data

# IMPORT THE DATA USED
data = load_data(sys.argv[1])

#setattr(data, 'qmin', 0.0)
#setattr(data, 'qmax', 10.0)

# DEFINE THE MODEL
kernel = load_model('ellipsoid@hayter_msa')

pars = dict(scale=6.4, background=0.06, sld=0.33, sld_solvent=2.15, radius_polar=14.0,
            radius_equatorial=24.0, volfraction=0.075, charge=66.373, temperature=298.0,
            concentration_salt=0.001, dielectconst=71.0)

model = Model(kernel, **pars)

# PARAMETER RANGES (ONLY THOSE PARAMETERS ARE FITTED)
model.scale.range(0, np.inf)
model.background.range(-np.inf, np.inf)
#model.sld.range(-np.inf, np.inf)
model.sld_solvent.range(-np.inf, np.inf)
#model.radius_polar.range(0, np.inf)
#model.radius_equatorial.range(0, np.inf)
#model.volfraction.range(0,0.74)
#model.charge.range(0, np.inf)
#model.temperature.range(0,1000)
#model.concentration_salt.range(0, 1)
#model.dielectconst.range(0,np.inf)

M = Experiment(data=data, model=model)

problem = FitProblem(M)
