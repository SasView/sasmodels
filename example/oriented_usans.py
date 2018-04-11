from bumps.names import *

from sasmodels.data import load_data
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment

# Spherical particle data, not ellipsoids
sans, usans = load_data('latex_smeared.xml', index='all')
usans.qmin, usans.qmax = np.min(usans.x), np.max(usans.x)
usans.mask = (usans.x < 0.0)
usans.oriented = True
#print sans.dxl, usans.dxl
#import pprint; pprint.pprint(sans.__dict__)

kernel = load_model("ellipsoid")

model = Model(
    kernel,
    scale=0.08, background=0,
    sld=.291, sld_solvent=7.105,
    radius_polar=1800, radius_polar_pd=0.222296, radius_polar_pd_n=0,
    radius_equatorial=2600, radius_equatorial_pd=0.28, radius_equatorial_pd_n=0,
    theta=60, theta_pd=0, theta_pd_n=0,
    phi=60, phi_pd=0, phi_pd_n=0,
    )

# SET THE FITTING PARAMETERS
model.radius_polar.range(1000, 10000)
model.radius_equatorial.range(1000, 10000)
model.theta.range(0, 360)
model.phi.range(0, 360)
model.background.range(0,1000)
model.scale.range(0, 10)


#cutoff = 0     # no cutoff on polydisperisity loops
#cutoff = 1e-5  # default cutoff
cutoff = 1e-3  # low precision cutoff
M = Experiment(data=usans, model=model, cutoff=cutoff)
problem = FitProblem(M)

