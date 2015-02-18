import numpy as np
from bumps.names import *

from sasmodels import bumps_model as sas
kernel = sas.load_model('sphere', dtype='single')


radius = 1000
if False: # fix when data loader exists
    from sasmodels.dataloader import load_sesans, load_sans
    data = load_sesans('mydatfile.pz')
    sans_data = load_sans('mysansfile.xml')
else:
    SElength = np.linspace(0, 2400, 61) # [A]
    data = np.ones_like(SElength)
    err_data = np.ones_like(SElength)*0.03

    class SesansData:
        #q_zmax = 0.23 # [A^-1]
        q_zmax = 0.1 # [A^-1]
        SElength = np.linspace(0, 2400, 61) # [A]
        Rmax = 3*radius # [A]
        wavelength = 2e-10 # [m]
        thickness = 0.2 # [cm]
        SElength = SElength
        data = data
        err_data = err_data
    data = SesansData()

phi = Parameter(0.1, name="phi")
model = sas.BumpsModel(data, kernel,
    scale=phi*(1-phi), sld=7.0, solvent_sld=1.0, radius=radius)
phi.pmp(10)
model.radius.pmp(40)
model.sld.pm(2)
model.background.range(0,5)


if False: # have sans data
    sansmodel = sas.BumpsModel(sans_data, kernel, **model.parameters())
    problem = FitProblem([model, sansmodel])
else:
    problem = FitProblem(model)

