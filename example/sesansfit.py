import numpy as np
from bumps.names import *

from sasmodels import bumps_model as sas
kernel = sas.load_model('sphere', dtype='single')


if True: # fix when data loader exists
#    from sas.dataloader.readers\
    from sas.dataloader.loader import Loader
    loader=Loader()
    data=loader.load('testsasview1.ses')
    
#    data = load_sesans('mydatfile.pz')
#    sans_data = load_sans('mysansfile.xml')

else:
    SElength = np.linspace(0, 2400, 61) # [A]
    data = np.ones_like(SElength)
    err_data = np.ones_like(SElength)*0.03

    class SESANSData1D:
        #q_zmax = 0.23 # [A^-1]
        zacceptance = 0.1 # [A^-1]
        lam = 2e-10 # [m]
        thickness = 0.2 # [cm]
        x = SElength
        y = data
        dy = err_data
    data = SesansData()
print dir(data)

radius = 1000
data.Rmax = 3*radius # [A]

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

