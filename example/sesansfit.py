import numpy as np
from bumps.names import *

from sasmodels import bumps_model as sas
kernel = sas.load_model('sphere', dtype='single')
#kernel = sas.load_model('triaxial_ellipsoid', dtype='single')


if True: # fix when data loader exists
#    from sas.dataloader.readers\
    from sas.dataloader.loader import Loader
    loader=Loader()
    data=loader.load('testsasview1.ses')
    data.x /=10
    
#    data = load_sesans('mydatfile.pz')
#    sans_data = load_sans('mysansfile.xml')

else:
    SElength = np.linspace(0, 2400, 61) # [A]
    data = np.ones_like(SElength)
    err_data = np.ones_like(SElength)*0.03

    class Sample:
        zacceptance = 0.1 # [A^-1]
        thickness = 0.2 # [cm]
        
    class SESANSData1D:
        #q_zmax = 0.23 # [A^-1]
        lam = 0.2 # [nm]
        x = SElength
        y = data
        dy = err_data
        sample = Sample()
    data = SESANSData1D()

radius = 1000
data.Rmax = 3*radius # [A]

##  Sphere parameters

phi = Parameter(0.1, name="phi")
model = sas.BumpsModel(data, kernel,
    scale=phi*(1-phi), sld=7.0, solvent_sld=1.0, radius=radius)
phi.range(0.001,0.90)
#model.radius.pmp(40)
model.radius.range(100,10000)
#model.sld.pmp(5)
#model.background
#model.radius_pd=0
#model.radius_pd_n=0

if False: # have sans data
    sansmodel = sas.BumpsModel(sans_data, kernel, **model.parameters())
    problem = FitProblem([model, sansmodel])
else:
    problem = FitProblem(model)


### Tri-Axial Ellipsoid
#
#phi = Parameter(0.1, name='phi')
#model = sas.BumpsModel(data, kernel,
#    scale=phi*(1-phi), sld=7.0, solvent_sld=1.0, radius=radius)
#phi.range(0.001,0.90)
##model.radius.pmp(40)
#model.radius.range(100,10000)
##model.sld.pmp(5)
##model.background
##model.radius_pd=0
##model.radius_pd_n=0
#
#if False: # have sans data
#    sansmodel = sas.BumpsModel(sans_data, kernel, **model.parameters())
#    problem = FitProblem([model, sansmodel])
#else:
#    problem = FitProblem(model)
