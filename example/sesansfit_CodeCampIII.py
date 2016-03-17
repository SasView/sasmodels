from bumps.names import *

from sasmodels import core, bumps_model

if True: # fix when data loader exists
#    from sas.dataloader.readers\
    from sas.dataloader.loader import Loader
    loader = Loader()
    filename = 'sphere.ses'
    data = loader.load(filename)
    if data is None: raise IOError("Could not load file %r"%(filename,))
    data.x /= 10
#    print data
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

kernel = core.load_model("sphere", dtype='single')
phi = Parameter(0.1, name="phi")
model = bumps_model.Model(kernel,
    scale=phi*(1-phi), sld=7.0, solvent_sld=1.0, radius=radius,
    )
phi.range(0.001,0.5)
#model.radius.pmp(40)
model.radius.range(1,10000)
#model.sld.pm(5)
#model.background
#model.radius_pd=0
#model.radius_pd_n=0

### Tri-Axial Ellipsoid
#
#kernel = core.load_model("triaxial_ellipsoid", dtype='single')
#phi = Parameter(0.1, name='phi')
#model = bumps_model.Model(kernel,
#    scale=phi*(1-phi), sld=7.0, solvent_sld=1.0, radius=radius,
#    )
#phi.range(0.001,0.90)
##model.radius.pmp(40)
#model.radius.range(100,10000)
##model.sld.pmp(5)
##model.background
##model.radius_pd = 0
##model.radius_pd_n = 0

if False: # have sans data
    M_sesans = bumps_model.Experiment(data=data, model=model)
    M_sans = bumps_model.Experiment(data=sans_data, model=model)
    problem = FitProblem([M_sesans, M_sans])
else:
    M_sesans = bumps_model.Experiment(data=data, model=model)
    problem = FitProblem(M_sesans)