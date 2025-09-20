from pathlib import Path

from bumps.names import FitProblem, Parameter

from sasmodels import bumps_model, core
from sasmodels.data import load_data

path = Path(__file__).resolve().parent
data = load_data(str(path / 'sphere.ses'))

radius = 1000
data.Rmax = 3*radius # [A]

##  Sphere parameters

kernel = core.load_model("sphere")
phi = Parameter(0.1, name="phi")
model = bumps_model.Model(
    kernel,
    scale=phi*(1-phi), sld=7.0, sld_solvent=1.0, radius=radius,
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


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    problem.plot()
    plt.show()