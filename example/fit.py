#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from bumps.names import *
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasmodels.data import load_data, set_beam_stop, set_top

""" IMPORT THE DATA USED """
radial_data = load_data('DEC07267.DAT')
set_beam_stop(radial_data, 0.00669, outer=0.025)
set_top(radial_data, -.0185)

tan_data = load_data('DEC07266.DAT')
set_beam_stop(tan_data, 0.00669, outer=0.025)
set_top(tan_data, -.0185)
#sas.set_half(tan_data, 'right')

name = "ellipsoid" if len(sys.argv) < 2 else sys.argv[1]
section = "radial" if len(sys.argv) < 3 else sys.argv[2]
if section not in ("radial","tangential","both"):
    raise ValueError("section %r should be 'radial', 'tangential' or 'both'"
            % section)
data = radial_data if section != "tangential" else tan_data
theta = 89.9 if section != "tangential" else 0
phi = 90
kernel = load_model(name, dtype="single")
cutoff = 1e-3

if name == "ellipsoid":
    pars = dict(
        scale=0.08, background=35,
        radius_polar=15, radius_equatorial=800,
        sld=.291, sld_solvent=7.105,
        theta=theta, phi=phi,
        theta_pd=0, theta_pd_n=0, theta_pd_nsigma=3,
        phi_pd=0, phi_pd_n=20, phi_pd_nsigma=3,
        radius_polar_pd=0.222296, radius_polar_pd_n=1, radius_polar_pd_nsigma=0,
        radius_equatorial_pd=.000128, radius_equatorial_pd_n=1,
        radius_equatorial_pd_nsigma=0,
        )
    model = Model(kernel, **pars)

    # SET THE FITTING PARAMETERS
    model.radius_polar.range(15, 1000)
    model.radius_equatorial.range(15, 1000)
    #model.theta.range(0, 90)
    #model.theta_pd.range(0,10)
    model.phi_pd.range(0,20)
    model.phi.range(0, 180)
    model.background.range(0,1000)
    model.scale.range(0, 10)


elif name == "lamellar":
    pars = dict(
        scale=0.08, background=0.003,
        thickness=19.2946,
        sld=5.38,sld_sol=7.105,
        thickness_pd= 0.37765, thickness_pd_n=10, thickness_pd_nsigma=3,
        )
    model = Model(kernel, **pars)

    # SET THE FITTING PARAMETERS
    #model.thickness.range(0, 1000)
    #model.scale.range(0, 1)
    #model.thickness_pd.range(0, 1000)
    #model.background.range(0, 1000)
    model.sld.range(0, 1)


elif name == "cylinder":
    pars = dict(
        scale=.01, background=35,
        sld=.291, sld_solvent=5.77,
        radius=250, length=178,
        radius_pd=0.1, radius_pd_n=5, radius_pd_nsigma=3,
        length_pd=0.1,length_pd_n=5, length_pd_nsigma=3,
        theta=theta, phi=phi,
        theta_pd=0, theta_pd_n=0, theta_pd_nsigma=3,
        phi_pd=10, phi_pd_n=20, phi_pd_nsigma=3)
    model = Model(kernel, **pars)

    # SET THE FITTING PARAMETERS
    model.radius.range(1, 500)
    model.length.range(1, 5000)
    #model.theta.range(0, 90)
    model.phi.range(0, 180)
    model.phi_pd.range(0, 30)
    model.radius_pd.range(0, 1)
    model.length_pd.range(0, 1)
    model.scale.range(0, 10)
    model.background.range(0, 100)


elif name == "core_shell_cylinder":
    model = Model(kernel,
        scale= .031, background=0,
        radius=19.5, thickness=30, length=22,
        sld_core=7.105, sld_shell=.291, sld_solvent=7.105,
        radius_pd=0.26, radius_pd_n=10, radius_pd_nsigma=3,
        length_pd=0.26, length_pd_n=10, length_pd_nsigma=3,
        thickness_pd=1, thickness_pd_n=1, thickness_pd_nsigma=1,
        theta=theta, phi=phi,
        theta_pd=1, theta_pd_n=1, theta_pd_nsigma=3,
        phi_pd=0, phi_pd_n=20, phi_pd_nsigma=3,
        )

    # SET THE FITTING PARAMETERS
    model.radius.range(115, 1000)
    model.length.range(0, 2500)
    #model.thickness.range(18, 38)
    #model.thickness_pd.range(0, 1)
    #model.phi.range(0, 90)
    model.phi_pd.range(0,20)
    #model.radius_pd.range(0, 1)
    #model.length_pd.range(0, 1)
    #model.theta_pd.range(0, 360)
    #model.background.range(0,5)
    model.scale.range(0, 1)



elif name == "capped_cylinder":
    model = Model(kernel,
        scale=.08, background=35,
        radius=20, cap_radius=40, length=400,
        sld=1, sld_solvent=6.3,
        radius_pd=.1, radius_pd_n=5, radius_pd_nsigma=3,
        cap_radius_pd=.1, cap_radius_pd_n=5, cap_radius_pd_nsigma=3,
        length_pd=.1, length_pd_n=1, length_pd_nsigma=0,
        theta=theta, phi=phi,
        theta_pd=0, theta_pd_n=1, theta_pd_nsigma=0,
        phi_pd=10, phi_pd_n=20, phi_pd_nsigma=0,
        )

    model.radius.range(115, 1000)
    model.length.range(0, 2500)
    #model.thickness.range(18, 38)
    #model.thickness_pd.range(0, 1)
    #model.phi.range(0, 90)
    model.phi_pd.range(0,20)
    #model.radius_pd.range(0, 1)
    #model.length_pd.range(0, 1)
    #model.theta_pd.range(0, 360)
    #model.background.range(0,5)
    model.scale.range(0, 1)


elif name == "triaxial_ellipsoid":
    pars = dict(
        scale=0.08, background=35,
        radius_equat_minor=15, radius_equat_major=20, radius_polar=500,
        sld=7.105, solvent_sld=.291,
        radius_equat_minor_pd=.1, radius_equat_minor_pd_n=1, radius_equat_minor_pd_nsigma=0,
        radius_equat_major_pd=.1, radius_equat_major_pd_n=1, radius_equat_major_pd_nsigma=0,
        radius_polar_pd=.1, radius_polar_pd_n=1, radius_polar_pd_nsigma=0,
        theta=theta, phi=phi, psi=0,
        theta_pd=20, theta_pd_n=40, theta_pd_nsigma=3,
        phi_pd=.1, phi_pd_n=1, phi_pd_nsigma=0,
        psi_pd=30, psi_pd_n=1, psi_pd_nsigma=0,
        )
    model = Model(kernel, **pars)

    # SET THE FITTING PARAMETERS
    model.radius_equat_minor.range(15, 1000)
    model.radius_equat_major.range(15, 1000)
    #model.radius_polar.range(15, 1000)
    #model.background.range(0,1000)
    #model.theta_pd.range(0, 360)
    #model.phi_pd.range(0, 360)
    #model.psi_pd.range(0, 360)

else:
    print("No parameters for %s"%name)
    sys.exit(1)

model.cutoff = cutoff
M = Experiment(data=data, model=model)
if section == "both":
    tan_model = Model(model.sasmodel, **model.parameters())
    tan_model.phi = model.phi - 90
    tan_model.cutoff = cutoff
    tan_M = Experiment(data=tan_data, model=tan_model)
    problem = FitProblem([M, tan_M])
else:
    problem = FitProblem(M)

if __name__ == "__main__":
    problem.plot()
    import pylab; pylab.show()
