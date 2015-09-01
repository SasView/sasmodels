#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from bumps.names import *
from sasmodels.core import load_model
from sasmodels import bumps_model as sas

""" IMPORT THE DATA USED """
radial_data = sas.load_data('DEC07267.DAT')
sas.set_beam_stop(radial_data, 0.00669, outer=0.025)
sas.set_top(radial_data, -.0185)

tan_data = sas.load_data('DEC07266.DAT')
sas.set_beam_stop(tan_data, 0.00669, outer=0.025)
sas.set_top(tan_data, -.0185)
#sas.set_half(tan_data, 'right')

name = "ellipsoid" if len(sys.argv) < 2 else sys.argv[1]
section = "radial" if len(sys.argv) < 3 else sys.argv[2]
if section not in ("radial","tangential","both"):
    raise ValueError("section %r should be 'radial', 'tangential' or 'both'"
            % section)
data = radial_data if section != "tangential" else tan_data
phi = 0 if section != "tangential" else 90
kernel = load_model(name, dtype="single")
cutoff = 1e-3

if name == "ellipsoid":
    model = sas.Model(kernel,
        scale=0.08,
        rpolar=15, requatorial=800,
        sld=.291, solvent_sld=7.105,
        background=0,
        theta=90, phi=phi,
        theta_pd=15, theta_pd_n=40, theta_pd_nsigma=3,
        rpolar_pd=0.222296, rpolar_pd_n=1, rpolar_pd_nsigma=0,
        requatorial_pd=.000128, requatorial_pd_n=1, requatorial_pd_nsigma=0,
        phi_pd=0, phi_pd_n=20, phi_pd_nsigma=3,
        )


    # SET THE FITTING PARAMETERS
    model.rpolar.range(15, 1000)
    model.requatorial.range(15, 1000)
    model.theta_pd.range(0, 360)
    model.background.range(0,1000)
    model.scale.range(0, 10)



elif name == "lamellar":
    model = sas.Model(kernel,
        scale=0.08,
        thickness=19.2946,
        sld=5.38,sld_sol=7.105,
        background=0.003,
        thickness_pd= 0.37765, thickness_pd_n=10, thickness_pd_nsigma=3,
        )


    # SET THE FITTING PARAMETERS
    #model.thickness.range(0, 1000)
    #model.scale.range(0, 1)
    #model.thickness_pd.range(0, 1000)
    #model.background.range(0, 1000)
    model.sld.range(0, 1)


elif name == "cylinder":
    """
    pars = dict(scale=0.0023, radius=92.5, length=798.3,
        sld=.29, solvent_sld=7.105, background=5,
        theta=0, phi=phi,
        theta_pd=22.11, theta_pd_n=5, theta_pd_nsigma=3,
        radius_pd=.0084, radius_pd_n=10, radius_pd_nsigma=3,
        length_pd=0.493, length_pd_n=10, length_pd_nsigma=3,
        phi_pd=0, phi_pd_n=5, phi_pd_nsigma=3,)
        """
    pars = dict(
        scale=.01, background=35,
        sld=.291, solvent_sld=5.77, 
        radius=250, length=178, 
        theta=90, phi=phi,
        radius_pd=0.1, radius_pd_n=5, radius_pd_nsigma=3,
        length_pd=0.1,length_pd_n=5, length_pd_nsigma=3,
        theta_pd=10, theta_pd_n=50, theta_pd_nsigma=3,
        phi_pd=0, phi_pd_n=10, phi_pd_nsigma=3)
    model = sas.Model(kernel, **pars)

    # SET THE FITTING PARAMETERS
    model.radius.range(1, 500)
    model.length.range(1, 5000)
    model.theta.range(-90,100)
    model.theta_pd.range(0, 30)
    model.theta_pd_n = model.theta_pd + 5
    model.radius_pd.range(0, 1)
    model.length_pd.range(0, 2)
    model.scale.range(0, 10)
    model.background.range(0, 100)


elif name == "core_shell_cylinder":
    model = sas.Model(kernel,
        scale= .031, radius=19.5, thickness=30, length=22,
        core_sld=7.105, shell_sld=.291, solvent_sld=7.105,
        background=0, theta=0, phi=phi,

        radius_pd=0.26, radius_pd_n=10, radius_pd_nsigma=3,
        length_pd=0.26, length_pd_n=10, length_pd_nsigma=3,
        thickness_pd=1, thickness_pd_n=1, thickness_pd_nsigma=1,
        theta_pd=1, theta_pd_n=10, theta_pd_nsigma=3,
        phi_pd=0.1, phi_pd_n=1, phi_pd_nsigma=1,
        )

    # SET THE FITTING PARAMETERS
    #model.radius.range(115, 1000)
    #model.length.range(0, 2500)
    #model.thickness.range(18, 38)
    #model.thickness_pd.range(0, 1)
    #model.phi.range(0, 90)
    #model.radius_pd.range(0, 1)
    #model.length_pd.range(0, 1)
    #model.theta_pd.range(0, 360)
    #model.background.range(0,5)
    model.scale.range(0, 1)



elif name == "capped_cylinder":
    model = sas.Model(kernel,
        scale=.08, radius=20, cap_radius=40, length=400,
        sld_capcyl=1, sld_solv=6.3,
        background=0, theta=0, phi=phi,
        radius_pd=.1, radius_pd_n=5, radius_pd_nsigma=3,
        cap_radius_pd=.1, cap_radius_pd_n=5, cap_radius_pd_nsigma=3,
        length_pd=.1, length_pd_n=1, length_pd_nsigma=0,
        theta_pd=.1, theta_pd_n=1, theta_pd_nsigma=0,
        phi_pd=.1, phi_pd_n=1, phi_pd_nsigma=0,
        )

    model.scale.range(0, 1)


elif name == "triaxial_ellipsoid":
    model = sas.Model(kernel,
        scale=0.08, req_minor=15, req_major=20, rpolar=500,
        sldEll=7.105, solvent_sld=.291,
        background=5, theta=0, phi=phi, psi=0,
        theta_pd=20, theta_pd_n=40, theta_pd_nsigma=3,
        phi_pd=.1, phi_pd_n=1, phi_pd_nsigma=0,
        psi_pd=30, psi_pd_n=1, psi_pd_nsigma=0,
        req_minor_pd=.1, req_minor_pd_n=1, req_minor_pd_nsigma=0,
        req_major_pd=.1, req_major_pd_n=1, req_major_pd_nsigma=0,
        rpolar_pd=.1, rpolar_pd_n=1, rpolar_pd_nsigma=0,
        )

    # SET THE FITTING PARAMETERS
    model.req_minor.range(15, 1000)
    model.req_major.range(15, 1000)
    #model.rpolar.range(15, 1000)
    #model.background.range(0,1000)
    #model.theta_pd.range(0, 360)
    #model.phi_pd.range(0, 360)
    #model.psi_pd.range(0, 360)

else:
    print "No parameters for %s"%name
    sys.exit(1)

model.cutoff = cutoff
M = sas.Experiment(data=data, model=model)
if section == "both":
   tan_model = sas.Model(model.kernel, **model.parameters())
   tan_model.phi = model.phi - 90
   tan_model.cutoff = cutoff
   tan_M = sas.Experiment(data=tan_data, model=tan_model)
   problem = FitProblem([M, tan_M])
else:
   problem = FitProblem(M)
