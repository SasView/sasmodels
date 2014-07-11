#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bumps.names import *
from code_cylinder import GpuCylinder
from code_lamellar import GpuLamellar
from code_ellipse import GpuEllipse
from code_coreshellcyl import GpuCoreShellCylinder
from code_capcyl import GpuCapCylinder
from code_triaxialellipse import GpuTriEllipse
from sasmodel import SasModel, load_data, set_beam_stop


data = load_data('JUN03289.DAT')
set_beam_stop(data, 0.004)
"""
model = SasModel(data, GpuCylinder, scale=1, radius=64.1, length=266.96, sldCyl=.291e-6, sldSolv=5.77e-6, background=0,
                              cyl_theta=0, cyl_phi=0, radius_pd=0.1, radius_pd_n=10, radius_pd_nsigma=3,length_pd=0.1,
                              length_pd_n=5, length_pd_nsigma=3, cyl_theta_pd=0.1, cyl_theta_pd_n=5, cyl_theta_pd_nsigma=3,
                              cyl_phi_pd=0.1, cyl_phi_pd_n=10, cyl_phi_pd_nsigma=3, dtype='float')
model.radius.range(0,100)
model.length.range(0, 1000)
model.cyl_theta.range(0,90)
model.cyl_phi.range(0,90)
"""


model = SasModel(data, GpuEllipse, scale=.027, radius_a=60, radius_b=180, sldEll=.297e-6, sldSolv=5.773e-6, background=4.9,
                 axis_theta=0, axis_phi=90, radius_a_pd=0.1, radius_a_pd_n=10, radius_a_pd_nsigma=3, radius_b_pd=0.1, radius_b_pd_n=10,
                 radius_b_pd_nsigma=3, axis_theta_pd=0.1, axis_theta_pd_n=6, axis_theta_pd_nsigma=3, axis_phi_pd=0.1,
                 axis_phi_pd_n=6, axis_phi_pd_nsigma=3, dtype='float')


"""
model = SasModel(data, GpuLamellar, scale=1, bi_thick=100, sld_bi=.291e-6, sld_sol=5.77e-6, background=0,
                 bi_thick_pd=0.1, bi_thick_pd_n=35, bi_thick_pd_nsigma=3, dtype='float')
"""

"""
model = SasModel(data, GpuCoreShellCylinder, scale=1, radius=64.1, thickness=1, length=266.96, core_sld=1e-6, shell_sld=4e-6,
                 solvent_sld=1e-6, background=0, axis_theta=0, axis_phi=0, radius_pd=0.1, radius_pd_n=10, radius_pd_nsigma=3,
                 length_pd=0.1, length_pd_n=10, length_pd_nsigma=3, thickness_pd=0.1, thickness_pd_n=2, thickness_pd_nsigma=3,
                 axis_theta_pd=0.1, axis_theta_pd_n=2, axis_theta_pd_nsigma=3, axis_phi_pd=0.1, axis_phi_pd_n=2,
                 axis_phi_pd_nsigma=3, dtype='float')
"""


"""
model = SasModel(data, GpuCapCylinder, scale=1, rad_cyl=20, rad_cap=40, length=400, sld_capcyl=1e-6, sld_solv=6.3e-6,
                 background=0, theta=0, phi=0, rad_cyl_pd=.1, rad_cyl_pd_n=10, rad_cyl_nsigma=3, rad_cap_pd=.1, rad_cap_pd_n=1,
                 rad_cap_pd_nsigma=3, length_pd=.1, length_pd_n=10, length_pd_nsigma=3, theta_pd=.1, theta_pd_n=4,
                 theta_pd_nsigma=3, phi_pd=.1, phi_pd_n=4, phi_pd_nsigma=3, dtype='float')
"""

"""
model = SasModel(data, GpuTriEllipse, scale=1, axisA=35, axisB=100, axisC=400, sldEll=1e-6, sldSolv=6.3e-6,
                 background=0, theta=57, phi=57, psi=0, theta_pd=.1, theta_pd_n=4,
                 theta_pd_nsigma=3, phi_pd=.1, phi_pd_n=4, phi_pd_nsigma=3, psi_pd=.1, psi_pd_n=4, psi_pd_nsigma=3,
                 axisA_pd=.1, axisA_pd_n=4, axisA_pd_nsigma=3, axisB_pd=.1, axisB_pd_n=4, axisB_pd_nsigma=3,
                 axisC_pd=.1, axisC_pd_n=4, axisC_pd_nsigma=3, dtype='float')
"""


model.scale.range(0,1)

problem = FitProblem(model)


