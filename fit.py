#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bumps.names import *
from cylcode import GpuCylinder
from lamellarcode import GpuLamellar
from ellipsecode import GpuEllipse
from coreshellcylcode import GpuCoreShellCylinder
from sasmodel import SasModel, load_data, set_beam_stop


data = load_data('JUN03289.DAT')
set_beam_stop(data, 0.004)

"""
model = SasModel(data, GpuCylinder, scale=1, radius=64.1, length=266.96, sldCyl=.291e-6, sldSolv=5.77e-6, background=0,
                              cyl_theta=0, cyl_phi=0, radius_pd=0.1, radius_pd_n=10, radius_pd_nsigma=3,length_pd=0.1,
                              length_pd_n=5, length_pd_nsigma=3, cyl_theta_pd=0.1, cyl_theta_pd_n=5, cyl_theta_pd_nsigma=3,
                              cyl_phi_pd=0.1, cyl_phi_pd_n=10, cyl_phi_pd_nsigma=3,)

model = SasModel(data, GpuEllipse, scale=.027, radius_a=60, radius_b=180, sldEll=.297e-6, sldSolv=5.773e-6, background=4.9,
                 axis_theta=0, axis_phi=90, radius_a_pd=0.1, radius_a_pd_n=10, radius_a_pd_nsigma=3, radius_b_pd=0.1, radius_b_pd_n=10,
                 radius_b_pd_nsigma=3, axis_theta_pd=0.1, axis_theta_pd_n=6, axis_theta_pd_nsigma=3, axis_phi_pd=0.1,
                 axis_phi_pd_n=6, axis_phi_pd_nsigma=3)

model = SasModel(data, GpuLamellar, scale=1, bi_thick=100, sld_bi=.291e-6, sld_sol=5.77e-6, background=0,
                 bi_thick_pd=0.1, bi_thick_pd_n=35, bi_thick_pd_nsigma=3)

"""
model = SasModel(data, GpuCoreShellCylinder, scale=1, radius=64.1, thickness=1, length=266.96, core_sld=.251e-6, shell_sld=6.2e-6,
                 solvent_sld=5.77e-6, background=0, axis_theta=0, axis_phi=0, radius_pd=0.1, radius_pd_n=10, radius_pd_nsigma=3,
                 length_pd=0.1, length_pd_n=10, length_pd_nsigma=3, thickness_pd=0.1, thickness_pd_n=2, thickness_pd_nsigma=3,
                 axis_theta_pd=0.1, axis_theta_pd_n=2, axis_theta_pd_nsigma=3, axis_phi_pd=0.1, axis_phi_pd_n=2,
                 axis_phi_pd_nsigma=3)

model.scale.range(0,10)

problem = FitProblem(model)


