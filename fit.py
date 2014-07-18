#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bumps.names import *
from code_cylinder import GpuCylinder
from code_lamellar import GpuLamellar
from code_ellipse import GpuEllipse
from code_coreshellcyl import GpuCoreShellCylinder
from code_capcyl import GpuCapCylinder
from code_triaxialellipse import GpuTriEllipse
from sasmodel import SasModel, load_data, set_beam_stop, set_half


data = load_data('DEC07282.DAT')
set_beam_stop(data, 0.0048)#, outer=0.025)
#set_half(data, 'right')
"""
model = SasModel(data, GpuCylinder, scale=0.09402, radius=106, length=442.298,
                                    sldCyl=.291e-6, sldSolv=5.77e-6, background=0,
                                    cyl_theta=0, cyl_phi=0,
                                    radius_pd=0.1, radius_pd_n=10, radius_pd_nsigma=0,
                                    length_pd=0.1, length_pd_n=5, length_pd_nsigma=0,
                                    cyl_theta_pd=0.1, cyl_theta_pd_n=10, cyl_theta_pd_nsigma=3,
                                    cyl_phi_pd=0.1, cyl_phi_pd_n=4, cyl_phi_pd_nsigma=3,
                                    dtype='float')
model.radius.range(0, 1000)
model.length.range(0, 1000)
#model.cyl_theta.range(0,90)
#model.cyl_phi.range(0,90)
model.scale.range(0, 1)
"""
"""
model = SasModel(data, GpuEllipse,
                 scale=0.0007, radius_a=42.8251, radius_b=700, sldEll=.291e-6, sldSolv=7.105e-6,
                 background=10, axis_theta=0, axis_phi=0,
                 radius_a_pd=0.222296, radius_a_pd_n=1, radius_a_pd_nsigma=0,
                 radius_b_pd=.000128, radius_b_pd_n=1, radius_b_pd_nsigma=0,
                 axis_theta_pd=24.8059, axis_theta_pd_n=40, axis_theta_pd_nsigma=3,
                 axis_phi_pd=2.63698e-05, axis_phi_pd_n=20, axis_phi_pd_nsigma=3,
                 dtype='float')

model.radius_a.range(15, 1000)
#model.radius_b.range(15, 1000)
#model.axis_theta_pd.range(0, 360)
#model.background.range(0,1000)
model.scale.range(0, 1)
"""
"""
model = SasModel(data, GpuLamellar, scale=0.0128, bi_thick=35.1, sld_bi=.291e-6, sld_sol=5.77e-6, background=17.3,
                 bi_thick_pd= 0.000009, bi_thick_pd_n=35, bi_thick_pd_nsigma=3, dtype='float')
model.bi_thick.range(0, 1000)
model.scale.range(0, 1)
#model.bi_thick_pd.range(0, 1000)
#model.background.range(0, 1000)
"""
"""
model = SasModel(data, GpuCoreShellCylinder,
                 scale= 1.77881e-06, radius=325, thickness=25, length=34.2709,
                 core_sld=1e-6, shell_sld=.291e-6, solvent_sld=7.105e-6,
                 background=223.827, axis_theta=90, axis_phi=0,
                 axis_theta_pd=15.8,

                 radius_pd=0.1, radius_pd_n=1, radius_pd_nsigma=0,
                 length_pd=0.1, length_pd_n=1, length_pd_nsigma=0,
                 thickness_pd=0.1, thickness_pd_n=1, thickness_pd_nsigma=0,
                 axis_theta_pd_n=10, axis_theta_pd_nsigma=3,
                 axis_phi_pd=0.0008748, axis_phi_pd_n=10, axis_phi_pd_nsigma=3,
                 dtype='float')

model.radius.range(15, 1000)
#model.length.range(0, 1000)
#model.thickness.range(20, 50)
#model.axis_phi.range(0, 90)
#model.radius_pd.range(0, 1)
#model.radius_b_pd.range(0, 1)
#model.axis_theta_pd.range(0, 360)
#model.axis_phi_pd.range(0, 360)
#model.background.range(0,1000)
model.scale.range(0, 1)
"""


"""
model = SasModel(data, GpuCapCylinder, scale=1, rad_cyl=20, rad_cap=40, length=400, sld_capcyl=1e-6, sld_solv=6.3e-6,
                 background=0, theta=0, phi=0, rad_cyl_pd=.1, rad_cyl_pd_n=10, rad_cyl_nsigma=3, rad_cap_pd=.1, rad_cap_pd_n=1,
                 rad_cap_pd_nsigma=3, length_pd=.1, length_pd_n=10, length_pd_nsigma=3, theta_pd=.1, theta_pd_n=4,
                 theta_pd_nsigma=3, phi_pd=.1, phi_pd_n=4, phi_pd_nsigma=3, dtype='float')
"""


model = SasModel(data, GpuTriEllipse,
                 scale=0.0036, axisA=118, axisB=70, axisC=800,
                 sldEll=7.105e-6, sldSolv=.291e-6,
                 background=15, theta=90, phi=0, psi=0,
                 theta_pd=22, theta_pd_n=40, theta_pd_nsigma=3,
                 phi_pd=.1, phi_pd_n=1, phi_pd_nsigma=0,
                 psi_pd=30, psi_pd_n=1, psi_pd_nsigma=0,
                 axisA_pd=.1, axisA_pd_n=1, axisA_pd_nsigma=0,
                 axisB_pd=.1, axisB_pd_n=1, axisB_pd_nsigma=0,
                 axisC_pd=.1, axisC_pd_n=1, axisC_pd_nsigma=0, dtype='float')

#model.axisA.range(15, 1000)
#model.axisB.range(15, 1000)
#model.axisC.range(15, 1000)
#model.background.range(0,1000)
model.scale.range(0, 1)
#model.theta_pd.range(0, 360)
#model.phi_pd.range(0, 360)
#model.psi_pd.range(0, 360)





problem = FitProblem(model)


