#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bumps.names import *
from sasmodel import SasModel, load_data, set_beam_stop, set_half, set_top
from Models.code_capcyl import GpuCapCylinder
from Models.code_coreshellcyl import GpuCoreShellCylinder
from Models.code_cylinder import GpuCylinder, OneDGpuCylinder
from Models.code_ellipse import GpuEllipse
from Models.code_lamellar import GpuLamellar
from Models.code_triaxialellipse import GpuTriEllipse

""" IMPORT THE DATA USED """
#data = load_data('December/Tangential/Sector0/DEC07133.ABS')
data = load_data('December/DEC07102.DAT')

""" SET INNER BEAM STOP, OUTER RING, AND MASK HALF OF THE DATA """
set_beam_stop(data, 0.00669)#, outer=0.025)
set_top(data, -.018)
#set_half(data, 'right')


if 0:
    model = SasModel(data, OneDGpuCylinder,
    scale=0.0013, radius=105, length=1000,
    background=21, sldCyl=.291e-6,sldSolv=7.105e-6,
    radius_pd=0.1,radius_pd_n=10,radius_pd_nsigma=0,
    length_pd=0.1,length_pd_n=5,length_pd_nsigma=0,
    bolim=0.0, uplim=90) #bottom limit, upper limit of angle integral


if 0:
    model = SasModel(data, GpuEllipse,
    scale=0.08,
    radius_a=15, radius_b=800,
    sldEll=.291e-6, sldSolv=7.105e-6,
    background=0,
    axis_theta=90, axis_phi=0,
    axis_theta_pd=15, axis_theta_pd_n=40, axis_theta_pd_nsigma=3,
    radius_a_pd=0.222296, radius_a_pd_n=1, radius_a_pd_nsigma=0,
    radius_b_pd=.000128, radius_b_pd_n=1, radius_b_pd_nsigma=0,
    axis_phi_pd=2.63698e-05, axis_phi_pd_n=20, axis_phi_pd_nsigma=0,
    dtype='float')


    # SET THE FITTING PARAMETERS
    model.radius_a.range(15, 1000)
    model.radius_b.range(15, 1000)
    #model.axis_theta_pd.range(0, 360)
    #model.background.range(0,1000)
    #model.scale.range(0, 1)


if 0:
    model = SasModel(data, GpuLamellar,
    scale=0.08,
    bi_thick=19.2946,
    sld_bi=5.38e-6,sld_sol=7.105e-6,
    background=0.003,
    bi_thick_pd= 0.37765, bi_thick_pd_n=10, bi_thick_pd_nsigma=3,
    dtype='float')

    # SET THE FITTING PARAMETERS
    #model.bi_thick.range(0, 1000)
    #model.scale.range(0, 1)
    #model.bi_thick_pd.range(0, 1000)
    #model.background.range(0, 1000)
    model.sld_bi.range(0, 1)


if 1:
    model = SasModel(data, GpuCylinder,
    scale=0.0104, radius=92.5, length=798.3,
    sldCyl=.29e-6, sldSolv=7.105e-6, background=5,
    cyl_theta=0, cyl_phi=0,
    cyl_theta_pd=22.11, cyl_theta_pd_n=20, cyl_theta_pd_nsigma=3,
    radius_pd=.0084, radius_pd_n=10, radius_pd_nsigma=3,
    length_pd=0.493, length_pd_n=10, length_pd_nsigma=3,
    cyl_phi_pd=0, cyl_phi_pd_n=1, cyl_phi_pd_nsigma=3,
    dtype='float')



    # SET THE FITTING PARAMETERS
    #model.radius.range(1, 500)
    #model.length.range(1, 4000)
    #model.cyl_theta.range(-90,100)
    #model.cyl_theta_pd.range(0, 90)
    #model.cyl_theta_pd_n = model.cyl_theta_pd + 5
    #model.radius_pd.range(0, 90)
    #model.length_pd.range(0, 90)
    model.scale.range(0, 1)
    #model.background.range(0, 100)
    #model.sldCyl.range(0, 1)


if 0:
    model = SasModel(data, GpuCoreShellCylinder,
    scale= .00031, radius=19.5, thickness=30, length=22,
    core_sld=7.105e-6, shell_sld=.291e-6, solvent_sld=7.105e-6,
    background=0.2, axis_theta=0, axis_phi=0,

    radius_pd=0.26, radius_pd_n=10, radius_pd_nsigma=3,
    length_pd=0.26, length_pd_n=10, length_pd_nsigma=3,
    thickness_pd=1, thickness_pd_n=1, thickness_pd_nsigma=1,
    axis_theta_pd=1, axis_theta_pd_n=10, axis_theta_pd_nsigma=3,
    axis_phi_pd=0.1, axis_phi_pd_n=1, axis_phi_pd_nsigma=1,
    dtype='float')

    # SET THE FITTING PARAMETERS
    #model.radius.range(115, 1000)
    #model.length.range(0, 2500)
    #model.thickness.range(18, 38)
    #model.thickness_pd.range(0, 1)
    #model.axis_phi.range(0, 90)
    #model.radius_pd.range(0, 1)
    #model.length_pd.range(0, 1)
    #model.axis_theta_pd.range(0, 360)
    #model.background.range(0,5)
    model.scale.range(0, 1)



if 0:
    model = SasModel(data, GpuCapCylinder,
    scale=1, rad_cyl=20, rad_cap=40, length=400,
    sld_capcyl=1e-6, sld_solv=6.3e-6,
    background=0, theta=0, phi=0,
    rad_cyl_pd=.1, rad_cyl_pd_n=1, rad_cyl_pd_nsigma=0,
    rad_cap_pd=.1, rad_cap_pd_n=1, rad_cap_pd_nsigma=0,
    length_pd=.1, length_pd_n=1, length_pd_nsigma=0,
    theta_pd=.1, theta_pd_n=1, theta_pd_nsigma=0,
    phi_pd=.1, phi_pd_n=1, phi_pd_nsigma=0,
    dtype='float')

    model.scale.range(0, 1)


if 0:
    model = SasModel(data, GpuTriEllipse,
    scale=0.08, axisA=15, axisB=20, axisC=500,
    sldEll=7.105e-6, sldSolv=.291e-6,
    background=5, theta=0, phi=0, psi=0,
    theta_pd=20, theta_pd_n=40, theta_pd_nsigma=3,
    phi_pd=.1, phi_pd_n=1, phi_pd_nsigma=0,
    psi_pd=30, psi_pd_n=1, psi_pd_nsigma=0,
    axisA_pd=.1, axisA_pd_n=1, axisA_pd_nsigma=0,
    axisB_pd=.1, axisB_pd_n=1, axisB_pd_nsigma=0,
    axisC_pd=.1, axisC_pd_n=1, axisC_pd_nsigma=0, dtype='float')

    # SET THE FITTING PARAMETERS
    model.axisA.range(15, 1000)
    model.axisB.range(15, 1000)
    #model.axisC.range(15, 1000)
    #model.background.range(0,1000)
    #model.theta_pd.range(0, 360)
    #model.phi_pd.range(0, 360)
    #model.psi_pd.range(0, 360)




problem = FitProblem(model)

