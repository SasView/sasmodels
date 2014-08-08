#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sasmodel import SasModel, load_data, set_beam_stop, plot_data, tic, toc
import numpy as np


def sasview_model(modelname, **pars):
    modelname = modelname+"Model"
    sans = __import__('sans.models.'+modelname)
    ModelClass = getattr(getattr(sans.models,modelname,None),modelname,None)
    if ModelClass is None:
        raise ValueError("could not find model %r in sans.models"%modelname)
    model = ModelClass()

    for k,v in pars.items():
        if k.endswith("_pd"):
            model.dispersion[k[:-3]]['width'] = v
        elif k.endswith("_pd_n"):
            model.dispersion[k[:-5]]['npts'] = v
        elif k.endswith("_pd_nsigma"):
            model.dispersion[k[:-10]]['nsigmas'] = v
        else:
            model.setParam(k, v)
    return model


def sasview_eval(model, data):
    theory = model.evalDistribution([data.qx_data, data.qy_data])
    return theory

def plot(data, cpumodel, gpumodel, N, pars):

    model = SasModel(data, gpumodel, dtype='f', **pars)
    tic()
    for i in range(N):
        #pars['scale'] = np.random.rand()
        model.update()
        gpu = model.theory()
    gpu_time = toc()*1000./N
    print "ocl t=%.1f ms"%gpu_time
    tic()
    for i in range(N):
        cpu = sasview_eval(cpumodel, data)
    cpu_time = toc()*1000./N
    print "omp t=%.1f ms"%cpu_time

    import matplotlib.pyplot as plt

    print "gpu/cpu", max(gpu/cpu)
    cpu *= max(gpu/cpu)
    relerr = (gpu - cpu)/cpu
    print "max(|(ocl-omp)/ocl|)", max(abs(relerr))


    plt.subplot(131); plot_data(data, cpu); plt.title("omp t=%.1f ms"%cpu_time)
    plt.subplot(132); plot_data(data, gpu); plt.title("ocl t=%.1f ms"%gpu_time)
    plt.subplot(133); plot_data(data, relerr); plt.title("relerr x 10^8"); plt.colorbar()
    plt.show()

def newlen(N):
    import sys

    if len(sys.argv) > 1:
        N = int(sys.argv[1])
    data = load_data('December/DEC07098.DAT')
    set_beam_stop(data, 0.004)

    return data, N

def cyl(N=1):

    dtype = "float"
    pars = dict(
        scale=.003, radius=64.1, length=66.96, sldCyl=.291e-6, sldSolv=5.77e-6, background=.1,
        cyl_theta=90, cyl_phi=0,
        radius_pd=0.1, radius_pd_n=10, radius_pd_nsigma=3,
        length_pd=0.1,length_pd_n=5, length_pd_nsigma=3,
        cyl_theta_pd=0.1, cyl_theta_pd_n=5, cyl_theta_pd_nsigma=3,
        cyl_phi_pd=0.1, cyl_phi_pd_n=10, cyl_phi_pd_nsigma=3)

    from Models.code_cylinder_f import GpuCylinder as gpumodel
    model = sasview_model('Cylinder', **pars)
    data, N = newlen(N)

    plot(data, model, gpumodel, N, pars)


def ellipse(N=1):

    pars = dict(scale=.027, radius_a=60, radius_b=180, sldEll=.297e-6, sldSolv=5.773e-6, background=4.9,
                axis_theta=0, axis_phi=90,
                radius_a_pd=0.1, radius_a_pd_n=10, radius_a_pd_nsigma=3,
                radius_b_pd=0.1, radius_b_pd_n=10, radius_b_pd_nsigma=3,
                axis_theta_pd=0.1, axis_theta_pd_n=6, axis_theta_pd_nsigma=3,
                axis_phi_pd=0.1, axis_phi_pd_n=6, axis_phi_pd_nsigma=3,)

    from Models.code_ellipse import GpuEllipse as gpumodel
    model = sasview_model('Ellipsoid', **pars)

    data, N = newlen(N)
    plot(data, model, gpumodel, N, pars)

def coreshell(N=1):

    data, N = newlen(N)

    pars = dict(scale= 1.77881e-06, radius=325, thickness=25, length=34.2709,
                 core_sld=1e-6, shell_sld=.291e-6, solvent_sld=7.105e-6,
                 background=223.827, axis_theta=90, axis_phi=0,
                 axis_theta_pd=15.8,
                 radius_pd=0.1, radius_pd_n=10, radius_pd_nsigma=3,
                 length_pd=0.1, length_pd_n=10, length_pd_nsigma=3,
                 thickness_pd=0.1, thickness_pd_n=5, thickness_pd_nsigma=3,
                 axis_theta_pd_n=20, axis_theta_pd_nsigma=5,
                 axis_phi_pd=0.0008748, axis_phi_pd_n=5, axis_phi_pd_nsigma=3,)

    model = sasview_model('CoreShellCylinder', **pars)
    from Models.code_coreshellcyl import GpuCoreShellCylinder as gpumodel

    plot(data, model, gpumodel, N, pars)

def trellipse(N=1):

    data, N = newlen(N)

    pars = dict(scale=0.08, semi_axisA=15, semi_axisB=20, semi_axisC=500,
                sldEll=7.105e-6, sldSolv=.291e-6,
                background=5, axis_theta=0, axis_phi=0, axis_psi=0,
                axis_theta_pd=20, axis_theta_pd_n=10, axis_theta_pd_nsigma=3,
                axis_phi_pd=.1, axis_phi_pd_n=10, axis_phi_pd_nsigma=3,
                axis_psi_pd=30, axis_psi_pd_n=5, axis_psi_pd_nsigma=3,
                semi_axisA_pd=.1, semi_axisA_pd_n=5, semi_axisA_pd_nsigma=3,
                semi_axisB_pd=.1, semi_axisB_pd_n=5, semi_axisB_pd_nsigma=3,
                semi_axisC_pd=.1, semi_axisC_pd_n=5, semi_axisC_pd_nsigma=3,)

    model = sasview_model('TriaxialEllipsoid', **pars)
    from Models.code_triaxialellipse import GpuTriEllipse as gpumodel

    plot(data, model, gpumodel, N, pars)

def lamellar(N=1):

    data, N = newlen(N)

    pars = dict(scale=0.08,
                bi_thick=19.2946,
                sld_bi=5.38e-6,sld_sol=7.105e-6,
                background=0.003,
                bi_thick_pd= 0.37765, bi_thick_pd_n=40, bi_thick_pd_nsigma=3)

    model = sasview_model('Lamellar', **pars)
    from Models.code_lamellar import GpuLamellar as gpumodel

    plot(data, model, gpumodel, N, pars)

def cap(N=1):

    data, N = newlen(N)

    pars = dict(scale=.08, rad_cyl=20, rad_cap=40, len_cyl=400,
                sld_capcyl=1e-6, sld_solv=6.3e-6,
                background=0, theta=0, phi=0,
                rad_cyl_pd=.1, rad_cyl_pd_n=10, rad_cyl_pd_nsigma=3,
                rad_cap_pd=.1, rad_cap_pd_n=10, rad_cap_pd_nsigma=3,
                len_cyl_pd=.1, len_cyl_pd_n=3, len_cyl_pd_nsigma=3,
                theta_pd=.1, theta_pd_n=3, theta_pd_nsigma=3,
                phi_pd=.1, phi_pd_n=3, phi_pd_nsigma=3)


    model = sasview_model('CappedCylinder', **pars)
    from Models.code_capcyl import GpuCapCylinder as gpumodel

    plot(data, model, gpumodel, N, pars)


if __name__ == "__main__":
    cyl(2)





























