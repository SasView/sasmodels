#!/usr/bin/env python
# -*- coding: utf-8 -*-

import datetime
from sasmodel import SasModel, load_data, set_beam_stop, plot_data

TIC = None
def tic():
    global TIC
    TIC = datetime.datetime.now()

def toc():
    now = datetime.datetime.now()
    return (now-TIC).total_seconds()

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

def cyl(N=1):
    import sys
    import matplotlib.pyplot as plt

    if len(sys.argv) > 1:
        N = int(sys.argv[1])
    data = load_data('JUN03289.DAT')
    set_beam_stop(data, 0.004)

    pars = dict(
        scale=1, radius=64.1, length=266.96, sldCyl=.291e-6, sldSolv=5.77e-6, background=0,
        cyl_theta=0, cyl_phi=0, radius_pd=0.1, radius_pd_n=10, radius_pd_nsigma=3,length_pd=0.1,
        length_pd_n=5, length_pd_nsigma=3, cyl_theta_pd=0.1, cyl_theta_pd_n=5, cyl_theta_pd_nsigma=3,
        cyl_phi_pd=0.1, cyl_phi_pd_n=10, cyl_phi_pd_nsigma=3,
        )

    model = sasview_model('Cylinder', **pars)
    tic()
    for i in range(N):
        cpu = sasview_eval(model, data)
    cpu_time = toc()*1000./N

    from code_cylinder import GpuCylinder
    model = SasModel(data, GpuCylinder, dtype='f', **pars)
    tic()
    for i in range(N):
        gpu = model.theory()
    gpu_time = toc()*1000./N

    relerr = (gpu - cpu)/cpu
    print "max(|(ocl-omp)/ocl|)", max(abs(relerr))
    print "omp t=%.1f ms"%cpu_time
    print "ocl t=%.1f ms"%gpu_time

    plt.subplot(131); plot_data(data, cpu); plt.title("omp t=%.1f ms"%cpu_time)
    plt.subplot(132); plot_data(data, gpu); plt.title("ocl t=%.1f ms"%gpu_time)
    plt.subplot(133); plot_data(data, 1e8*relerr); plt.title("relerr x 10^8"); plt.colorbar()
    plt.show()

def ellipse(N=4):
    import sys
    import matplotlib.pyplot as plt

    if len(sys.argv) > 1:
        N = int(sys.argv[1])
    data = load_data('DEC07133.DAT')
    set_beam_stop(data, 0.004)

    pars = dict(scale=.027, radius_a=60, radius_b=180, sldEll=.297e-6, sldSolv=5.773e-6, background=4.9,
                axis_theta=0, axis_phi=90, radius_a_pd=0.1, radius_a_pd_n=10, radius_a_pd_nsigma=3, radius_b_pd=0.1, radius_b_pd_n=10,
                radius_b_pd_nsigma=3, axis_theta_pd=0.1, axis_theta_pd_n=6, axis_theta_pd_nsigma=3, axis_phi_pd=0.1,
                axis_phi_pd_n=6, axis_phi_pd_nsigma=3,)

    model = sasview_model('Ellipsoid', **pars)
    tic()
    for i in range(N):
        cpu = sasview_eval(model, data)
    cpu_time = toc()*1000./N

    from code_ellipse import GpuEllipse
    model = SasModel(data, GpuEllipse, dtype='f', **pars)
    tic()
    for i in range(N):
        gpu = model.theory()
    gpu_time = toc()*1000./N

    relerr = (gpu - cpu)/cpu
    print "max(|(ocl-omp)/ocl|)", max(abs(relerr))
    print "omp t=%.1f ms"%cpu_time
    print "ocl t=%.1f ms"%gpu_time

    plt.subplot(131); plot_data(data, cpu); plt.title("omp t=%.1f ms"%cpu_time)
    plt.subplot(132); plot_data(data, gpu); plt.title("ocl t=%.1f ms"%gpu_time)
    plt.subplot(133); plot_data(data, 1e8*relerr); plt.title("relerr x 10^8"); plt.colorbar()
    plt.show()

def coreshell(N=4):
    import sys
    import matplotlib.pyplot as plt

    if len(sys.argv) > 1:
        N = int(sys.argv[1])
    data = load_data('DEC07133.DAT')
    set_beam_stop(data, 0.004)

    pars = dict(scale= 1.77881e-06, radius=325, thickness=25, length=34.2709,
                 core_sld=1e-6, shell_sld=.291e-6, solvent_sld=7.105e-6,
                 background=223.827, axis_theta=90, axis_phi=0,
                 axis_theta_pd=15.8,
                 radius_pd=0.1, radius_pd_n=1, radius_pd_nsigma=0,
                 length_pd=0.1, length_pd_n=1, length_pd_nsigma=0,
                 thickness_pd=0.1, thickness_pd_n=1, thickness_pd_nsigma=0,
                 axis_theta_pd_n=20, axis_theta_pd_nsigma=3,
                 axis_phi_pd=0.0008748, axis_phi_pd_n=60, axis_phi_pd_nsigma=3,)

    model = sasview_model('CoreShellCylinder', **pars)
    tic()
    for i in range(N):
        cpu = sasview_eval(model, data)
    cpu_time = toc()*1000./N

    from code_coreshellcyl import GpuCoreShellCylinder
    model = SasModel(data, GpuCoreShellCylinder, dtype='f', **pars)
    tic()
    for i in range(N):
        gpu = model.theory()
    gpu_time = toc()*1000./N

    relerr = (gpu - cpu)/cpu
    print "max(|(ocl-omp)/ocl|)", max(abs(relerr))
    print "omp t=%.1f ms"%cpu_time
    print "ocl t=%.1f ms"%gpu_time

    plt.subplot(131); plot_data(data, cpu); plt.title("omp t=%.1f ms"%cpu_time)
    plt.subplot(132); plot_data(data, gpu); plt.title("ocl t=%.1f ms"%gpu_time)
    plt.subplot(133); plot_data(data, 1e8*relerr); plt.title("relerr x 10^8"); plt.colorbar()
    plt.show()

if __name__ == "__main__":
    coreshell()





























