#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import numpy as np

from sasmodels.core import BumpsModel, plot_data, tic, opencl_model, dll_model

def sasview_model(modelname, **pars):
    """
    Load a sasview model given the model name.
    """
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
        elif k.endswith("_pd_type"):
            model.dispersion[k[:-8]]['type'] = v
        else:
            model.setParam(k, v)
    return model

def load_opencl(modelname, dtype='single'):
    sasmodels = __import__('sasmodels.models.'+modelname)
    module = getattr(sasmodels.models, modelname, None)
    kernel = opencl_model(module, dtype=dtype)
    return kernel

def load_dll(modelname, dtype='single'):
    sasmodels = __import__('sasmodels.models.'+modelname)
    module = getattr(sasmodels.models, modelname, None)
    kernel = dll_model(module, dtype=dtype)
    return kernel


def compare(Ncpu, cpuname, cpupars, Ngpu, gpuname, gpupars):

    #from sasmodels.core import load_data
    #data = load_data('December/DEC07098.DAT')
    from sasmodels.core import empty_data1D
    data = empty_data1D(np.logspace(-4, -1, 128))
    #from sasmodels.core import empty_2D, set_beam_stop
    #data = empty_data2D(np.linspace(-0.05, 0.05, 128))
    #set_beam_stop(data, 0.004)
    is2D = hasattr(data, 'qx_data')

    if Ngpu > 0:
        gpumodel = load_opencl(gpuname, dtype='single')
        model = BumpsModel(data, gpumodel, **gpupars)
        toc = tic()
        for i in range(Ngpu):
            #pars['scale'] = np.random.rand()
            model.update()
            gpu = model.theory()
        gpu_time = toc()*1000./Ngpu
        print "ocl t=%.1f ms, intensity=%.0f"%(gpu_time, sum(gpu[~np.isnan(gpu)]))
        #print max(gpu), min(gpu)

    if 0 and Ncpu > 0: # Hack to compare ctypes vs. opencl
        dllmodel = load_dll(gpuname, dtype='double')
        model = BumpsModel(data, dllmodel, **gpupars)
        toc = tic()
        for i in range(Ncpu):
            model.update()
            cpu = model.theory()
        cpu_time = toc()*1000./Ncpu
        print "dll t=%.1f ms"%cpu_time

    elif 0: # Hack to check new vs old for GpuCylinder
        from Models.code_cylinder_f import GpuCylinder as oldgpu
        from sasmodel import SasModel
        oldmodel = SasModel(data, oldgpu, dtype='single', **cpupars)
        toc = tic()
        for i in range(Ngpu):
            oldmodel.update()
            cpu = oldmodel.theory()
        cpu_time = toc()*1000./Ngpu
        print "old t=%.1f ms"%cpu_time

    elif Ncpu > 0:
        cpumodel = sasview_model(cpuname, **cpupars)
        toc = tic()
        if is2D:
            for i in range(Ncpu):
                cpu = cpumodel.evalDistribution([data.qx_data, data.qy_data])
        else:
            for i in range(Ncpu):
                cpu = cpumodel.evalDistribution(data.x)
        cpu_time = toc()*1000./Ncpu
        print "sasview t=%.1f ms, intensity=%.0f"%(cpu_time, sum(cpu[model.index]))

    if Ngpu > 0 and Ncpu > 0:
        print "gpu/cpu", max(abs(gpu/cpu)), "%.15g"%max(abs(gpu)), "%.15g"%max(abs(cpu))
        #cpu *= max(gpu/cpu)
        abserr = (gpu - cpu)
        relerr = (gpu - cpu)/cpu
        print "max(|ocl-omp|)", max(abs(abserr[model.index]))
        print "max(|(ocl-omp)/ocl|)", max(abs(relerr[model.index]))
    #return

    import matplotlib.pyplot as plt
    if Ncpu > 0:
        if Ngpu > 0: plt.subplot(131)
        plot_data(data, cpu, scale='log')
        plt.title("omp t=%.1f ms"%cpu_time)
    if Ngpu > 0:
        if Ncpu > 0: plt.subplot(132)
        plot_data(data, gpu, scale='log')
        plt.title("ocl t=%.1f ms"%gpu_time)
    if Ncpu > 0 and Ngpu > 0:
        plt.subplot(133)
        plot_data(data, 1e8*relerr, scale='linear')
        plt.title("max rel err = %.3g"%max(abs(relerr)))
        if is2D: plt.colorbar()
    plt.show()

def rename(pars, **names):
    newpars = pars.copy()
    for new,old in names.items():
        for variant in ("", "_pd", "_pd_n", "_pd_nsigma"):
            if old+variant in newpars:
                newpars[new+variant] = pars[old+variant]
                del newpars[old+variant]
    return newpars

def rescale_sld(pars, names):
    newpars = pars.copy()
    for p in names:
        newpars[p] *= 1e6
    return newpars

# ===========================================================================
#
MODELS = {}
def model(name):
    def gather_function(fn):
        MODELS[name] = fn
        return fn
    return gather_function

USAGE="""
usage: compare model [Nopencl] [Nsasview]

Compare the speed and value for a model between the SasView original and the
OpenCL rewrite.

* Nopencl is the number of times to run the OpenCL model (default=5)

* Nsasview is the number of times to run the Sasview model (default=1)

* model is the name of the model to compare:

    %s
"""

def main():
    if len(sys.argv) == 1:
        models = "\n    ".join("%-7s: %s"%(k,v.__name__.replace('_',' '))
                               for k,v in sorted(MODELS.items()))
        print(USAGE%models)
        sys.exit(1)

    cpuname, cpupars, gpuname, gpupars = MODELS[sys.argv[1]]()
    Nopencl = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    Nsasview = int(sys.argv[3]) if len(sys.argv) > 3 else 1

    compare(Nsasview, cpuname, cpupars, Nopencl, gpuname, gpupars)

@model('cyl')
def cylinder():
    cpupars = dict(
        scale=.003, background=.1,
        sldCyl=.291e-6, sldSolv=5.77e-6,
        radius=264.1, length=66.96,
        cyl_theta=85, cyl_phi=0,
        radius_pd=0.1, radius_pd_n=10, radius_pd_nsigma=3,
        length_pd=0.1,length_pd_n=1, length_pd_nsigma=3,
        cyl_theta_pd=45, cyl_theta_pd_n=50, cyl_theta_pd_nsigma=3,
        cyl_phi_pd=0.1, cyl_phi_pd_n=5, cyl_phi_pd_nsigma=3,
        )
    cpuname = 'Cylinder'

    gpupars = rename(cpupars, theta='cyl_theta', phi='cyl_phi', sld='sldCyl', solvent_sld='sldSolv')
    gpupars = rescale_sld(gpupars, ['sld', 'solvent_sld'])
    gpuname = 'cylinder'
    return cpuname, cpupars, gpuname, gpupars


@model('ell')
def ellipse():
    pars = dict(
        scale=.027, background=4.9,
        sldEll=.297e-6, sldSolv=5.773e-6,
        radius_a=60, radius_b=180,
        axis_theta=0, axis_phi=90,
        radius_a_pd=0.1, radius_a_pd_n=10, radius_a_pd_nsigma=3,
        radius_b_pd=0.1, radius_b_pd_n=10, radius_b_pd_nsigma=3,
        axis_theta_pd=0.1, axis_theta_pd_n=6, axis_theta_pd_nsigma=3,
        axis_phi_pd=0.1, axis_phi_pd_n=6, axis_phi_pd_nsigma=3,
        )

    from Models.code_ellipse import GpuEllipse as gpumodel
    model = sasview_model('Ellipsoid', **pars)

    pars = rename(pars, theta='axis_theta', phi='axis_phi', sld='sldEll', solvent_sld='sldSolv')
    pars = rescale_sld(pars, ['sld', 'solvent_sld'])
    return model, gpumodel, pars


@model('cscyl')
def core_shell_cylinder(N=1):
    pars = dict(
        scale= 1.77881e-06, background=223.827,
        core_sld=1e-6, shell_sld=.291e-6, solvent_sld=7.105e-6,
        radius=325, thickness=25, length=34.2709,
        axis_theta=90, axis_phi=0,
        radius_pd=0.1, radius_pd_n=10, radius_pd_nsigma=3,
        length_pd=0.1, length_pd_n=10, length_pd_nsigma=3,
        thickness_pd=0.1, thickness_pd_n=5, thickness_pd_nsigma=3,
        axis_theta_pd=15.8, axis_theta_pd_n=20, axis_theta_pd_nsigma=5,
        axis_phi_pd=0.0008748, axis_phi_pd_n=5, axis_phi_pd_nsigma=3,
        )

    model = sasview_model('CoreShellCylinder', **pars)
    from Models.code_coreshellcyl_f import GpuCoreShellCylinder as gpumodel

    pars = rename(pars, theta='axis_theta', phi='axis_phi')
    pars = rescale_sld(pars, ['core_sld', 'shell_sld', 'solvent_sld'])
    return model, gpumodel, pars


@model('ell3')
def triaxial_ellipse(N=1):
    pars = dict(
        scale=0.08, background=5,
        sldEll=7.105e-6, sldSolv=.291e-6,
        axis_theta=0, axis_phi=0, axis_psi=0,
        semi_axisA=15, semi_axisB=20, semi_axisC=500,
        axis_theta_pd=20, axis_theta_pd_n=10, axis_theta_pd_nsigma=3,
        axis_phi_pd=.1, axis_phi_pd_n=10, axis_phi_pd_nsigma=3,
        axis_psi_pd=30, axis_psi_pd_n=5, axis_psi_pd_nsigma=3,
        semi_axisA_pd=.1, semi_axisA_pd_n=5, semi_axisA_pd_nsigma=3,
        semi_axisB_pd=.1, semi_axisB_pd_n=5, semi_axisB_pd_nsigma=3,
        semi_axisC_pd=.1, semi_axisC_pd_n=5, semi_axisC_pd_nsigma=3,
        )

    model = sasview_model('TriaxialEllipsoid', **pars)
    from Models.code_triaxialellipse import GpuTriEllipse as gpumodel

    pars = rename(pars,
                  theta='axis_theta', phi='axis_phi', psi='axis_psi',
                  sld='sldEll', solvent_sld='sldSolv',
                  radius_a='semi_axisA', radius_b='semi_axisB',
                  radius_c='semi_axisC',
                  )
    pars = rescale_sld(pars, ['sld', 'solvent_sld'])
    return model, gpumodel, pars

@model('lam')
def lamellar(N=1):
    pars = dict(
        scale=0.08, background=0.003,
        sld_bi=5.38e-6,sld_sol=7.105e-6,
        bi_thick=19.2946,
        bi_thick_pd= 0.37765, bi_thick_pd_n=40, bi_thick_pd_nsigma=3,
        )

    model = sasview_model('Lamellar', **pars)
    from Models.code_lamellar import GpuLamellar as gpumodel

    pars = rename(pars, sld='sld_bi', solvent_sld='sld_sol', thickness='bi_thick')
    pars = rescale_sld(pars, ['sld', 'solvent_sld'])
    return model, gpumodel, pars

@model('capcyl')
def capped_cylinder(N=1):
    pars = dict(
        scale=.08, background=0,
        sld_capcyl=1e-6, sld_solv=6.3e-6,
        rad_cyl=20, rad_cap=40, len_cyl=400,
        theta=0, phi=0,
        rad_cyl_pd=.1, rad_cyl_pd_n=10, rad_cyl_pd_nsigma=3,
        rad_cap_pd=.1, rad_cap_pd_n=10, rad_cap_pd_nsigma=3,
        len_cyl_pd=.1, len_cyl_pd_n=3, len_cyl_pd_nsigma=3,
        theta_pd=.1, theta_pd_n=3, theta_pd_nsigma=3,
        phi_pd=.1, phi_pd_n=3, phi_pd_nsigma=3,
        )


    model = sasview_model('CappedCylinder', **pars)
    from Models.code_capcyl import GpuCapCylinder as gpumodel

    pars = rename(pars,
                  sld='sld_capcyl', solvent_sld='sld_solv',
                  length='len_cyl', radius='rad_cyl',
                  cap_radius='rad_cap')
    pars = rescale_sld(pars, ['sld', 'solvent_sld'])
    return model, gpumodel, pars


if __name__ == "__main__":
    main()
