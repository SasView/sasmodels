#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import numpy as np

from sasmodels.bumps_model import BumpsModel, plot_data, tic
from sasmodels import gpu, dll
from sasmodels.convert import revert_model


def sasview_model(modelname, **pars):
    """
    Load a sasview model given the model name.
    """
    modelname = modelname
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
    kernel = gpu.load_model(module, dtype=dtype)
    return kernel

def load_dll(modelname, dtype='single'):
    sasmodels = __import__('sasmodels.models.'+modelname)
    module = getattr(sasmodels.models, modelname, None)
    kernel = dll.load_model(module, dtype=dtype)
    return kernel

def randomize(p, v):
    """
    Randomizing pars using sasview names.
    Guessing at angles and slds.
    """
    if any(p.endswith(s) for s in ('_pd_n','_pd_nsigmas','_pd_type')):
        return v
    if any(s in p for s in ('theta','phi','psi')):
        return np.random.randint(0,90)
    elif any(s in p for s in ('sld','Sld','SLD')):
        return np.random.rand()*1e-5
    else:
        return np.random.rand()

def parlist(pars):
    return "\n".join("%s: %s"%(p,v) for p,v in sorted(pars.items()))

def suppress_pd(pars):
    """
    Suppress theta_pd for now until the normalization is resolved.

    May also suppress complete polydispersity of the model to test
    models more quickly.
    """
    for p in pars:
        if p.endswith("_pd"): pars[p] = 0

def compare(gpuname, gpupars, Ncpu, Ngpu, opts):
    #from sasmodels.core import load_data
    #data = load_data('December/DEC07098.DAT')
    qmax = 1.0 if '-highq' in opts else (0.2 if '-midq' in opts else 0.05)
    if "-1d" in opts:
        from sasmodels.bumps_model import empty_data1D
        qmax = np.log10(qmax)
        data = empty_data1D(np.logspace(qmax-3, qmax, 128))
    else:
        from sasmodels.bumps_model import empty_data2D, set_beam_stop
        data = empty_data2D(np.linspace(-qmax, qmax, 128))
        set_beam_stop(data, 0.004)
    is2D = hasattr(data, 'qx_data')
    dtype = 'double' if '-double' in opts else 'single'
    cutoff_opts = [s for s in opts if s.startswith('-cutoff')]
    cutoff = float(cutoff_opts[0].split('=')[1]) if cutoff_opts else 1e-5


    if '-random' in opts:
        gpupars = dict((p,randomize(p,v)) for p,v in gpupars.items())
    if '-pars' in opts: print "pars",parlist(gpupars)
    if '-mono' in opts: suppress_pd(gpupars)

    cpuname, cpupars = revert_model(gpuname, gpupars)

    try:
        gpumodel = load_opencl(gpuname, dtype=dtype)
    except Exception,exc:
        print exc
        print "... trying again with single precision"
        gpumodel = load_opencl(gpuname, dtype='single')
    model = BumpsModel(data, gpumodel, cutoff=cutoff, **gpupars)
    if Ngpu > 0:
        toc = tic()
        for i in range(Ngpu):
            #pars['scale'] = np.random.rand()
            model.update()
            gpu = model.theory()
        gpu_time = toc()*1000./Ngpu
        print "ocl t=%.1f ms, intensity=%f"%(gpu_time, sum(gpu[~np.isnan(gpu)]))
        #print max(gpu), min(gpu)

    comp = None
    if Ncpu > 0 and "-dll" in opts:
        dllmodel = load_dll(gpuname, dtype=dtype)
        model = BumpsModel(data, dllmodel, cutoff=cutoff, **gpupars)
        toc = tic()
        for i in range(Ncpu):
            model.update()
            cpu = model.theory()
        cpu_time = toc()*1000./Ncpu
        comp = "dll"

    elif 0: # Hack to check new vs old for GpuCylinder
        from Models.code_cylinder_f import GpuCylinder as oldgpu
        from sasmodel import SasModel
        oldmodel = SasModel(data, oldgpu, dtype=dtype, **cpupars)
        toc = tic()
        for i in range(Ngpu):
            oldmodel.update()
            cpu = oldmodel.theory()
        cpu_time = toc()*1000./Ngpu
        comp = "old"

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
        comp = 'sasview'
        #print "sasview t=%.1f ms, intensity=%.0f"%(cpu_time, sum(cpu[model.index]))

    if comp:
        print "%s t=%.1f ms, intensity=%f"%(comp, cpu_time, sum(cpu[model.index]))
    if Ngpu > 0 and Ncpu > 0:
        #print "speedup %.2g"%(cpu_time/gpu_time)
        #print "max |gpu/cpu|", max(abs(gpu/cpu)), "%.15g"%max(abs(gpu)), "%.15g"%max(abs(cpu))
        #cpu *= max(gpu/cpu)
        resid, relerr = np.zeros_like(gpu), np.zeros_like(gpu)
        resid[model.index] = (gpu - cpu)[model.index]
        relerr[model.index] = resid[model.index]/cpu[model.index]
        print "max(|ocl-%s|)"%comp, max(abs(resid[model.index]))
        print "max(|(ocl-%s)/ocl|)"%comp, max(abs(relerr[model.index]))

    if '-noplot' in opts: return

    import matplotlib.pyplot as plt
    if Ncpu > 0:
        if Ngpu > 0: plt.subplot(131)
        plot_data(data, cpu, scale='log')
        plt.title("%s t=%.1f ms"%(comp,cpu_time))
    if Ngpu > 0:
        if Ncpu > 0: plt.subplot(132)
        plot_data(data, gpu, scale='log')
        plt.title("ocl t=%.1f ms"%gpu_time)
    if Ncpu > 0 and Ngpu > 0:
        plt.subplot(133)
        err = resid if '-abs' in opts else relerr
        plot_data(data, err, scale='linear')
        plt.title("max rel err = %.3g"%max(abs(err[model.index])))
        if is2D: plt.colorbar()
    plt.show()

# ===========================================================================
#
USAGE="""
usage: compare.py model [Nopencl] [Nsasview] [options...] [key=val]

Compare the speed and value for a model between the SasView original and the
OpenCL rewrite.

model is the name of the model to compare (see below).
Nopencl is the number of times to run the OpenCL model (default=5)
Nsasview is the number of times to run the Sasview model (default=1)

Options (* for default):

    -plot*/-noplot plots or suppress the plot of the model
    -single/-double uses double precision for comparison
    -lowq/-midq/-highq use q values up to 0.05, 0.2 or 1.0
    -2d/-1d uses 1d or 2d random data
    -preset/-random randomizes the parameters
    -poly/-mono force monodisperse/polydisperse
    -sasview/-dll whether cpu is tested using sasview or dll
    -cutoff=1e-5/value cutoff for including a point in polydispersity
    -nopars/-pars prints the parameter set or not
    -rel/-abs plot relative or absolute error

Key=value pairs allow you to set specific values to any of the model
parameters.

Available models:

    %s
"""

def main():
    opts = [arg for arg in sys.argv[1:] if arg.startswith('-')]
    args = [arg for arg in sys.argv[1:] if not arg.startswith('-')]
    models = "\n    ".join("%-7s: %s"%(k,v.__name__.replace('_',' '))
                           for k,v in sorted(MODELS.items()))
    if len(args) == 0:
        print(USAGE%models)
        sys.exit(1)
    if args[0] not in MODELS:
        print "Model %r not available. Use one of:\n    %s"%(args[0],models)
        sys.exit(1)

    name, pars = MODELS[args[0]]()
    Nopencl = int(args[1]) if len(args) > 1 else 5
    Nsasview = int(args[2]) if len(args) > 3 else 1

    # Fill in default polydispersity parameters
    pds = set(p.split('_pd')[0] for p in pars if p.endswith('_pd'))
    for p in pds:
        if p+"_pd_nsigma" not in pars: pars[p+"_pd_nsigma"] = 3
        if p+"_pd_type" not in pars: pars[p+"_pd_type"] = "gaussian"

    # Fill in parameters given on the command line
    for arg in args[3:]:
        k,v = arg.split('=')
        if k not in pars:
            # extract base name without distribution
            # style may be either a.d or a_pd_d
            s = set((p.split('_pd')[0]).split('.')[0] for p in pars)
            print "%r invalid; parameters are: %s"%(k,", ".join(sorted(s)))
            sys.exit(1)
        pars[k] = float(v) if not v.endswith('type') else v

    compare(name, pars, Nsasview, Nopencl, opts)

# ===========================================================================
#

MODELS = {}
def model(name):
    def gather_function(fn):
        MODELS[name] = fn
        return fn
    return gather_function


@model('cyl')
def cylinder():
    pars = dict(
        scale=1, background=0,
        sld=6e-6, solvent_sld=1e-6,
        #radius=5, length=20,
        radius=260, length=290,
        theta=30, phi=15,
        radius_pd=.2, radius_pd_n=1,
        length_pd=.2,length_pd_n=1,
        theta_pd=15, theta_pd_n=25,
        phi_pd=15, phi_pd_n=1,
        )
    return 'cylinder', pars

@model('capcyl')
def capped_cylinder():
    pars = dict(
        scale=1, background=0,
        sld=6e-6, solvent_sld=1e-6,
        radius=260, cap_radius=80000, length=290,
        theta=30, phi=15,
        radius_pd=.2, radius_pd_n=1,
        cap_radius_pd=.2, cap_radius_pd_n=1,
        length_pd=.2, length_pd_n=1,
        theta_pd=15, theta_pd_n=25,
        phi_pd=15, phi_pd_n=1,
        )
    return 'capped_cylinder', pars


@model('cscyl')
def core_shell_cylinder():
    pars = dict(
        scale=1, background=0,
        core_sld=6e-6, shell_sld=8e-6, solvent_sld=1e-6,
        radius=325, thickness=25, length=34.2709,
        theta=90, phi=0,
        radius_pd=0.1, radius_pd_n=10,
        length_pd=0.1, length_pd_n=5,
        thickness_pd=0.1, thickness_pd_n=5,
        theta_pd=15.8, theta_pd_n=5,
        phi_pd=0.0008748, phi_pd_n=5,
        )
    return 'core_shell_cylinder', pars


@model('ell')
def ellipsoid():
    pars = dict(
        scale=1, background=0,
        sld=6e-6, solvent_sld=1e-6,
        rpolar=50, requatorial=30,
        theta=0, phi=0,
        rpolar_pd=0.3, rpolar_pd_n=10,
        requatorial_pd=0, requatorial_pd_n=10,
        theta_pd=0, theta_pd_n=45,
        phi_pd=0, phi_pd_n=45,
        )
    return 'ellipsoid', pars


@model('ell3')
def triaxial_ellipsoid():
    pars = dict(
        scale=1, background=0,
        sld=6e-6, solvent_sld=1e-6,
        theta=0, phi=0, psi=0,
        req_minor=25, req_major=36, rpolar=50,
        theta_pd=0, theta_pd_n=5,
        phi_pd=0, phi_pd_n=5,
        psi_pd=0, psi_pd_n=5,
        req_minor_pd=0, req_minor_pd_n=5,
        req_major_pd=0, req_major_pd_n=5,
        rpolar_pd=.3, rpolar_pd_n=25,
        )
    return 'triaxial_ellipsoid', pars

@model('sph')
def sphere():
    pars = dict(
        scale=1, background=0,
        sld=6e-6, solvent_sld=1e-6,
        radius=120,
        radius_pd=.3, radius_pd_n=5,
        )
    return 'sphere', pars

@model('lam')
def lamellar():
    pars = dict(
        scale=1, background=0,
        sld=6e-6, solvent_sld=1e-6,
        thickness=40,
        thickness_pd= 0.3, thickness_pd_n=40,
        )
    return 'lamellar', pars

if __name__ == "__main__":
    main()
