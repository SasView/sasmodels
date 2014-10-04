#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import math

import numpy as np

from sasmodels.bumps_model import BumpsModel, plot_data, tic
from sasmodels import gpu, dll
from sasmodels.convert import revert_model


def sasview_model(modelname, **pars):
    """
    Load a sasview model given the model name.
    """
    # convert model parameters from sasmodel form to sasview form
    #print "old",sorted(pars.items())
    modelname, pars = revert_model(modelname, pars)
    #print "new",sorted(pars.items())
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

def load_ctypes(modelname, dtype='single'):
    sasmodels = __import__('sasmodels.models.'+modelname)
    module = getattr(sasmodels.models, modelname, None)
    kernel = dll.load_model(module, dtype=dtype)
    return kernel

def randomize(p, v):
    """
    Randomizing parameter.

    Guess the parameter type from name.
    """
    if any(p.endswith(s) for s in ('_pd_n','_pd_nsigma','_pd_type')):
        return v
    elif any(s in p for s in ('theta','phi','psi')):
        # orientation in [-180,180], orientation pd in [0,45]
        if p.endswith('_pd'):
            return 45*np.random.rand()
        else:
            return 360*np.random.rand() - 180
    elif 'sld' in p:
        # sld in in [-0.5,10]
        return 10.5*np.random.rand() - 0.5
    elif p.endswith('_pd'):
        # length pd in [0,1]
        return np.random.rand()
    else:
        # length, scale, background in [0,200]
        return 200*np.random.rand()

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

def compare(name, pars, Ncpu, Ngpu, opts, set_pars):

    # Sort out data
    qmax = 1.0 if '-highq' in opts else (0.2 if '-midq' in opts else 0.05)
    if "-1d" in opts:
        from sasmodels.bumps_model import empty_data1D
        qmax = math.log10(qmax)
        data = empty_data1D(np.logspace(qmax-3, qmax, 128))
        index = slice(None, None)
    else:
        from sasmodels.bumps_model import empty_data2D, set_beam_stop
        data = empty_data2D(np.linspace(-qmax, qmax, 128))
        set_beam_stop(data, 0.004)
        index = ~data.mask
    is2D = hasattr(data, 'qx_data')

    # modelling accuracy is determined by dtype and cutoff
    dtype = 'double' if '-double' in opts else 'single'
    cutoff_opts = [s for s in opts if s.startswith('-cutoff')]
    cutoff = float(cutoff_opts[0].split('=')[1]) if cutoff_opts else 1e-5

    # randomize parameters
    random_opts = [s for s in opts if s.startswith('-random')]
    if random_opts:
        if '=' in random_opts[0]:
            seed = int(random_opts[0].split('=')[1])
        else:
            seed = int(np.random.rand()*10000)
        np.random.seed(seed)
        print "Randomize using -random=%i"%seed
        # Note: the sort guarantees order of calls to random number generator
        pars = dict((p,randomize(p,v)) for p,v in sorted(pars.items()))
        # The capped cylinder model has a constraint on its parameters
        if name == 'capped_cylinder' and pars['cap_radius'] < pars['radius']:
            pars['radius'],pars['cap_radius'] = pars['cap_radius'],pars['radius']
    pars.update(set_pars)

    # parameter selection
    if '-mono' in opts:
        suppress_pd(pars)
    if '-pars' in opts:
        print "pars",parlist(pars)

    # OpenCl calculation
    if Ngpu > 0:
        try:
            model = load_opencl(name, dtype=dtype)
        except Exception,exc:
            print exc
            print "... trying again with single precision"
            model = load_opencl(name, dtype='single')
        problem = BumpsModel(data, model, cutoff=cutoff, **pars)
        toc = tic()
        for _ in range(Ngpu):
            #pars['scale'] = np.random.rand()
            problem.update()
            gpu = problem.theory()
        gpu_time = toc()*1000./Ngpu
        print "opencl t=%.1f ms, intensity=%.0f"%(gpu_time, sum(gpu[index]))
        #print max(gpu), min(gpu)

    # ctypes/sasview calculation
    if Ncpu > 0 and "-ctypes" in opts:
        model = load_ctypes(name, dtype=dtype)
        problem = BumpsModel(data, model, cutoff=cutoff, **pars)
        toc = tic()
        for _ in range(Ncpu):
            problem.update()
            cpu = problem.theory()
        cpu_time = toc()*1000./Ncpu
        comp = "ctypes"
        print "ctypes t=%.1f ms, intensity=%.0f"%(cpu_time, sum(cpu[index]))
    elif Ncpu > 0:
        model = sasview_model(name, **pars)
        toc = tic()
        for _ in range(Ncpu):
            if is2D:
                cpu = model.evalDistribution([data.qx_data, data.qy_data])
            else:
                cpu = model.evalDistribution(data.x)
        cpu_time = toc()*1000./Ncpu
        comp = "sasview"
        print "sasview t=%.1f ms, intensity=%.0f"%(cpu_time, sum(cpu[index]))

    # Compare, but only if computing both forms
    if Ngpu > 0 and Ncpu > 0:
        #print "speedup %.2g"%(cpu_time/gpu_time)
        #print "max |gpu/cpu|", max(abs(gpu/cpu)), "%.15g"%max(abs(gpu)), "%.15g"%max(abs(cpu))
        #cpu *= max(gpu/cpu)
        resid, relerr = np.zeros_like(gpu), np.zeros_like(gpu)
        resid[index] = (gpu - cpu)[index]
        relerr[index] = resid[index]/cpu[index]
        print "max(|ocl-%s|)"%comp, max(abs(resid[index]))
        print "max(|(ocl-%s)/ocl|)"%comp, max(abs(relerr[index]))

    # Plot if requested
    if '-noplot' in opts: return
    import matplotlib.pyplot as plt
    if Ncpu > 0:
        if Ngpu > 0: plt.subplot(131)
        plot_data(data, cpu, scale='log')
        plt.title("%s t=%.1f ms"%(comp,cpu_time))
    if Ngpu > 0:
        if Ncpu > 0: plt.subplot(132)
        plot_data(data, gpu, scale='log')
        plt.title("opencl t=%.1f ms"%gpu_time)
    if Ncpu > 0 and Ngpu > 0:
        plt.subplot(133)
        err = resid if '-abs' in opts else relerr
        errstr = "abs err" if '-abs' in opts else "rel err"
        #err,errstr = gpu/cpu,"ratio"
        plot_data(data, err, scale='linear')
        plt.title("max %s = %.3g"%(errstr, max(abs(err[index]))))
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
    -single*/-double uses double precision for comparison
    -lowq*/-midq/-highq use q values up to 0.05, 0.2 or 1.0
    -1d/-2d* uses 1d or 2d random data
    -preset*/-random[=seed] preset or random parameters
    -mono/-poly* force monodisperse/polydisperse
    -ctypes/-sasview* whether cpu is tested using sasview or ctypes
    -cutoff=1e-5*/value cutoff for including a point in polydispersity
    -pars/-nopars* prints the parameter set or not
    -abs/-rel* plot relative or absolute error

Key=value pairs allow you to set specific values to any of the model
parameters.

Available models:

    %s
"""

VALID_OPTIONS = [
    'plot','noplot',
    'single','double',
    'lowq','midq','highq',
    '2d','1d',
    'preset','random',
    'poly','mono',
    'sasview','ctypes',
    'nopars','pars',
    'rel','abs',
    ]

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

    valid_opts = set(VALID_OPTIONS)
    invalid = [o[1:] for o in opts
               if o[1:] not in valid_opts
                  and not o.startswith('-cutoff=')
                  and not o.startswith('-random=')]
    if invalid:
        print "Invalid options: %s"%(", ".join(invalid))
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
    set_pars = {}
    for arg in args[3:]:
        k,v = arg.split('=')
        if k not in pars:
            # extract base name without distribution
            s = set(p.split('_pd')[0] for p in pars)
            print "%r invalid; parameters are: %s"%(k,", ".join(sorted(s)))
            sys.exit(1)
        set_pars[k] = float(v) if not v.endswith('type') else v

    compare(name, pars, Nsasview, Nopencl, opts, set_pars)

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
        sld=6, solvent_sld=1,
        #radius=5, length=20,
        radius=260, length=290,
        theta=30, phi=0,
        radius_pd=.2, radius_pd_n=1,
        length_pd=.2,length_pd_n=1,
        theta_pd=15, theta_pd_n=45,
        phi_pd=15, phi_pd_n=1,
        )
    return 'cylinder', pars

@model('capcyl')
def capped_cylinder():
    pars = dict(
        scale=1, background=0,
        sld=6, solvent_sld=1,
        radius=260, cap_radius=290, length=290,
        theta=30, phi=15,
        radius_pd=.2, radius_pd_n=1,
        cap_radius_pd=.2, cap_radius_pd_n=1,
        length_pd=.2, length_pd_n=1,
        theta_pd=15, theta_pd_n=45,
        phi_pd=15, phi_pd_n=1,
        )
    return 'capped_cylinder', pars


@model('cscyl')
def core_shell_cylinder():
    pars = dict(
        scale=1, background=0,
        core_sld=6, shell_sld=8, solvent_sld=1,
        radius=45, thickness=25, length=340,
        theta=30, phi=15,
        radius_pd=.2, radius_pd_n=1,
        length_pd=.2, length_pd_n=1,
        thickness_pd=.2, thickness_pd_n=1,
        theta_pd=15, theta_pd_n=45,
        phi_pd=15, phi_pd_n=1,
        )
    return 'core_shell_cylinder', pars


@model('ell')
def ellipsoid():
    pars = dict(
        scale=1, background=0,
        sld=6, solvent_sld=1,
        rpolar=50, requatorial=30,
        theta=30, phi=15,
        rpolar_pd=.2, rpolar_pd_n=1,
        requatorial_pd=.2, requatorial_pd_n=1,
        theta_pd=15, theta_pd_n=45,
        phi_pd=15, phi_pd_n=1,
        )
    return 'ellipsoid', pars


@model('ell3')
def triaxial_ellipsoid():
    pars = dict(
        scale=1, background=0,
        sld=6, solvent_sld=1,
        theta=30, phi=15, psi=5,
        req_minor=25, req_major=36, rpolar=50,
        req_minor_pd=0, req_minor_pd_n=1,
        req_major_pd=0, req_major_pd_n=1,
        rpolar_pd=.2, rpolar_pd_n=1,
        theta_pd=15, theta_pd_n=45,
        phi_pd=15, phi_pd_n=1,
        psi_pd=15, psi_pd_n=1,
        )
    return 'triaxial_ellipsoid', pars

@model('sph')
def sphere():
    pars = dict(
        scale=1, background=0,
        sld=6, solvent_sld=1,
        radius=120,
        radius_pd=.2, radius_pd_n=45,
        )
    return 'sphere', pars

@model('lam')
def lamellar():
    pars = dict(
        scale=1, background=0,
        sld=6, solvent_sld=1,
        thickness=40,
        thickness_pd= 0.2, thickness_pd_n=40,
        )
    return 'lamellar', pars

if __name__ == "__main__":
    main()
