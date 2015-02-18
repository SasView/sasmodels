#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import math
from os.path import basename, dirname, join as joinpath
import glob

import numpy as np

from sasmodels.bumps_model import BumpsModel, plot_data, tic
try: from sasmodels import kernelcl
except: from sasmodels import kerneldll as kernelcl
from sasmodels import kerneldll
from sasmodels.convert import revert_model

# List of available models
ROOT = dirname(__file__)
MODELS = [basename(f)[:-3]
          for f in sorted(glob.glob(joinpath(ROOT,"sasmodels","models","[a-zA-Z]*.py")))]


def sasview_model(modelname, **pars):
    """
    Load a sasview model given the model name.
    """
    # convert model parameters from sasmodel form to sasview form
    #print "old",sorted(pars.items())
    modelname, pars = revert_model(modelname, pars)
    #print "new",sorted(pars.items())
    sas = __import__('sas.models.'+modelname)
    ModelClass = getattr(getattr(sas.models,modelname,None),modelname,None)
    if ModelClass is None:
        raise ValueError("could not find model %r in sas.models"%modelname)
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
    kernel = kernelcl.load_model(module, dtype=dtype)
    return kernel

def load_ctypes(modelname, dtype='single'):
    sasmodels = __import__('sasmodels.models.'+modelname)
    module = getattr(sasmodels.models, modelname, None)
    kernel = kerneldll.load_model(module, dtype=dtype)
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

def randomize_model(name, pars, seed=None):
    if seed is None:
        seed = np.random.randint(1e9)
    np.random.seed(seed)
    # Note: the sort guarantees order of calls to random number generator
    pars = dict((p,randomize(p,v)) for p,v in sorted(pars.items()))
    # The capped cylinder model has a constraint on its parameters
    if name == 'capped_cylinder' and pars['cap_radius'] < pars['radius']:
        pars['radius'],pars['cap_radius'] = pars['cap_radius'],pars['radius']
    return pars, seed

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

def eval_sasview(name, pars, data, Nevals=1):
    model = sasview_model(name, **pars)
    toc = tic()
    for _ in range(Nevals):
        if hasattr(data, 'qx_data'):
            value = model.evalDistribution([data.qx_data, data.qy_data])
        else:
            value = model.evalDistribution(data.x)
    average_time = toc()*1000./Nevals
    return value, average_time

def eval_opencl(name, pars, data, dtype='single', Nevals=1, cutoff=0):
    try:
        model = load_opencl(name, dtype=dtype)
    except Exception,exc:
        print exc
        print "... trying again with single precision"
        model = load_opencl(name, dtype='single')
    problem = BumpsModel(data, model, cutoff=cutoff, **pars)
    toc = tic()
    for _ in range(Nevals):
        #pars['scale'] = np.random.rand()
        problem.update()
        value = problem.theory()
    average_time = toc()*1000./Nevals
    return value, average_time

def eval_ctypes(name, pars, data, dtype='double', Nevals=1, cutoff=0):
    model = load_ctypes(name, dtype=dtype)
    problem = BumpsModel(data, model, cutoff=cutoff, **pars)
    toc = tic()
    for _ in range(Nevals):
        problem.update()
        value = problem.theory()
    average_time = toc()*1000./Nevals
    return value, average_time

def make_data(qmax, is2D, Nq=128):
    if is2D:
        from sasmodels.bumps_model import empty_data2D, set_beam_stop
        data = empty_data2D(np.linspace(-qmax, qmax, Nq))
        set_beam_stop(data, 0.004)
        index = ~data.mask
    else:
        from sasmodels.bumps_model import empty_data1D
        qmax = math.log10(qmax)
        data = empty_data1D(np.logspace(qmax-3, qmax, Nq))
        index = slice(None, None)
    return data, index

def compare(name, pars, Ncpu, Nocl, opts, set_pars):
    opt_values = dict(split
                      for s in opts for split in ((s.split('='),))
                      if len(split) == 2)
    # Sort out data
    qmax = 1.0 if '-highq' in opts else (0.2 if '-midq' in opts else 0.05)
    Nq = int(opt_values.get('-Nq', '128'))
    is2D = not "-1d" in opts
    data, index = make_data(qmax, is2D, Nq)


    # modelling accuracy is determined by dtype and cutoff
    dtype = 'double' if '-double' in opts else 'single'
    cutoff = float(opt_values.get('-cutoff','1e-5'))

    # randomize parameters
    if '-random' in opts or '-random' in opt_values:
        seed = int(opt_values['-random']) if '-random' in opt_values else None
        pars, seed = randomize_model(name, pars, seed=seed)
        print "Randomize using -random=%i"%seed
    pars.update(set_pars)

    # parameter selection
    if '-mono' in opts:
        suppress_pd(pars)
    if '-pars' in opts:
        print "pars",parlist(pars)

    # OpenCl calculation
    if Nocl > 0:
        ocl, ocl_time = eval_opencl(name, pars, data, dtype, Nocl)
        print "opencl t=%.1f ms, intensity=%.0f"%(ocl_time, sum(ocl[index]))
        #print max(ocl), min(ocl)

    # ctypes/sasview calculation
    if Ncpu > 0 and "-ctypes" in opts:
        cpu, cpu_time = eval_ctypes(name, pars, data, dtype=dtype, cutoff=cutoff, Nevals=Ncpu)
        comp = "ctypes"
        print "ctypes t=%.1f ms, intensity=%.0f"%(cpu_time, sum(cpu[index]))
    elif Ncpu > 0:
        cpu, cpu_time = eval_sasview(name, pars, data, Ncpu)
        comp = "sasview"
        print "sasview t=%.1f ms, intensity=%.0f"%(cpu_time, sum(cpu[index]))

    # Compare, but only if computing both forms
    if Nocl > 0 and Ncpu > 0:
        #print "speedup %.2g"%(cpu_time/ocl_time)
        #print "max |ocl/cpu|", max(abs(ocl/cpu)), "%.15g"%max(abs(ocl)), "%.15g"%max(abs(cpu))
        #cpu *= max(ocl/cpu)
        resid, relerr = np.zeros_like(ocl), np.zeros_like(ocl)
        resid[index] = (ocl - cpu)[index]
        relerr[index] = resid[index]/cpu[index]
        #bad = (relerr>1e-4)
        #print relerr[bad],cpu[bad],ocl[bad],data.qx_data[bad],data.qy_data[bad]
        print "max(|ocl-%s|)"%comp, max(abs(resid[index]))
        print "max(|(ocl-%s)/%s|)"%(comp,comp), max(abs(relerr[index]))
        p98 = int(len(relerr[index])*0.98)
        print "98%% (|(ocl-%s)/%s|) <"%(comp,comp), np.sort(abs(relerr[index]))[p98]


    # Plot if requested
    if '-noplot' in opts: return
    import matplotlib.pyplot as plt
    if Ncpu > 0:
        if Nocl > 0: plt.subplot(131)
        plot_data(data, cpu, scale='log')
        plt.title("%s t=%.1f ms"%(comp,cpu_time))
    if Nocl > 0:
        if Ncpu > 0: plt.subplot(132)
        plot_data(data, ocl, scale='log')
        plt.title("opencl t=%.1f ms"%ocl_time)
    if Ncpu > 0 and Nocl > 0:
        plt.subplot(133)
        err = resid if '-abs' in opts else relerr
        errstr = "abs err" if '-abs' in opts else "rel err"
        #err,errstr = ocl/cpu,"ratio"
        plot_data(data, err, scale='linear')
        plt.title("max %s = %.3g"%(errstr, max(abs(err[index]))))
    if is2D: plt.colorbar()

    if Ncpu > 0 and Nocl > 0 and '-hist' in opts:
        plt.figure()
        v = relerr[index]
        v[v==0] = 0.5*np.min(np.abs(v[v!=0]))
        plt.hist(np.log10(np.abs(v)), normed=1, bins=50);
        plt.xlabel('log10(err), err = | F(q) single - F(q) double| / | F(q) double |');
        plt.ylabel('P(err)')
        plt.title('Comparison of single and double precision models for %s'%name)

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
    -Nq=128 sets the number of Q points in the data set
    -1d/-2d* computes 1d or 2d data
    -preset*/-random[=seed] preset or random parameters
    -mono/-poly* force monodisperse/polydisperse
    -ctypes/-sasview* whether cpu is tested using sasview or ctypes
    -cutoff=1e-5*/value cutoff for including a point in polydispersity
    -pars/-nopars* prints the parameter set or not
    -abs/-rel* plot relative or absolute error
    -hist/-nohist* plot histogram of relative error

Key=value pairs allow you to set specific values to any of the model
parameters.

Available models:

    %s
"""

NAME_OPTIONS = set([
    'plot','noplot',
    'single','double',
    'lowq','midq','highq',
    '2d','1d',
    'preset','random',
    'poly','mono',
    'sasview','ctypes',
    'nopars','pars',
    'rel','abs',
    'hist','nohist',
    ])
VALUE_OPTIONS = [
    # Note: random is both a name option and a value option
    'cutoff', 'random', 'Nq',
    ]

def get_demo_pars(name):
    import sasmodels.models
    __import__('sasmodels.models.'+name)
    model = getattr(sasmodels.models, name)
    pars = getattr(model, 'demo', None)
    if pars is None: pars = dict((p[0],p[2]) for p in model.parameters)
    return pars

def main():
    opts = [arg for arg in sys.argv[1:] if arg.startswith('-')]
    args = [arg for arg in sys.argv[1:] if not arg.startswith('-')]
    models = "\n    ".join("%-15s"%v for v in MODELS)
    if len(args) == 0:
        print(USAGE%models)
        sys.exit(1)
    if args[0] not in MODELS:
        print "Model %r not available. Use one of:\n    %s"%(args[0],models)
        sys.exit(1)

    invalid = [o[1:] for o in opts
               if o[1:] not in NAME_OPTIONS
                  and not any(o.startswith('-%s='%t) for t in VALUE_OPTIONS)]
    if invalid:
        print "Invalid options: %s"%(", ".join(invalid))
        sys.exit(1)

    # Get demo parameters from model definition, or use default parameters
    # if model does not define demo parameters
    name = args[0]
    pars = get_demo_pars(name)

    Nopencl = int(args[1]) if len(args) > 1 else 5
    Nsasview = int(args[2]) if len(args) > 2 else 1

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

if __name__ == "__main__":
    main()
