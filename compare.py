#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import math
from os.path import basename, dirname, join as joinpath
import glob

import numpy as np

from sasmodels.bumps_model import Model, Experiment, plot_theory, tic
from sasmodels import core
from sasmodels import kerneldll
from sasmodels.convert import revert_model
kerneldll.ALLOW_SINGLE_PRECISION_DLLS = True

# List of available models
ROOT = dirname(__file__)
MODELS = [basename(f)[:-3]
          for f in sorted(glob.glob(joinpath(ROOT,"sasmodels","models","[a-zA-Z]*.py")))]


def sasview_model(model_definition, **pars):
    """
    Load a sasview model given the model name.
    """
    # convert model parameters from sasmodel form to sasview form
    #print "old",sorted(pars.items())
    modelname, pars = revert_model(model_definition, pars)
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
        # values from 0 to 2*x for all other parameters
        return 2*np.random.rand()*(v if v != 0 else 1)

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
    from sas.models.qsmearing import smear_selection
    model = sasview_model(name, **pars)
    smearer = smear_selection(data, model=model)
    value = None  # silence the linter
    toc = tic()
    for _ in range(max(Nevals, 1)):  # make sure there is at least one eval
        if hasattr(data, 'qx_data'):
            q = np.sqrt(data.qx_data**2 + data.qy_data**2)
            index = ((~data.mask) & (~np.isnan(data.data))
                     & (q >= data.qmin) & (q <= data.qmax))
            if smearer is not None:
                smearer.model = model  # because smear_selection has a bug
                smearer.accuracy = data.accuracy
                smearer.set_index(index)
                value = smearer.get_value()
            else:
                value = model.evalDistribution([data.qx_data[index], data.qy_data[index]])
        else:
            value = model.evalDistribution(data.x)
            if smearer is not None:
                value = smearer(value)
    average_time = toc()*1000./Nevals
    return value, average_time

def eval_opencl(model_definition, pars, data, dtype='single', Nevals=1, cutoff=0.):
    try:
        model = core.load_model(model_definition, dtype=dtype, platform="ocl")
    except Exception,exc:
        print exc
        print "... trying again with single precision"
        model = core.load_model(model_definition, dtype='single', platform="ocl")
    problem = Experiment(data, Model(model, **pars), cutoff=cutoff)
    value = None  # silence the linter
    toc = tic()
    for _ in range(max(Nevals, 1)):  # force at least one eval
        #pars['scale'] = np.random.rand()
        problem.update()
        value = problem.theory()
    average_time = toc()*1000./Nevals
    return value, average_time

def eval_ctypes(model_definition, pars, data, dtype='double', Nevals=1, cutoff=0.):
    model = core.load_model(model_definition, dtype=dtype, platform="dll")
    problem = Experiment(data, Model(model, **pars), cutoff=cutoff)
    value = None  # silence the linter
    toc = tic()
    for _ in range(max(Nevals, 1)):  # force at least one eval
        problem.update()
        value = problem.theory()
    average_time = toc()*1000./Nevals
    return value, average_time

def make_data(qmax, is2D, Nq=128, resolution=0.0, accuracy='Low', view='log'):
    if is2D:
        from sasmodels.bumps_model import empty_data2D, set_beam_stop
        data = empty_data2D(np.linspace(-qmax, qmax, Nq), resolution=resolution)
        data.accuracy = accuracy
        set_beam_stop(data, 0.004)
        index = ~data.mask
    else:
        from sasmodels.bumps_model import empty_data1D
        if view == 'log':
            qmax = math.log10(qmax)
            q = np.logspace(qmax-3, qmax, Nq)
        else:
            q = np.linspace(0.001*qmax, qmax, Nq)
        data = empty_data1D(q, resolution=resolution)
        index = slice(None, None)
    return data, index

def compare(name, pars, Ncpu, Nocl, opts, set_pars):
    view = 'linear' if '-linear' in opts else 'log' if '-log' in opts else 'q4' if '-q4' in opts else 'log'

    opt_values = dict(split
                      for s in opts for split in ((s.split('='),))
                      if len(split) == 2)
    # Sort out data
    qmax = 10.0 if '-exq' in opts else 1.0 if '-highq' in opts else 0.2 if '-midq' in opts else 0.05
    Nq = int(opt_values.get('-Nq', '128'))
    res = float(opt_values.get('-res', '0'))
    accuracy = opt_values.get('-accuracy', 'Low')
    is2D = not "-1d" in opts
    data, index = make_data(qmax, is2D, Nq, res, accuracy, view=view)


    # modelling accuracy is determined by dtype and cutoff
    dtype = 'double' if '-double' in opts else 'single'
    cutoff = float(opt_values.get('-cutoff','1e-5'))

    # randomize parameters
    pars.update(set_pars)
    if '-random' in opts or '-random' in opt_values:
        seed = int(opt_values['-random']) if '-random' in opt_values else None
        pars, seed = randomize_model(name, pars, seed=seed)
        print "Randomize using -random=%i"%seed

    # parameter selection
    if '-mono' in opts:
        suppress_pd(pars)
    if '-pars' in opts:
        print "pars",parlist(pars)

    model_definition = core.load_model_definition(name)
    # OpenCl calculation
    if Nocl > 0:
        ocl, ocl_time = eval_opencl(model_definition, pars, data,
                                    dtype=dtype, cutoff=cutoff, Nevals=Nocl)
        print "opencl t=%.1f ms, intensity=%.0f"%(ocl_time, sum(ocl))
        #print max(ocl), min(ocl)

    # ctypes/sasview calculation
    if Ncpu > 0 and "-ctypes" in opts:
        cpu, cpu_time = eval_ctypes(model_definition, pars, data,
                                    dtype=dtype, cutoff=cutoff, Nevals=Ncpu)
        comp = "ctypes"
        print "ctypes t=%.1f ms, intensity=%.0f"%(cpu_time, sum(cpu))
    elif Ncpu > 0:
        cpu, cpu_time = eval_sasview(model_definition, pars, data, Ncpu)
        comp = "sasview"
        print "sasview t=%.1f ms, intensity=%.0f"%(cpu_time, sum(cpu))

    # Compare, but only if computing both forms
    if Nocl > 0 and Ncpu > 0:
        #print "speedup %.2g"%(cpu_time/ocl_time)
        #print "max |ocl/cpu|", max(abs(ocl/cpu)), "%.15g"%max(abs(ocl)), "%.15g"%max(abs(cpu))
        #cpu *= max(ocl/cpu)
        resid = (ocl - cpu)
        relerr = resid/cpu
        #bad = (relerr>1e-4)
        #print relerr[bad],cpu[bad],ocl[bad],data.qx_data[bad],data.qy_data[bad]
        _print_stats("|ocl-%s|"%comp+(" "*(3+len(comp))), resid)
        _print_stats("|(ocl-%s)/%s|"%(comp,comp), relerr)

    # Plot if requested
    if '-noplot' in opts: return
    import matplotlib.pyplot as plt
    if Ncpu > 0:
        if Nocl > 0: plt.subplot(131)
        plot_theory(data, cpu, view=view)
        plt.title("%s t=%.1f ms"%(comp,cpu_time))
        cbar_title = "log I"
    if Nocl > 0:
        if Ncpu > 0: plt.subplot(132)
        plot_theory(data, ocl, view=view)
        plt.title("opencl t=%.1f ms"%ocl_time)
        cbar_title = "log I"
    if Ncpu > 0 and Nocl > 0:
        plt.subplot(133)
        if '-abs' in opts:
            err,errstr,errview = resid, "abs err", "linear"
        else:
            err,errstr,errview = abs(relerr), "rel err", "log"
        #err,errstr = ocl/cpu,"ratio"
        plot_theory(data, err, view=errview)
        plt.title("max %s = %.3g"%(errstr, max(abs(err))))
        cbar_title = errstr if errview=="linear" else "log "+errstr
    if is2D:
        h = plt.colorbar()
        h.ax.set_title(cbar_title)

    if Ncpu > 0 and Nocl > 0 and '-hist' in opts:
        plt.figure()
        v = relerr
        v[v==0] = 0.5*np.min(np.abs(v[v!=0]))
        plt.hist(np.log10(np.abs(v)), normed=1, bins=50);
        plt.xlabel('log10(err), err = | F(q) single - F(q) double| / | F(q) double |');
        plt.ylabel('P(err)')
        plt.title('Comparison of single and double precision models for %s'%name)

    plt.show()

def _print_stats(label, err):
    sorted_err = np.sort(abs(err))
    p50 = int((len(err)-1)*0.50)
    p98 = int((len(err)-1)*0.98)
    data = [
        "max:%.3e"%sorted_err[-1],
        "median:%.3e"%sorted_err[p50],
        "98%%:%.3e"%sorted_err[p98],
        "rms:%.3e"%np.sqrt(np.mean(err**2)),
        "zero-offset:%+.3e"%np.mean(err),
        ]
    print label,"  ".join(data)



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
    -lowq*/-midq/-highq/-exq use q values up to 0.05, 0.2, 1.0, 10.0
    -Nq=128 sets the number of Q points in the data set
    -1d/-2d* computes 1d or 2d data
    -preset*/-random[=seed] preset or random parameters
    -mono/-poly* force monodisperse/polydisperse
    -ctypes/-sasview* whether cpu is tested using sasview or ctypes
    -cutoff=1e-5* cutoff value for including a point in polydispersity
    -pars/-nopars* prints the parameter set or not
    -abs/-rel* plot relative or absolute error
    -linear/-log/-q4 intensity scaling
    -hist/-nohist* plot histogram of relative error
    -res=0 sets the resolution width dQ/Q if calculating with resolution
    -accuracy=Low resolution accuracy Low, Mid, High, Xhigh

Key=value pairs allow you to set specific values to any of the model
parameters.

Available models:

    %s
"""

NAME_OPTIONS = set([
    'plot','noplot',
    'single','double',
    'lowq','midq','highq','exq',
    '2d','1d',
    'preset','random',
    'poly','mono',
    'sasview','ctypes',
    'nopars','pars',
    'rel','abs',
    'linear', 'log', 'q4',
    'hist','nohist',
    ])
VALUE_OPTIONS = [
    # Note: random is both a name option and a value option
    'cutoff', 'random', 'Nq', 'res', 'accuracy',
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
