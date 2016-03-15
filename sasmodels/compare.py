#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Program to compare models using different compute engines.

This program lets you compare results between OpenCL and DLL versions
of the code and between precision (half, fast, single, double, quad),
where fast precision is single precision using native functions for
trig, etc., and may not be completely IEEE 754 compliant.  This lets
make sure that the model calculations are stable, or if you need to
tag the model as double precision only.

Run using ./compare.sh (Linux, Mac) or compare.bat (Windows) in the
sasmodels root to see the command line options.

Note that there is no way within sasmodels to select between an
OpenCL CPU device and a GPU device, but you can do so by setting the
PYOPENCL_CTX environment variable ahead of time.  Start a python
interpreter and enter::

    import pyopencl as cl
    cl.create_some_context()

This will prompt you to select from the available OpenCL devices
and tell you which string to use for the PYOPENCL_CTX variable.
On Windows you will need to remove the quotes.
"""

from __future__ import print_function

import sys
import math
import datetime
import traceback

import numpy as np

from . import core
from . import kerneldll
from . import product
from .data import plot_theory, empty_data1D, empty_data2D
from .direct_model import DirectModel
from .convert import revert_pars, constrain_new_to_old

USAGE = """
usage: compare.py model N1 N2 [options...] [key=val]

Compare the speed and value for a model between the SasView original and the
sasmodels rewrite.

model is the name of the model to compare (see below).
N1 is the number of times to run sasmodels (default=1).
N2 is the number times to run sasview (default=1).

Options (* for default):

    -plot*/-noplot plots or suppress the plot of the model
    -lowq*/-midq/-highq/-exq use q values up to 0.05, 0.2, 1.0, 10.0
    -nq=128 sets the number of Q points in the data set
    -1d*/-2d computes 1d or 2d data
    -preset*/-random[=seed] preset or random parameters
    -mono/-poly* force monodisperse/polydisperse
    -cutoff=1e-5* cutoff value for including a point in polydispersity
    -pars/-nopars* prints the parameter set or not
    -abs/-rel* plot relative or absolute error
    -linear/-log*/-q4 intensity scaling
    -hist/-nohist* plot histogram of relative error
    -res=0 sets the resolution width dQ/Q if calculating with resolution
    -accuracy=Low accuracy of the resolution calculation Low, Mid, High, Xhigh
    -edit starts the parameter explorer

Any two calculation engines can be selected for comparison:

    -single/-double/-half/-fast sets an OpenCL calculation engine
    -single!/-double!/-quad! sets an OpenMP calculation engine
    -sasview sets the sasview calculation engine

The default is -single -sasview.  Note that the interpretation of quad
precision depends on architecture, and may vary from 64-bit to 128-bit,
with 80-bit floats being common (1e-19 precision).

Key=value pairs allow you to set specific values for the model parameters.
"""

# Update docs with command line usage string.   This is separate from the usual
# doc string so that we can display it at run time if there is an error.
# lin
__doc__ = (__doc__  # pylint: disable=redefined-builtin
           + """
Program description
-------------------

"""
           + USAGE)

kerneldll.ALLOW_SINGLE_PRECISION_DLLS = True

MODELS = core.list_models()

# CRUFT python 2.6
if not hasattr(datetime.timedelta, 'total_seconds'):
    def delay(dt):
        """Return number date-time delta as number seconds"""
        return dt.days * 86400 + dt.seconds + 1e-6 * dt.microseconds
else:
    def delay(dt):
        """Return number date-time delta as number seconds"""
        return dt.total_seconds()


class push_seed(object):
    """
    Set the seed value for the random number generator.

    When used in a with statement, the random number generator state is
    restored after the with statement is complete.

    :Parameters:

    *seed* : int or array_like, optional
        Seed for RandomState

    :Example:

    Seed can be used directly to set the seed::

        >>> from numpy.random import randint
        >>> push_seed(24)
        <...push_seed object at...>
        >>> print(randint(0,1000000,3))
        [242082    899 211136]

    Seed can also be used in a with statement, which sets the random
    number generator state for the enclosed computations and restores
    it to the previous state on completion::

        >>> with push_seed(24):
        ...    print(randint(0,1000000,3))
        [242082    899 211136]

    Using nested contexts, we can demonstrate that state is indeed
    restored after the block completes::

        >>> with push_seed(24):
        ...    print(randint(0,1000000))
        ...    with push_seed(24):
        ...        print(randint(0,1000000,3))
        ...    print(randint(0,1000000))
        242082
        [242082    899 211136]
        899

    The restore step is protected against exceptions in the block::

        >>> with push_seed(24):
        ...    print(randint(0,1000000))
        ...    try:
        ...        with push_seed(24):
        ...            print(randint(0,1000000,3))
        ...            raise Exception()
        ...    except:
        ...        print("Exception raised")
        ...    print(randint(0,1000000))
        242082
        [242082    899 211136]
        Exception raised
        899
    """
    def __init__(self, seed=None):
        self._state = np.random.get_state()
        np.random.seed(seed)

    def __enter__(self):
        return None

    def __exit__(self, *args):
        np.random.set_state(self._state)

def tic():
    """
    Timer function.

    Use "toc=tic()" to start the clock and "toc()" to measure
    a time interval.
    """
    then = datetime.datetime.now()
    return lambda: delay(datetime.datetime.now() - then)


def set_beam_stop(data, radius, outer=None):
    """
    Add a beam stop of the given *radius*.  If *outer*, make an annulus.

    Note: this function does not use the sasview package
    """
    if hasattr(data, 'qx_data'):
        q = np.sqrt(data.qx_data**2 + data.qy_data**2)
        data.mask = (q < radius)
        if outer is not None:
            data.mask |= (q >= outer)
    else:
        data.mask = (data.x < radius)
        if outer is not None:
            data.mask |= (data.x >= outer)


def parameter_range(p, v):
    """
    Choose a parameter range based on parameter name and initial value.
    """
    if p.endswith('_pd_n'):
        return [0, 100]
    elif p.endswith('_pd_nsigma'):
        return [0, 5]
    elif p.endswith('_pd_type'):
        return v
    elif any(s in p for s in ('theta', 'phi', 'psi')):
        # orientation in [-180,180], orientation pd in [0,45]
        if p.endswith('_pd'):
            return [0, 45]
        else:
            return [-180, 180]
    elif 'sld' in p:
        return [-0.5, 10]
    elif p.endswith('_pd'):
        return [0, 1]
    elif p == 'background':
        return [0, 10]
    elif p == 'scale':
        return [0, 1e3]
    elif p == 'case_num':
        # RPA hack
        return [0, 10]
    elif v < 0:
        # Kxy parameters in rpa model can be negative
        return [2*v, -2*v]
    else:
        return [0, (2*v if v > 0 else 1)]


def _randomize_one(p, v):
    """
    Randomize a single parameter.
    """
    if any(p.endswith(s) for s in ('_pd_n', '_pd_nsigma', '_pd_type')):
        return v
    else:
        return np.random.uniform(*parameter_range(p, v))


def randomize_pars(pars, seed=None):
    """
    Generate random values for all of the parameters.

    Valid ranges for the random number generator are guessed from the name of
    the parameter; this will not account for constraints such as cap radius
    greater than cylinder radius in the capped_cylinder model, so
    :func:`constrain_pars` needs to be called afterward..
    """
    with push_seed(seed):
        # Note: the sort guarantees order `of calls to random number generator
        pars = dict((p, _randomize_one(p, v))
                    for p, v in sorted(pars.items()))
    return pars

def constrain_pars(model_info, pars):
    """
    Restrict parameters to valid values.

    This includes model specific code for models such as capped_cylinder
    which need to support within model constraints (cap radius more than
    cylinder radius in this case).
    """
    name = model_info['id']
    # if it is a product model, then just look at the form factor since
    # none of the structure factors need any constraints.
    if '*' in name:
        name = name.split('*')[0]

    if name == 'capped_cylinder' and pars['cap_radius'] < pars['radius']:
        pars['radius'], pars['cap_radius'] = pars['cap_radius'], pars['radius']
    if name == 'barbell' and pars['bell_radius'] < pars['radius']:
        pars['radius'], pars['bell_radius'] = pars['bell_radius'], pars['radius']

    # Limit guinier to an Rg such that Iq > 1e-30 (single precision cutoff)
    if name == 'guinier':
        #q_max = 0.2  # mid q maximum
        q_max = 1.0  # high q maximum
        rg_max = np.sqrt(90*np.log(10) + 3*np.log(pars['scale']))/q_max
        pars['rg'] = min(pars['rg'], rg_max)

    if name == 'rpa':
        # Make sure phi sums to 1.0
        if pars['case_num'] < 2:
            pars['Phia'] = 0.
            pars['Phib'] = 0.
        elif pars['case_num'] < 5:
            pars['Phia'] = 0.
        total = sum(pars['Phi'+c] for c in 'abcd')
        for c in 'abcd':
            pars['Phi'+c] /= total

def parlist(pars):
    """
    Format the parameter list for printing.
    """
    active = None
    fields = {}
    lines = []
    for k, v in sorted(pars.items()):
        parts = k.split('_pd')
        #print(k, active, parts)
        if len(parts) == 1:
            if active: lines.append(_format_par(active, **fields))
            active = k
            fields = {'value': v}
        else:
            assert parts[0] == active
            if parts[1]:
                fields[parts[1][1:]] = v
            else:
                fields['pd'] = v
    if active: lines.append(_format_par(active, **fields))
    return "\n".join(lines)

    #return "\n".join("%s: %s"%(p, v) for p, v in sorted(pars.items()))

def _format_par(name, value=0., pd=0., n=0, nsigma=3., type='gaussian'):
    line = "%s: %g"%(name, value)
    if pd != 0.  and n != 0:
        line += " +/- %g  (%d points in [-%g,%g] sigma %s)"\
                % (pd, n, nsigma, nsigma, type)
    return line

def suppress_pd(pars):
    """
    Suppress theta_pd for now until the normalization is resolved.

    May also suppress complete polydispersity of the model to test
    models more quickly.
    """
    pars = pars.copy()
    for p in pars:
        if p.endswith("_pd_n"): pars[p] = 0
    return pars

def eval_sasview(model_info, data):
    """
    Return a model calculator using the SasView fitting engine.
    """
    # importing sas here so that the error message will be that sas failed to
    # import rather than the more obscure smear_selection not imported error
    import sas
    from sas.models.qsmearing import smear_selection

    def get_model(name):
        #print("new",sorted(_pars.items()))
        sas = __import__('sas.models.' + name)
        ModelClass = getattr(getattr(sas.models, name, None), name, None)
        if ModelClass is None:
            raise ValueError("could not find model %r in sas.models"%name)
        return ModelClass()

    # grab the sasview model, or create it if it is a product model
    if model_info['composition']:
        composition_type, parts = model_info['composition']
        if composition_type == 'product':
            from sas.models import MultiplicationModel
            P, S = [get_model(p) for p in model_info['oldname']]
            model = MultiplicationModel(P, S)
        else:
            raise ValueError("mixture models not handled yet")
    else:
        model = get_model(model_info['oldname'])

    # build a smearer with which to call the model, if necessary
    smearer = smear_selection(data, model=model)
    if hasattr(data, 'qx_data'):
        q = np.sqrt(data.qx_data**2 + data.qy_data**2)
        index = ((~data.mask) & (~np.isnan(data.data))
                 & (q >= data.qmin) & (q <= data.qmax))
        if smearer is not None:
            smearer.model = model  # because smear_selection has a bug
            smearer.accuracy = data.accuracy
            smearer.set_index(index)
            theory = lambda: smearer.get_value()
        else:
            theory = lambda: model.evalDistribution([data.qx_data[index],
                                                     data.qy_data[index]])
    elif smearer is not None:
        theory = lambda: smearer(model.evalDistribution(data.x))
    else:
        theory = lambda: model.evalDistribution(data.x)

    def calculator(**pars):
        """
        Sasview calculator for model.
        """
        # paying for parameter conversion each time to keep life simple, if not fast
        pars = revert_pars(model_info, pars)
        for k, v in pars.items():
            parts = k.split('.')  # polydispersity components
            if len(parts) == 2:
                model.dispersion[parts[0]][parts[1]] = v
            else:
                model.setParam(k, v)
        return theory()

    calculator.engine = "sasview"
    return calculator

DTYPE_MAP = {
    'half': '16',
    'fast': 'fast',
    'single': '32',
    'double': '64',
    'quad': '128',
    'f16': '16',
    'f32': '32',
    'f64': '64',
    'longdouble': '128',
}
def eval_opencl(model_info, data, dtype='single', cutoff=0.):
    """
    Return a model calculator using the OpenCL calculation engine.
    """
    def builder(model_info):
        try:
            return core.build_model(model_info, dtype=dtype, platform="ocl")
        except Exception as exc:
            print(exc)
            print("... trying again with single precision")
            dtype = 'single'
            return core.build_model(model_info, dtype=dtype, platform="ocl")
    if model_info['composition']:
        composition_type, parts = model_info['composition']
        if composition_type == 'product':
            P, S = [builder(p) for p in parts]
            model = product.ProductModel(P, S)
        else:
            raise ValueError("mixture models not handled yet")
    else:
        model = builder(model_info)
    calculator = DirectModel(data, model, cutoff=cutoff)
    calculator.engine = "OCL%s"%DTYPE_MAP[dtype]
    return calculator

def eval_ctypes(model_info, data, dtype='double', cutoff=0.):
    """
    Return a model calculator using the DLL calculation engine.
    """
    if dtype == 'quad':
        dtype = 'longdouble'
    def builder(model_info):
        return core.build_model(model_info, dtype=dtype, platform="dll")

    if model_info['composition']:
        composition_type, parts = model_info['composition']
        if composition_type == 'product':
            P, S = [builder(p) for p in parts]
            model = product.ProductModel(P, S)
        else:
            raise ValueError("mixture models not handled yet")
    else:
        model = builder(model_info)
    calculator = DirectModel(data, model, cutoff=cutoff)
    calculator.engine = "OMP%s"%DTYPE_MAP[dtype]
    return calculator

def time_calculation(calculator, pars, Nevals=1):
    """
    Compute the average calculation time over N evaluations.

    An additional call is generated without polydispersity in order to
    initialize the calculation engine, and make the average more stable.
    """
    # initialize the code so time is more accurate
    value = calculator(**suppress_pd(pars))
    toc = tic()
    for _ in range(max(Nevals, 1)):  # make sure there is at least one eval
        value = calculator(**pars)
    average_time = toc()*1000./Nevals
    return value, average_time

def make_data(opts):
    """
    Generate an empty dataset, used with the model to set Q points
    and resolution.

    *opts* contains the options, with 'qmax', 'nq', 'res',
    'accuracy', 'is2d' and 'view' parsed from the command line.
    """
    qmax, nq, res = opts['qmax'], opts['nq'], opts['res']
    if opts['is2d']:
        data = empty_data2D(np.linspace(-qmax, qmax, nq), resolution=res)
        data.accuracy = opts['accuracy']
        set_beam_stop(data, 0.004)
        index = ~data.mask
    else:
        if opts['view'] == 'log':
            qmax = math.log10(qmax)
            q = np.logspace(qmax-3, qmax, nq)
        else:
            q = np.linspace(0.001*qmax, qmax, nq)
        data = empty_data1D(q, resolution=res)
        index = slice(None, None)
    return data, index

def make_engine(model_info, data, dtype, cutoff):
    """
    Generate the appropriate calculation engine for the given datatype.

    Datatypes with '!' appended are evaluated using external C DLLs rather
    than OpenCL.
    """
    if dtype == 'sasview':
        return eval_sasview(model_info, data)
    elif dtype.endswith('!'):
        return eval_ctypes(model_info, data, dtype=dtype[:-1], cutoff=cutoff)
    else:
        return eval_opencl(model_info, data, dtype=dtype, cutoff=cutoff)

def compare(opts, limits=None):
    """
    Preform a comparison using options from the command line.

    *limits* are the limits on the values to use, either to set the y-axis
    for 1D or to set the colormap scale for 2D.  If None, then they are
    inferred from the data and returned. When exploring using Bumps,
    the limits are set when the model is initially called, and maintained
    as the values are adjusted, making it easier to see the effects of the
    parameters.
    """
    Nbase, Ncomp = opts['n1'], opts['n2']
    pars = opts['pars']
    data = opts['data']

    # Base calculation
    if Nbase > 0:
        base = opts['engines'][0]
        try:
            base_value, base_time = time_calculation(base, pars, Nbase)
            print("%s t=%.1f ms, intensity=%.0f"
                  % (base.engine, base_time, sum(base_value)))
        except ImportError:
            traceback.print_exc()
            Nbase = 0

    # Comparison calculation
    if Ncomp > 0:
        comp = opts['engines'][1]
        try:
            comp_value, comp_time = time_calculation(comp, pars, Ncomp)
            print("%s t=%.1f ms, intensity=%.0f"
                  % (comp.engine, comp_time, sum(comp_value)))
        except ImportError:
            traceback.print_exc()
            Ncomp = 0

    # Compare, but only if computing both forms
    if Nbase > 0 and Ncomp > 0:
        resid = (base_value - comp_value)
        relerr = resid/comp_value
        _print_stats("|%s-%s|"
                     % (base.engine, comp.engine) + (" "*(3+len(comp.engine))),
                     resid)
        _print_stats("|(%s-%s)/%s|"
                     % (base.engine, comp.engine, comp.engine),
                     relerr)

    # Plot if requested
    if not opts['plot'] and not opts['explore']: return
    view = opts['view']
    import matplotlib.pyplot as plt
    if limits is None:
        vmin, vmax = np.Inf, -np.Inf
        if Nbase > 0:
            vmin = min(vmin, min(base_value))
            vmax = max(vmax, max(base_value))
        if Ncomp > 0:
            vmin = min(vmin, min(comp_value))
            vmax = max(vmax, max(comp_value))
        limits = vmin, vmax

    if Nbase > 0:
        if Ncomp > 0: plt.subplot(131)
        plot_theory(data, base_value, view=view, use_data=False, limits=limits)
        plt.title("%s t=%.1f ms"%(base.engine, base_time))
        #cbar_title = "log I"
    if Ncomp > 0:
        if Nbase > 0: plt.subplot(132)
        plot_theory(data, comp_value, view=view, use_data=False, limits=limits)
        plt.title("%s t=%.1f ms"%(comp.engine, comp_time))
        #cbar_title = "log I"
    if Ncomp > 0 and Nbase > 0:
        plt.subplot(133)
        if not opts['rel_err']:
            err, errstr, errview = resid, "abs err", "linear"
        else:
            err, errstr, errview = abs(relerr), "rel err", "log"
        #err,errstr = base/comp,"ratio"
        plot_theory(data, None, resid=err, view=errview, use_data=False)
        if view == 'linear':
            plt.xscale('linear')
        plt.title("max %s = %.3g"%(errstr, max(abs(err))))
        #cbar_title = errstr if errview=="linear" else "log "+errstr
    #if is2D:
    #    h = plt.colorbar()
    #    h.ax.set_title(cbar_title)

    if Ncomp > 0 and Nbase > 0 and '-hist' in opts:
        plt.figure()
        v = relerr
        v[v == 0] = 0.5*np.min(np.abs(v[v != 0]))
        plt.hist(np.log10(np.abs(v)), normed=1, bins=50)
        plt.xlabel('log10(err), err = |(%s - %s) / %s|'
                   % (base.engine, comp.engine, comp.engine))
        plt.ylabel('P(err)')
        plt.title('Distribution of relative error between calculation engines')

    if not opts['explore']:
        plt.show()

    return limits

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
    print(label+"  "+"  ".join(data))



# ===========================================================================
#
NAME_OPTIONS = set([
    'plot', 'noplot',
    'half', 'fast', 'single', 'double',
    'single!', 'double!', 'quad!', 'sasview',
    'lowq', 'midq', 'highq', 'exq',
    '2d', '1d',
    'preset', 'random',
    'poly', 'mono',
    'nopars', 'pars',
    'rel', 'abs',
    'linear', 'log', 'q4',
    'hist', 'nohist',
    'edit',
    ])
VALUE_OPTIONS = [
    # Note: random is both a name option and a value option
    'cutoff', 'random', 'nq', 'res', 'accuracy',
    ]

def columnize(L, indent="", width=79):
    """
    Format a list of strings into columns.

    Returns a string with carriage returns ready for printing.
    """
    column_width = max(len(w) for w in L) + 1
    num_columns = (width - len(indent)) // column_width
    num_rows = len(L) // num_columns
    L = L + [""] * (num_rows*num_columns - len(L))
    columns = [L[k*num_rows:(k+1)*num_rows] for k in range(num_columns)]
    lines = [" ".join("%-*s"%(column_width, entry) for entry in row)
             for row in zip(*columns)]
    output = indent + ("\n"+indent).join(lines)
    return output


def get_demo_pars(model_info):
    """
    Extract demo parameters from the model definition.
    """
    # Get the default values for the parameters
    pars = dict((p[0], p[2]) for p in model_info['parameters'])

    # Fill in default values for the polydispersity parameters
    for p in model_info['parameters']:
        if p[4] in ('volume', 'orientation'):
            pars[p[0]+'_pd'] = 0.0
            pars[p[0]+'_pd_n'] = 0
            pars[p[0]+'_pd_nsigma'] = 3.0
            pars[p[0]+'_pd_type'] = "gaussian"

    # Plug in values given in demo
    pars.update(model_info['demo'])
    return pars


def parse_opts():
    """
    Parse command line options.
    """
    MODELS = core.list_models()
    flags = [arg for arg in sys.argv[1:]
             if arg.startswith('-')]
    values = [arg for arg in sys.argv[1:]
              if not arg.startswith('-') and '=' in arg]
    args = [arg for arg in sys.argv[1:]
            if not arg.startswith('-') and '=' not in arg]
    models = "\n    ".join("%-15s"%v for v in MODELS)
    if len(args) == 0:
        print(USAGE)
        print("\nAvailable models:")
        print(columnize(MODELS, indent="  "))
        sys.exit(1)
    if len(args) > 3:
        print("expected parameters: model N1 N2")

    def load_model(name):
        try:
            model_info = core.load_model_info(name)
        except ImportError, exc:
            print(str(exc))
            print("Use one of:\n    " + models)
            sys.exit(1)
        return model_info

    name = args[0]
    if '*' in name:
        parts = [load_model(k) for k in name.split('*')]
        model_info = product.make_product_info(*parts)
    else:
        model_info = load_model(name)

    invalid = [o[1:] for o in flags
               if o[1:] not in NAME_OPTIONS
               and not any(o.startswith('-%s='%t) for t in VALUE_OPTIONS)]
    if invalid:
        print("Invalid options: %s"%(", ".join(invalid)))
        sys.exit(1)


    # pylint: disable=bad-whitespace
    # Interpret the flags
    opts = {
        'plot'      : True,
        'view'      : 'log',
        'is2d'      : False,
        'qmax'      : 0.05,
        'nq'        : 128,
        'res'       : 0.0,
        'accuracy'  : 'Low',
        'cutoff'    : 1e-5,
        'seed'      : -1,  # default to preset
        'mono'      : False,
        'show_pars' : False,
        'show_hist' : False,
        'rel_err'   : True,
        'explore'   : False,
    }
    engines = []
    for arg in flags:
        if arg == '-noplot':    opts['plot'] = False
        elif arg == '-plot':    opts['plot'] = True
        elif arg == '-linear':  opts['view'] = 'linear'
        elif arg == '-log':     opts['view'] = 'log'
        elif arg == '-q4':      opts['view'] = 'q4'
        elif arg == '-1d':      opts['is2d'] = False
        elif arg == '-2d':      opts['is2d'] = True
        elif arg == '-exq':     opts['qmax'] = 10.0
        elif arg == '-highq':   opts['qmax'] = 1.0
        elif arg == '-midq':    opts['qmax'] = 0.2
        elif arg == '-lowq':    opts['qmax'] = 0.05
        elif arg.startswith('-nq='):       opts['nq'] = int(arg[4:])
        elif arg.startswith('-res='):      opts['res'] = float(arg[5:])
        elif arg.startswith('-accuracy='): opts['accuracy'] = arg[10:]
        elif arg.startswith('-cutoff='):   opts['cutoff'] = float(arg[8:])
        elif arg.startswith('-random='):   opts['seed'] = int(arg[8:])
        elif arg == '-random':  opts['seed'] = np.random.randint(1e6)
        elif arg == '-preset':  opts['seed'] = -1
        elif arg == '-mono':    opts['mono'] = True
        elif arg == '-poly':    opts['mono'] = False
        elif arg == '-pars':    opts['show_pars'] = True
        elif arg == '-nopars':  opts['show_pars'] = False
        elif arg == '-hist':    opts['show_hist'] = True
        elif arg == '-nohist':  opts['show_hist'] = False
        elif arg == '-rel':     opts['rel_err'] = True
        elif arg == '-abs':     opts['rel_err'] = False
        elif arg == '-half':    engines.append(arg[1:])
        elif arg == '-fast':    engines.append(arg[1:])
        elif arg == '-single':  engines.append(arg[1:])
        elif arg == '-double':  engines.append(arg[1:])
        elif arg == '-single!': engines.append(arg[1:])
        elif arg == '-double!': engines.append(arg[1:])
        elif arg == '-quad!':   engines.append(arg[1:])
        elif arg == '-sasview': engines.append(arg[1:])
        elif arg == '-edit':    opts['explore'] = True
    # pylint: enable=bad-whitespace

    if len(engines) == 0:
        engines.extend(['single', 'sasview'])
    elif len(engines) == 1:
        if engines[0][0] != 'sasview':
            engines.append('sasview')
        else:
            engines.append('single')
    elif len(engines) > 2:
        del engines[2:]

    n1 = int(args[1]) if len(args) > 1 else 1
    n2 = int(args[2]) if len(args) > 2 else 1

    # Get demo parameters from model definition, or use default parameters
    # if model does not define demo parameters
    pars = get_demo_pars(model_info)

    # Fill in parameters given on the command line
    presets = {}
    for arg in values:
        k, v = arg.split('=', 1)
        if k not in pars:
            # extract base name without polydispersity info
            s = set(p.split('_pd')[0] for p in pars)
            print("%r invalid; parameters are: %s"%(k, ", ".join(sorted(s))))
            sys.exit(1)
        presets[k] = float(v) if not k.endswith('type') else v

    # randomize parameters
    #pars.update(set_pars)  # set value before random to control range
    if opts['seed'] > -1:
        pars = randomize_pars(pars, seed=opts['seed'])
        print("Randomize using -random=%i"%opts['seed'])
    if opts['mono']:
        pars = suppress_pd(pars)
    pars.update(presets)  # set value after random to control value
    constrain_pars(model_info, pars)
    constrain_new_to_old(model_info, pars)
    if opts['show_pars']:
        print(str(parlist(pars)))

    # Create the computational engines
    data, _ = make_data(opts)
    if n1:
        base = make_engine(model_info, data, engines[0], opts['cutoff'])
    else:
        base = None
    if n2:
        comp = make_engine(model_info, data, engines[1], opts['cutoff'])
    else:
        comp = None

    # pylint: disable=bad-whitespace
    # Remember it all
    opts.update({
        'name'      : name,
        'def'       : model_info,
        'n1'        : n1,
        'n2'        : n2,
        'presets'   : presets,
        'pars'      : pars,
        'data'      : data,
        'engines'   : [base, comp],
    })
    # pylint: enable=bad-whitespace

    return opts

def explore(opts):
    """
    Explore the model using the Bumps GUI.
    """
    import wx
    from bumps.names import FitProblem
    from bumps.gui.app_frame import AppFrame

    problem = FitProblem(Explore(opts))
    is_mac = "cocoa" in wx.version()
    app = wx.App()
    frame = AppFrame(parent=None, title="explore")
    if not is_mac: frame.Show()
    frame.panel.set_model(model=problem)
    frame.panel.Layout()
    frame.panel.aui.Split(0, wx.TOP)
    if is_mac: frame.Show()
    app.MainLoop()

class Explore(object):
    """
    Bumps wrapper for a SAS model comparison.

    The resulting object can be used as a Bumps fit problem so that
    parameters can be adjusted in the GUI, with plots updated on the fly.
    """
    def __init__(self, opts):
        from bumps.cli import config_matplotlib
        from . import bumps_model
        config_matplotlib()
        self.opts = opts
        model_info = opts['def']
        pars, pd_types = bumps_model.create_parameters(model_info, **opts['pars'])
        if not opts['is2d']:
            active = [base + ext
                      for base in model_info['partype']['pd-1d']
                      for ext in ['', '_pd', '_pd_n', '_pd_nsigma']]
            active.extend(model_info['partype']['fixed-1d'])
            for k in active:
                v = pars[k]
                v.range(*parameter_range(k, v.value))
        else:
            for k, v in pars.items():
                v.range(*parameter_range(k, v.value))

        self.pars = pars
        self.pd_types = pd_types
        self.limits = None

    def numpoints(self):
        """
        Return the number of points.
        """
        return len(self.pars) + 1  # so dof is 1

    def parameters(self):
        """
        Return a dictionary of parameters.
        """
        return self.pars

    def nllf(self):
        """
        Return cost.
        """
        # pylint: disable=no-self-use
        return 0.  # No nllf

    def plot(self, view='log'):
        """
        Plot the data and residuals.
        """
        pars = dict((k, v.value) for k, v in self.pars.items())
        pars.update(self.pd_types)
        self.opts['pars'] = pars
        limits = compare(self.opts, limits=self.limits)
        if self.limits is None:
            vmin, vmax = limits
            vmax = 1.3*vmax
            vmin = vmax*1e-7
            self.limits = vmin, vmax


def main():
    """
    Main program.
    """
    opts = parse_opts()
    if opts['explore']:
        explore(opts)
    else:
        compare(opts)

if __name__ == "__main__":
    main()
