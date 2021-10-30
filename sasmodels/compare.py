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

Run using "./sascomp -h" in the sasmodels root to see the command
line options. To run from from an installed version of sasmodels,
use "python -m sasmodels.compare -h".

Note that there is no way within sasmodels to select between an
OpenCL CPU device and a GPU device, but you can do so by setting the
SAS_OPENCL environment variable. Start a python interpreter and enter::

    import pyopencl as cl
    cl.create_some_context()

This will prompt you to select from the available OpenCL devices
and tell you which string to use for the SAS_OPENCL variable.
On Windows you will need to remove the quotes.
"""

from __future__ import print_function, division

import sys
import os
import math
import datetime
import traceback
import re

import numpy as np  # type: ignore

from . import core
from . import weights
from . import kerneldll
from . import kernelcl
from . import kernelcuda
from .data import plot_theory, empty_data1D, empty_data2D, load_data
from .direct_model import DirectModel, get_mesh
from .generate import FLOAT_RE, set_integration_size

# pylint: disable=unused-import
from typing import Optional, Dict, Any, Callable, Tuple, List
from .modelinfo import ModelInfo, Parameter, ParameterSet
from .data import Data

try:
    # With python 3.8+ we can indicate that calculator takes floats.
    from typing import Protocol
    class Calculator(Protocol):
        """Kernel calculator takes *par=value* keyword arguments."""
        def __call__(self, **par: float) -> np.ndarray: ...
except ImportError:
    #: Kernel calculator takes *par=value* keyword arguments.
    Calculator = Callable[..., np.ndarray]
# pylint: enable=unused-import

USAGE = """
usage: sascomp model [options...] [key=val]

Generate and compare SAS models.  If a single model is specified it shows
a plot of that model.  Different models can be compared, or the same model
with different parameters.  The same model with the same parameters can
be compared with different calculation engines to see the effects of precision
on the resultant values.

model or model1,model2 are the names of the models to compare (see below).

Options (* for default):

    === data generation ===
    -data="path" uses q, dq from the data file
    -noise=0 sets the measurement error dI/I
    -res=0 sets the resolution width dQ/Q if calculating with resolution
    -lowq*/-midq/-highq/-exq use q values up to 0.05, 0.2, 1.0, 10.0
    -q=min:max alternative specification of qrange
    -nq=128 sets the number of Q points in the data set
    -1d*/-2d computes 1d or 2d data
    -zero indicates that q=0 should be included

    === model parameters ===
    -preset*/-random[=seed] preset or random parameters
    -sets=n generates n random datasets with the seed given by -random=seed
    -pars/-nopars* prints the parameter set or not
    -sphere[=150] set up spherical integration over theta/phi using n points
    -mono*/-poly suppress or allow polydispersity on generated parameters
    -magnetic/-nonmagnetic* suppress or allow magnetism on generated parameters
    -maxdim[=inf] limit randomly generate particle dimensions to maxdim

    === calculation options ===
    -cutoff=1e-5* cutoff value for including a point in polydispersity
    -accuracy=Low accuracy of the resolution calculation Low, Mid, High, Xhigh
    -neval=1 sets the number of evals for more accurate timing
    -ngauss=0 overrides the number of points in the 1-D gaussian quadrature

    === precision options ===
    -engine=default uses the default calcution precision
    -single/-double/-half/-fast sets an OpenCL calculation engine
    -single!/-double!/-quad! sets an OpenMP calculation engine

    === plotting ===
    -plot*/-noplot plots or suppress the plot of the model
    -linear/-log*/-q4 intensity scaling on plots
    -hist/-nohist* plot histogram of relative error
    -abs/-rel* plot relative or absolute error
    -title="note" adds note to the plot title, after the model name
    -weights shows weights plots for the polydisperse parameters
    -profile shows the sld profile if the model has a plottable sld profile

    === output options ===
    -edit starts the parameter explorer
    -help/-html shows the model docs instead of running the model

    === help ===
    -h/-? print this help
    -models[=all] show all builtin models of a given type:
        all, py, c, double, single, opencl, 1d, 2d, magnetic

    === environment variables ===
    -DSAS_MODELPATH=~/.sasmodels/custom_models sets path to custom models
    -DSAS_WEIGHTS_PATH=~/.sasview/weights sets path to custom distributions
    -DSAS_OPENCL=vendor:device|cuda:device|none sets the target GPU device
    -DXDG_CACHE_HOME=~/.cache sets the pyopencl cache root (linux only)
    -DSAS_COMPILER=tinycc|msvc|mingw|unix sets the DLL compiler
    -DSAS_OPENMP=0 set to 1 to turn on OpenMP for the DLLs
    -DSAS_DLL_PATH=~/.sasmodels/compiled_models sets the DLL cache
    -DPYOPENCL_NO_CACHE=1 turns off the PyOpenCL cache

The interpretation of quad precision depends on architecture, and may
vary from 64-bit to 128-bit, with 80-bit floats being common (1e-19 precision).
On unix and mac you may need single quotes around the DLL computation
engines, such as -engine='single!,double!' since !, is treated as a history
expansion request in the shell.

Key=value pairs allow you to set specific values for the model parameters.
Key=value1,value2 to compare different values of the same parameter. The
value can be an expression including other parameters.

Items later on the command line override those that appear earlier.

Examples:

    # compare single and double precision calculation for a barbell
    sascomp barbell -engine=single,double

    # generate 10 random lorentz models, with seed=27
    sascomp lorentz -sets=10 -seed=27

    # compare ellipsoid with R = R_polar = R_equatorial to sphere of radius R
    sascomp sphere,ellipsoid radius_polar=radius radius_equatorial=radius

    # model timing test requires multiple evals to perform the estimate
    sascomp pringle -engine=single,double -timing=100,100 -noplot
"""

kerneldll.ALLOW_SINGLE_PRECISION_DLLS = True

def build_math_context():
    # type: () -> Dict[str, Callable]
    """build dictionary of functions from math module"""
    return dict((k, getattr(math, k))
                for k in dir(math) if not k.startswith('_'))

#: list of math functions for use in evaluating parameters
MATH = build_math_context()

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
        ...    except Exception:
        ...        print("Exception raised")
        ...    print(randint(0,1000000))
        242082
        [242082    899 211136]
        Exception raised
        899
    """
    def __init__(self, seed=None):
        # type: (Optional[int]) -> None
        self._state = np.random.get_state()
        np.random.seed(seed)

    def __enter__(self):
        # type: () -> None
        pass

    def __exit__(self, exc_type, exc_value, trace):
        # type: (Any, BaseException, Any) -> None
        np.random.set_state(self._state)

def tic():
    # type: () -> Callable[[], float]
    """
    Timer function.

    Use "toc=tic()" to start the clock and "toc()" to measure
    a time interval.
    """
    then = datetime.datetime.now()
    return lambda: delay(datetime.datetime.now() - then)


def set_beam_stop(data, radius, outer=None):
    # type: (Data, float, float) -> None
    """
    Add a beam stop of the given *radius*.  If *outer*, make an annulus.
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
    # type: (str, float) -> Tuple[float, float]
    """
    Choose a parameter range based on parameter name and initial value.
    """
    # process the polydispersity options
    if p.endswith('_pd_n'):
        return 0., 100.
    elif p.endswith('_pd_nsigma'):
        return 0., 5.
    elif p.endswith('_pd_type'):
        raise ValueError("Cannot return a range for a string value")
    elif any(s in p for s in ('theta', 'phi', 'psi')):
        # orientation in [-180,180], orientation pd in [0,45]
        if p.endswith('_pd'):
            return 0., 180.
        else:
            return -180., 180.
    elif p.endswith('_pd'):
        return 0., 1.
    elif 'sld' in p:
        return -0.5, 10.
    elif p == 'background':
        return 0., 10.
    elif p == 'scale':
        return 0., 1.e3
    elif v < 0.:
        return 2.*v, -2.*v
    else:
        return 0., (2.*v if v > 0. else 1.)


def _randomize_one(model_info, name, value):
    # type: (ModelInfo, str, float) -> float
    """
    Randomize a single parameter.
    """
    # Set the amount of polydispersity/angular dispersion, but by default pd_n
    # is zero so there is no polydispersity.  This allows us to turn on/off
    # pd by setting pd_n, and still have randomly generated values
    if name.endswith('_pd'):
        par = model_info.parameters[name[:-3]]
        if par.type == 'orientation':
            # Let oriention variation peak around 13 degrees; 95% < 42 degrees
            return 180*np.random.beta(2.5, 20)
        else:
            # Let polydispersity peak around 15%; 95% < 0.4; max=100%
            return np.random.beta(1.5, 7)

    # pd is selected globally rather than per parameter, so set to 0 for no pd
    # In particular, when multiple pd dimensions, want to decrease the number
    # of points per dimension for faster computation
    if name.endswith('_pd_n'):
        return 0

    # Don't mess with distribution type for now
    if name.endswith('_pd_type'):
        return 'gaussian'

    # type-dependent value of number of sigmas; for gaussian use 3.
    if name.endswith('_pd_nsigma'):
        return 3.

    # background in the range [0.01, 1]
    if name == 'background':
        return 10**np.random.uniform(-2, 0)

    # scale defaults to 0.1% to 30% volume fraction
    if name == 'scale':
        return 10**np.random.uniform(-3, -0.5)

    # If it is a list of choices, pick one at random with equal probability
    par = model_info.parameters[name]
    if par.choices:  # choice list
        return np.random.randint(len(par.choices))

    # If it is a fixed range, pick from it with equal probability.
    # For logarithmic ranges, the model will have to override.
    if np.isfinite(par.limits).all():
        return np.random.uniform(*par.limits)

    # If the paramter is marked as an sld use the range of neutron slds
    if par.type == 'sld':
        return np.random.uniform(-0.5, 12)

    # Limit magnetic SLDs to a smaller range, from zero to iron=5/A^2
    if par.name.endswith('_M0'):
        return np.random.uniform(0, 5)

    # Guess at the random length/radius/thickness.  In practice, all models
    # are going to set their own reasonable ranges.
    if par.type == 'volume':
        if ('length' in par.name or
                'radius' in par.name or
                'thick' in par.name):
            return 10**np.random.uniform(2, 4)

    # In the absence of any other info, select a value in [0, 2v], or
    # [-2|v|, 2|v|] if v is negative, or [0, 1] if v is zero.  Mostly the
    # model random parameter generators will override this default.
    low, high = parameter_range(par.name, value)
    limits = (max(par.limits[0], low), min(par.limits[1], high))
    return np.random.uniform(*limits)

def _random_pd(model_info, pars, is2d):
    # type: (ModelInfo, Dict[str, float], bool) -> None
    """
    Generate a random dispersity distribution for the model.

    1% no shape dispersity
    85% single shape parameter
    13% two shape parameters
    1% three shape parameters

    If oriented, then put dispersity in theta, add phi and psi dispersity
    with 10% probability for each.
    """
    # Find the polydisperse parameters.
    pd = [p for p in model_info.parameters.kernel_parameters if p.polydisperse]

    # If the sample is oriented then add polydispersity to the orientation.
    oriented = any(p.type == 'orientation' for p in pd)
    num_oriented_pd = 0
    if oriented:
        if np.random.rand() < 0.8:
            # 80% change of pd on long axis (20x cost)
            pars['theta_pd_n'] = 20
            num_oriented_pd += 1
        if np.random.rand() < 0.1:
            # 10% change of pd on short axis (5x cost)
            pars['phi_pd_n'] = 5
            num_oriented_pd += 1
        if any(p.name == 'psi' for p in pd) and np.random.rand() < 0.1:
            # 10% change of pd on spin axis (5x cost)
            #print("generating psi_pd_n")
            pars['psi_pd_n'] = 5
            num_oriented_pd += 1

    # Process non-orientation parameters
    pd = [p for p in pd if p.type != 'orientation']

    # Find the remaining pd parameters, which are all volume parameters.
    # Use the parameter value as the weight on the choice function for
    # the polydispersity parameter. The I(Q) curve is more sensitive to
    # pd on larger dimensions, so they should be preferred.
    # TODO: choose better weights for parameters like num_pearls or num_disks.
    name = [] # type: List[str]  # name of the next volume parameter
    default = [] # type: List[float]  # default val for that volume parameter
    for p in pd:
        if p.length_control is not None:
            slots = int(pars.get(p.length_control, 1) + 0.5)
            name.extend(p.name+str(k+1) for k in range(slots))
            default.extend(p.default for k in range(slots))
        elif p.length > 1:
            slots = p.length
            name.extend(p.name+str(k+1) for k in range(slots))
            default.extend(p.default for k in range(slots))
        else:
            name.append(p.name)
            default.append(p.default)
    p = [pars.get(k, v) for k, v in zip(name, default)]  # relative weight
    p = np.array(p)/sum(p) if p else []  # normalize to probability

    # Select number of pd parameters to use.  The selection is biased
    # toward fewer pd parameters if there is already orientational pd
    # (effectively allowing only one volume pd) and the number of pd steps
    # is scaled down.  Ignore oriented if it is not 2d data.
    if not is2d:
        num_oriented_pd = 0
    n = len(name)
    u = np.random.rand()
    if u < (1 - 1/(1+num_oriented_pd)):
        # if lots of orientation dispersity then reject shape dispersity
        pass
    elif u < 0.01 or n < 1:
        # 1% chance of no polydispersity (1x cost)
        pass
    elif u < 0.66 or n < 2:
        # 65% chance of pd on one value (35x cost)
        choice = np.random.choice(n, size=1, replace=False, p=p)
        pars[name[choice[0]]+"_pd_n"] = 35 // (1 + num_oriented_pd)
    elif u < 0.99 or n < 3:
        # 33% chance of pd on two values (250x cost)
        choice = np.random.choice(n, size=2, replace=False, p=p)
        pars[name[choice[0]]+"_pd_n"] = 25 // (1 + num_oriented_pd)
        pars[name[choice[1]]+"_pd_n"] = 10 // (1 + num_oriented_pd)
    else:
        # 1% chance of pd on three values (1250x cost)
        choice = np.random.choice(n, size=3, replace=False, p=p)
        pars[name[choice[0]]+"_pd_n"] = 25
        pars[name[choice[1]]+"_pd_n"] = 10
        pars[name[choice[2]]+"_pd_n"] = 5


    ## Show selected polydispersity
    #for name, value in pars.items():
    #    if name.endswith('_pd_n') and value > 0:
    #        print(name, value, pars.get(name[:-5], 0), pars.get(name[:-2], 0))


def randomize_pars(model_info, pars, maxdim=np.inf, is2d=False):
    # type: (ModelInfo, ParameterSet, float, bool) -> ParameterSet
    """
    Generate random values for all of the parameters.

    Valid ranges for the random number generator are guessed from the name of
    the parameter; this will not account for constraints such as cap radius
    greater than cylinder radius in the capped_cylinder model, so
    :func:`constrain_pars` needs to be called afterward..
    """
    # Note: the sort guarantees order of calls to random number generator
    random_pars = dict((p, _randomize_one(model_info, p, v))
                       for p, v in sorted(pars.items()))
    if model_info.random is not None:
        random_pars.update(model_info.random())
    _random_pd(model_info, random_pars, is2d)
    limit_dimensions(model_info, random_pars, maxdim)
    return random_pars


def limit_dimensions(model_info, pars, maxdim):
    # type: (ModelInfo, ParameterSet, float) -> None
    """
    Limit parameters of units of Ang to maxdim.
    """
    for p in model_info.parameters.call_parameters:
        value = pars[p.name]
        if p.units == 'Ang' and value > maxdim:
            pars[p.name] = maxdim*10**np.random.uniform(-3, 0)

def _swap_pars(pars, a, b):
    # type: (ModelInfo, str, str) -> None
    """
    Swap polydispersity and magnetism when swapping parameters.

    Assume the parameters are of the same basic type (volume, sld, or other),
    so that if, for example, radius_pd is in pars but radius_bell_pd is not,
    then after the swap radius_bell_pd will be the old radius_pd and radius_pd
    will be removed.
    """
    for ext in ("", "_pd", "_pd_n", "_pd_nsigma", "_pd_type", "_M0", "_mphi", "_mtheta"):
        ax, bx = a+ext, b+ext
        if ax in pars and bx in pars:
            pars[ax], pars[bx] = pars[bx], pars[ax]
        elif ax in pars:
            pars[bx] = pars[ax]
            del pars[ax]
        elif bx in pars:
            pars[ax] = pars[bx]
            del pars[bx]

def constrain_pars(model_info, pars):
    # type: (ModelInfo, ParameterSet) -> None
    """
    Restrict parameters to valid values.

    This includes model specific code for models such as capped_cylinder
    which need to support within model constraints (cap radius more than
    cylinder radius in this case).

    Warning: this updates the *pars* dictionary in place.
    """
    # TODO: move the model specific code to the individual models
    name = model_info.id
    # if it is a product model, then just look at the form factor since
    # none of the structure factors need any constraints.
    if '*' in name:
        name = name.split('*')[0]

    # Suppress magnetism for python models (not yet implemented)
    if callable(model_info.Iq):
        pars.update(suppress_magnetism(pars))

    if name == 'barbell':
        if pars['radius_bell'] < pars['radius']:
            _swap_pars(pars, 'radius_bell', 'radius')

    elif name == 'capped_cylinder':
        if pars['radius_cap'] < pars['radius']:
            _swap_pars(pars, 'radius_cap', 'radius')

    elif name == 'guinier':
        # Limit guinier to an Rg such that Iq > 1e-30 (single precision cutoff)
        # I(q) = A e^-(Rg^2 q^2/3) > e^-(30 ln 10)
        # => ln A - (Rg^2 q^2/3) > -30 ln 10
        # => Rg^2 q^2/3 < 30 ln 10 + ln A
        # => Rg < sqrt(90 ln 10 + 3 ln A)/q
        #q_max = 0.2  # mid q maximum
        q_max = 1.0  # high q maximum
        rg_max = np.sqrt(90*np.log(10) + 3*np.log(pars['scale']))/q_max
        pars['rg'] = min(pars['rg'], rg_max)

    elif name == 'pearl_necklace':
        if pars['radius'] < pars['thick_string']:
            _swap_pars(pars, 'thick_string', 'radius')

    elif name == 'rpa':
        # Make sure phi sums to 1.0
        if pars['case_num'] < 2:
            pars['Phi1'] = 0.
            pars['Phi2'] = 0.
        elif pars['case_num'] < 5:
            pars['Phi1'] = 0.
        total = sum(pars['Phi'+c] for c in '1234')
        for c in '1234':
            pars['Phi'+c] /= total

def parlist(model_info, pars, is2d):
    # type: (ModelInfo, ParameterSet, bool) -> str
    """
    Format the parameter list for printing.
    """
    lines = []
    parameters = model_info.parameters
    magnetic = False
    magnetic_pars = []
    for p in parameters.user_parameters(pars, True):
        if any(p.id.endswith(x) for x in ('_M0', '_mtheta', '_mphi')):
            continue
        if p.id in set(('up_frac_i', 'up_frac_f', 'up_angle', 'up_phi')):
            magnetic_pars.append("%s=%s"%(p.id, pars.get(p.id, p.default)))
            continue
        if not is2d and p.id in ('theta', 'phi', 'psi'):
            continue
        fields = dict(
            value=pars.get(p.id, p.default),
            pd=pars.get(p.id+"_pd", 0.),
            n=int(pars.get(p.id+"_pd_n", 0)),
            nsigma=pars.get(p.id+"_pd_nsgima", 3.),
            pdtype=pars.get(p.id+"_pd_type", 'gaussian'),
            relative_pd=p.relative_pd,
            M0=pars.get(p.id+'_M0', 0.),
            mphi=pars.get(p.id+'_mphi', 0.),
            mtheta=pars.get(p.id+'_mtheta', 0.),
        )
        lines.append(_format_par(p.name, **fields))
        magnetic = magnetic or fields['M0'] != 0.
    if magnetic and magnetic_pars:
        lines.append(" ".join(magnetic_pars))
    return "\n".join(lines)

    #return "\n".join("%s: %s"%(p, v) for p, v in sorted(pars.items()))

def _format_par(name, value=0., pd=0., n=0, nsigma=3., pdtype='gaussian',
                relative_pd=False, M0=0., mphi=0., mtheta=0.):
    # type: (str, float, float, int, float, str, bool, float, float, float) -> str
    line = "%s: %g"%(name, value)
    if pd != 0.  and n != 0:
        if relative_pd:
            pd *= value
        line += " +/- %g  (%d points in [-%g,%g] sigma %s)"\
                % (pd, n, nsigma, nsigma, pdtype)
    if M0 != 0.:
        line += "  M0:%.3f  mtheta:%.1f  mphi:%.1f" % (M0, mtheta, mphi)
    return line

def suppress_pd(pars):
    # type: (ParameterSet) -> ParameterSet
    """
    Complete eliminate polydispersity of the model to test models more quickly.
    """
    pars = pars.copy()
    for p in pars:
        if p.endswith("_pd_n"):
            pars[p] = 0
    return pars

def suppress_magnetism(pars):
    # type: (ParameterSet) -> ParameterSet
    """
    Complete eliminate magnetism of the model to test models more quickly.
    """
    pars = pars.copy()
    for p in pars:
        if p.endswith("_M0"):
            pars[p] = 0
    return pars

def time_calculation(calculator: Calculator, pars: ParameterSet, evals: int=1):
    # not type: (Calculator, ParameterSet, int) -> Tuple[np.ndarray, float]
    """
    Compute the average calculation time over N evaluations.

    An additional call is generated without polydispersity in order to
    initialize the calculation engine, and make the average more stable.
    """
    # initialize the code so time is more accurate
    if evals > 1:
        calculator(**suppress_pd(pars))
    toc = tic()
    # make sure there is at least one eval
    value = calculator(**pars)
    for _ in range(evals-1):
        value = calculator(**pars)
    average_time = toc()*1000. / evals
    #print("I(q)",value)
    return value, average_time

def make_data(opts):
    # type: (Dict[str, Any]) -> Tuple[Data, np.ndarray]
    """
    Generate an empty dataset, used with the model to set Q points
    and resolution.

    *opts* contains the options, with 'qmax', 'nq', 'res',
    'accuracy', 'is2d' and 'view' parsed from the command line.
    """
    qmin, qmax, nq, res = opts['qmin'], opts['qmax'], opts['nq'], opts['res']
    if opts['is2d']:
        q = np.linspace(-qmax, qmax, nq)  # type: np.ndarray
        data = empty_data2D(q, resolution=res)
        data.accuracy = opts['accuracy']
        set_beam_stop(data, qmin)
        index = ~data.mask
    else:
        if opts['view'] == 'log' and not opts['zero']:
            q = np.logspace(math.log10(qmin), math.log10(qmax), nq)
        else:
            q = np.linspace(qmin, qmax, nq)
        if opts['zero']:
            q = np.hstack((0, q))
        # TODO: provide command line control of lambda and Delta lambda/lambda
        #L, dLoL = 5, 0.14/np.sqrt(6)  # wavelength and 14% triangular FWHM
        L, dLoL = 0, 0
        data = empty_data1D(q, resolution=res, L=L, dL=L*dLoL)
        index = slice(None, None)
    return data, index

def make_engine(
    model_info: ModelInfo, data: Data, dtype: str, cutoff: float, ngauss: int=0,
    ) -> Calculator:
    # not type: (ModelInfo, Data, str, float, int) -> Calculator
    """
    Generate the appropriate calculation engine for the given datatype.

    Datatypes with '!' appended are evaluated using external C DLLs rather
    than OpenCL.
    """
    if ngauss:
        set_integration_size(model_info, ngauss)

    if (dtype != "default" and not dtype.endswith('!')
            and not (kernelcl.use_opencl() or kernelcuda.use_cuda())):
        raise RuntimeError("OpenCL not available " + kernelcl.OPENCL_ERROR)

    model = core.build_model(model_info, dtype=dtype, platform="ocl")
    calculator = DirectModel(data, model, cutoff=cutoff)
    engine_type = calculator._model.__class__.__name__.replace('Model', '').upper()
    bits = calculator._model.dtype.itemsize*8
    precision = "fast" if getattr(calculator._model, 'fast', False) else str(bits)
    calculator.engine = "%s[%s]" % (engine_type, precision)
    return calculator

def _show_invalid(data, theory):
    # type: (Data, np.ma.ndarray) -> None
    """
    Display a list of the non-finite values in theory.
    """
    if not theory.mask.any():
        return

    if hasattr(data, 'x'):
        bad = zip(data.x[theory.mask], theory[theory.mask])
        print("   *** ", ", ".join("I(%g)=%g"%(x, y) for x, y in bad))


def compare(opts, limits=None, maxdim=None):
    # type: (Dict[str, Any], Optional[Tuple[float, float]], Optional[float]) -> Tuple[float, float]
    """
    Preform a comparison using options from the command line.

    *limits* are the display limits on the graph, either to set the y-axis
    for 1D or to set the colormap scale for 2D.  If None, then they are
    inferred from the data and returned. When exploring using Bumps,
    the limits are set when the model is initially called, and maintained
    as the values are adjusted, making it easier to see the effects of the
    parameters.

    *maxdim* **DEPRECATED** Use opts['maxdim'] instead.
    """
    # CRUFT: remove maxdim parameter
    if maxdim is not None:
        opts['maxdim'] = maxdim
    for k in range(opts['sets']):
        if k > 0:
            # print a separate seed for each dataset for better reproducibility
            new_seed = np.random.randint(1000000)  # type: int
            print("=== Set %d uses -random=%d ===" % (k+1, new_seed))
            np.random.seed(new_seed)
        opts['pars'] = parse_pars(opts, maxdim=maxdim)
        if opts['pars'] is None:
            return limits
        result = run_models(opts, verbose=True)
        if opts['plot']:
            if opts['is2d'] and k > 0:
                import matplotlib.pyplot as plt
                plt.figure()
            limits = plot_models(opts, result, limits=limits, setnum=k)
        if opts['show_weights']:
            base, _ = opts['engines']
            base_pars, _ = opts['pars']
            model_info = base._kernel.info
            dim = base._kernel.dim
            weights.plot_weights(model_info, get_mesh(model_info, base_pars, dim=dim))
        if opts['show_profile']:
            import pylab
            base, comp = opts['engines']
            base_pars, comp_pars = opts['pars']
            have_base = base._kernel.info.profile is not None
            have_comp = (
                comp is not None
                and comp._kernel.info.profile is not None
                and base_pars != comp_pars
            )
            if have_base or have_comp:
                pylab.figure()
            if have_base:
                plot_profile(base._kernel.info, **base_pars)
            if have_comp:
                plot_profile(comp._kernel.info, label='comp', **comp_pars)
                pylab.legend()
    if opts['plot']:
        import matplotlib.pyplot as plt
        plt.show()
    return limits

def plot_profile(model_info, label='base', **args):
    # type: (ModelInfo, List[Tuple[float, np.ndarray, np.ndarray]], float) -> None
    """
    Plot the profile returned by the model profile method.

    *model_info* defines model parameters, etc.

    *label* is the legend label for the plotted line.

    *args* are *parameter=value* pairs for the model profile function.
    """
    import pylab

    args = dict((k, v) for k, v in args.items()
                if "_pd" not in k
                and ":" not in k
                and k not in ("background", "scale", "theta", "phi", "psi"))
    args = args.copy()

    args.pop('scale', 1.)
    args.pop('background', 0.)
    z, rho = model_info.profile(**args)
    #pylab.interactive(True)
    pylab.plot(z, rho, '-', label=label)
    pylab.grid(True)
    #pylab.show()



def run_models(opts, verbose=False):
    # type: (Dict[str, Any], bool) -> Dict[str, Any]
    """
    Process a parameter set, return calculation results and times.
    """

    base, comp = opts['engines']
    base_n, comp_n = opts['count']
    base_pars, comp_pars = opts['pars']
    base_data, comp_data = opts['data']

    comparison = comp is not None

    base_time = comp_time = None
    base_value = comp_value = resid = relerr = None

    # Base calculation
    try:
        base_raw, base_time = time_calculation(base, base_pars, base_n)
        base_value = np.ma.masked_invalid(base_raw)
        if verbose:
            print("%s t=%.2f ms, intensity=%.0f"
                  % (base.engine, base_time, base_value.sum()))
        _show_invalid(base_data, base_value)
        #if base.results is not None: print(base.results())
    except ImportError:
        traceback.print_exc()

    # Comparison calculation
    if comparison:
        try:
            comp_raw, comp_time = time_calculation(comp, comp_pars, comp_n)
            comp_value = np.ma.masked_invalid(comp_raw)
            if verbose:
                print("%s t=%.2f ms, intensity=%.0f"
                      % (comp.engine, comp_time, comp_value.sum()))
            _show_invalid(base_data, comp_value)
        except ImportError:
            traceback.print_exc()

    # Compare, but only if computing both forms
    if comparison:
        resid = (base_value - comp_value)
        relerr = resid/np.where(comp_value != 0., abs(comp_value), 1.0)
        if verbose:
            _print_stats("|%s-%s|"
                         % (base.engine, comp.engine) + (" "*(3+len(comp.engine))),
                         resid)
            _print_stats("|(%s-%s)/%s|"
                         % (base.engine, comp.engine, comp.engine),
                         relerr)

    return dict(base_value=base_value, comp_value=comp_value,
                base_time=base_time, comp_time=comp_time,
                resid=resid, relerr=relerr)


def _print_stats(label, err):
    # type: (str, np.ma.ndarray) -> None
    # work with trimmed data, not the full set
    sorted_err = np.sort(abs(err.compressed()))
    if sorted_err.size == 0:
        print(label + "  no valid values")
        return

    p50 = int((len(sorted_err)-1)*0.50)
    p98 = int((len(sorted_err)-1)*0.98)
    data = [
        "max:%.3e"%sorted_err[-1],
        "median:%.3e"%sorted_err[p50],
        "98%%:%.3e"%sorted_err[p98],
        "rms:%.3e"%np.sqrt(np.mean(sorted_err**2)),
        "zero-offset:%+.3e"%np.mean(sorted_err),
        ]
    print(label+"  "+"  ".join(data))


def plot_models(opts, result, limits=None, setnum=0):
    # type: (Dict[str, Any], Dict[str, Any], Optional[Tuple[float, float]], int) -> Tuple[float, float]
    """
    Plot the results from :func:`run_models`.
    """
    import matplotlib.pyplot as plt

    base_value, comp_value = result['base_value'], result['comp_value']
    base_time, comp_time = result['base_time'], result['comp_time']
    resid, relerr = result['resid'], result['relerr']

    have_base, have_comp = (base_value is not None), (comp_value is not None)
    base, comp = opts['engines']
    base_data, comp_data = opts['data']
    use_data = (opts['datafile'] is not None) and (have_base ^ have_comp)

    # Plot if requested
    view = opts['view']
    #view = 'log'
    if limits is None:
        vmin, vmax = np.inf, -np.inf
        if have_base:
            vmin = min(vmin, base_value.min())
            vmax = max(vmax, base_value.max())
        if have_comp:
            vmin = min(vmin, comp_value.min())
            vmax = max(vmax, comp_value.max())
        limits = vmin, vmax

    if have_base:
        if have_comp:
            plt.subplot(131)
        plot_theory(base_data, base_value, view=view, use_data=use_data, limits=limits)
        plt.title("%s t=%.2f ms"%(base.engine, base_time))
        #cbar_title = "log I"
    if have_comp:
        if have_base:
            plt.subplot(132)
        if not opts['is2d'] and have_base:
            plot_theory(comp_data, base_value, view=view, use_data=use_data, limits=limits)
        plot_theory(comp_data, comp_value, view=view, use_data=use_data, limits=limits)
        plt.title("%s t=%.2f ms"%(comp.engine, comp_time))
        #cbar_title = "log I"
    if have_base and have_comp:
        plt.subplot(133)
        if not opts['rel_err']:
            err, errstr, errview = resid, "abs err", "linear"
        else:
            err, errstr, errview = abs(relerr), "rel err", "log"
            if (err == 0.).all():
                errview = 'linear'
        if 0:  # 95% cutoff
            sorted_err = np.sort(err.flatten())
            cutoff = sorted_err[int(sorted_err.size*0.95)]
            err[err > cutoff] = cutoff
        #err,errstr = base/comp,"ratio"
        # Note: base_data only since base and comp have same q values (though
        # perhaps different resolution), and we are plotting the difference
        # at each q
        plot_theory(base_data, None, resid=err, view=errview, use_data=use_data)
        plt.xscale('log' if view == 'log' and not opts['is2d'] else 'linear')
        plt.legend(['P%d'%(k+1) for k in range(setnum+1)], loc='best')
        plt.title("max %s = %.3g"%(errstr, abs(err).max()))
        #cbar_title = errstr if errview=="linear" else "log "+errstr
    #if is2D:
    #    h = plt.colorbar()
    #    h.ax.set_title(cbar_title)
    fig = plt.gcf()
    extra_title = ' '+opts['title'] if opts['title'] else ''
    fig.suptitle(":".join(opts['name']) + extra_title)

    if have_base and have_comp and opts['show_hist']:
        plt.figure()
        v = relerr
        v[v == 0] = 0.5*np.min(np.abs(v[v != 0]))
        plt.hist(np.log10(np.abs(v)), normed=1, bins=50)
        plt.xlabel('log10(err), err = |(%s - %s) / %s|'
                   % (base.engine, comp.engine, comp.engine))
        plt.ylabel('P(err)')
        plt.title('Distribution of relative error between calculation engines')

    return limits


# ===========================================================================
#

# Set of command line options.
# Normal options such as -plot/-noplot are specified as 'name'.
# For options such as -nq=500 which require a value use 'name='.
#
OPTIONS = [
    # Plotting
    'plot', 'noplot',
    'weights', 'profile',
    'linear', 'log', 'q4',
    'rel', 'abs',
    'hist', 'nohist',
    'title=',

    # Data generation
    'data=', 'noise=', 'res=', 'nq=', 'q=',
    'lowq', 'midq', 'highq', 'exq', 'zero',
    '2d', '1d',

    # Parameter set
    'preset', 'random', 'random=', 'sets=',
    'nopars', 'pars',
    'sphere', 'sphere=', # integrate over a sphere in 2d with n points
    'poly', 'mono',
    'magnetic', 'nonmagnetic',
    'maxdim=',

    # Calculation options
    'cutoff=', 'accuracy=', 'ngauss=',
    'neval=',  # for timing...

    # Precision options
    'engine=',
    'half', 'fast', 'single', 'double', 'single!', 'double!', 'quad!',

    # Output options
    'help', 'html', 'edit',

    # Help options
    'h', '?', 'models', 'models='
    ]

NAME_OPTIONS = (lambda: set(k for k in OPTIONS if not k.endswith('=')))()
VALUE_OPTIONS = (lambda: [k[:-1] for k in OPTIONS if k.endswith('=')])()


def columnize(items, indent="", width=None):
    # type: (List[str], str, int) -> str
    """
    Format a list of strings into columns.

    Returns a string with carriage returns ready for printing.
    """
    # Use the columnize package (pycolumize) if it is available
    try:
        from columnize import columnize as _columnize, default_opts
        if width is None:
            width = default_opts['displaywidth']
        return _columnize(list(items), displaywidth=width, lineprefix=indent)
    except ImportError:
        pass
    # Otherwise roll our own.
    if width is None:
        width = 120
    column_width = max(len(w) for w in items) + 1
    num_columns = (width - len(indent)) // column_width
    num_rows = len(items) // num_columns
    items = items + [""] * (num_rows * num_columns - len(items))
    columns = [items[k*num_rows:(k+1)*num_rows] for k in range(num_columns)]
    lines = [" ".join("%-*s"%(column_width, entry) for entry in row)
             for row in zip(*columns)]
    output = indent + ("\n"+indent).join(lines)
    return output


def get_pars(model_info):
    # type: (ModelInfo) -> ParameterSet
    """
    Extract default parameters from the model definition.
    """
    # Get the default values for the parameters
    pars = {}
    for p in model_info.parameters.call_parameters:
        parts = [('', p.default)]
        if p.polydisperse:
            parts.append(('_pd', 0.0))
            parts.append(('_pd_n', 0))
            parts.append(('_pd_nsigma', 3.0))
            parts.append(('_pd_type', "gaussian"))
        for ext, val in parts:
            if p.length > 1:
                dict(("%s%d%s" % (p.id, k, ext), val)
                     for k in range(1, p.length+1))
            else:
                pars[p.id + ext] = val
    return pars

INTEGER_RE = re.compile("^[+-]?[1-9][0-9]*$")
def isnumber(s):
    # type: (str) -> bool
    """Return True if string contains an int or float"""
    match = FLOAT_RE.match(s)
    isfloat = (match and not s[match.end():])
    return isfloat or INTEGER_RE.match(s)

def print_models(kind=None):
    """
    Print the list of available models in columns.
    """
    models = core.list_models(kind=kind)
    print(columnize(models, indent="  "))

# For distinguishing pairs of models for comparison
# key-value pair separator =
# shell characters  | & ; <> $ % ' " \ # `
# model and parameter names _
# parameter expressions - + * / . ( )
# path characters including tilde expansion and windows drive ~ / :
# not sure about brackets [] {}
# maybe one of the following @ ? ^ ! ,
PAR_SPLIT = ','
def parse_opts(argv):
    # type: (List[str]) -> Dict[str, Any]
    """
    Parse command line options.
    """

    flags = [arg for arg in argv
             if arg.startswith('-')]
    values = [arg for arg in argv
              if not arg.startswith('-') and '=' in arg]
    positional_args = [arg for arg in argv
                       if not arg.startswith('-') and '=' not in arg]

    # First check if help requested anywhere on line
    if '-h' in flags or '-?' in flags:
        print(USAGE)
        return None

    # Next check that all flags are valid.
    invalid = [o[1:] for o in flags
               if not (o[1:] in NAME_OPTIONS
                       or any(o.startswith('-%s='%t) for t in VALUE_OPTIONS)
                       or o.startswith('-D'))]
    if invalid:
        print("Invalid options: %s."%(", ".join(invalid)))
        print("usage: ./sasmodels [-?] [-models] model")
        return None

    # Check if requesting a list of models.  This is done after checking that
    # the flags are valid so we know it is -models or -models=.
    if any(v.startswith('-models') for v in flags):
        # grab last -models entry
        models = [v for v in flags if v.startswith('-models')][-1]
        if models == '-models':
            models = '-models=all'
        _, kind = models.split('=', 1)
        print_models(kind=kind)
        return None

    # Check that a model was given on the command line
    if not positional_args:
        print("usage: ./sascomp [-?] [-models] model")
        return None

    # Only the last model on the command line is used.
    name = positional_args[-1]

    # Interpret the flags
    # pylint: disable=bad-whitespace,C0321
    opts = {
        'plot'      : True,
        'view'      : 'log',
        'is2d'      : False,
        'qmin'      : None,
        'qmax'      : 0.05,
        'nq'        : 128,
        'res'       : '0.0',
        'noise'     : 0.0,
        'accuracy'  : 'Low',
        'cutoff'    : '0.0',
        'seed'      : -1,  # default to preset
        'mono'      : True,
        # Default to magnetic a magnetic moment is set on the command line
        'magnetic'  : False,
        'maxdim'    : np.inf,
        'show_pars' : False,
        'show_hist' : False,
        'rel_err'   : True,
        'explore'   : False,
        'zero'      : False,
        'html'      : False,
        'title'     : None,
        'datafile'  : None,
        'sets'      : 0,
        'engine'    : 'default',
        'count'     : '1',
        'show_weights' : False,
        'show_profile' : False,
        'sphere'    : 0,
        'ngauss'    : '0',
    }
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
        elif arg == '-zero':    opts['zero'] = True
        elif arg.startswith('-nq='):       opts['nq'] = int(arg[4:])
        elif arg.startswith('-q='):
            opts['qmin'], opts['qmax'] = [float(v) for v in arg[3:].split(':')]
        elif arg.startswith('-res='):      opts['res'] = arg[5:]
        elif arg.startswith('-noise='):    opts['noise'] = float(arg[7:])
        elif arg.startswith('-sets='):     opts['sets'] = int(arg[6:])
        elif arg.startswith('-accuracy='): opts['accuracy'] = arg[10:]
        elif arg.startswith('-cutoff='):   opts['cutoff'] = arg[8:]
        elif arg.startswith('-title='):    opts['title'] = arg[7:]
        elif arg.startswith('-data='):     opts['datafile'] = arg[6:]
        elif arg.startswith('-engine='):   opts['engine'] = arg[8:]
        elif arg.startswith('-neval='):    opts['count'] = arg[7:]
        elif arg.startswith('-ngauss='):   opts['ngauss'] = arg[8:]
        elif arg.startswith('-random='):
            opts['seed'] = int(arg[8:])
            opts['sets'] = 0
        elif arg == '-random':
            opts['seed'] = np.random.randint(1000000)
            opts['sets'] = 0
        elif arg.startswith('-sphere'):
            opts['sphere'] = int(arg[8:]) if len(arg) > 7 else 150
            opts['is2d'] = True
        elif arg.startswith('-maxdim'):
            opts['maxdim'] = float(arg[8:])
        elif arg == '-preset':  opts['seed'] = -1
        elif arg == '-mono':    opts['mono'] = True
        elif arg == '-poly':    opts['mono'] = False
        elif arg == '-magnetic':       opts['magnetic'] = True
        elif arg == '-nonmagnetic':    opts['magnetic'] = False
        elif arg == '-pars':    opts['show_pars'] = True
        elif arg == '-nopars':  opts['show_pars'] = False
        elif arg == '-hist':    opts['show_hist'] = True
        elif arg == '-nohist':  opts['show_hist'] = False
        elif arg == '-rel':     opts['rel_err'] = True
        elif arg == '-abs':     opts['rel_err'] = False
        elif arg == '-half':    opts['engine'] = 'half'
        elif arg == '-fast':    opts['engine'] = 'fast'
        elif arg == '-single':  opts['engine'] = 'single'
        elif arg == '-double':  opts['engine'] = 'double'
        elif arg == '-single!': opts['engine'] = 'single!'
        elif arg == '-double!': opts['engine'] = 'double!'
        elif arg == '-quad!':   opts['engine'] = 'quad!'
        elif arg == '-edit':    opts['explore'] = True
        elif arg == '-weights': opts['show_weights'] = True
        elif arg == '-profile': opts['show_profile'] = True
        elif arg == '-html':    opts['html'] = True
        elif arg == '-help':    opts['html'] = True
        elif arg.startswith('-D'):
            var, val = arg[2:].split('=')
            os.environ[var] = val
    # pylint: enable=bad-whitespace,C0321

    # Magnetism forces 2D for now
    if opts['magnetic']:
        opts['is2d'] = True

    # Force random if sets is used
    if opts['sets'] >= 1 and opts['seed'] < 0:
        opts['seed'] = np.random.randint(1000000)
    if opts['sets'] == 0:
        opts['sets'] = 1

    # Create the computational engines
    if opts['qmin'] is None:
        opts['qmin'] = 0.001*opts['qmax']

    comparison = any(PAR_SPLIT in v for v in values)

    if PAR_SPLIT in name:
        names = name.split(PAR_SPLIT, 2)
        comparison = True
    else:
        names = [name]*2
    try:
        model_info = [core.load_model_info(k) for k in names]
    except ImportError as exc:
        print(str(exc), "while loading", names)
        print("usage: ./sasmodels [-?] [-models] model")
        return None

    if PAR_SPLIT in opts['ngauss']:
        opts['ngauss'] = [int(k) for k in opts['ngauss'].split(PAR_SPLIT, 2)]
        comparison = True
    else:
        opts['ngauss'] = [int(opts['ngauss'])]*2

    if PAR_SPLIT in opts['engine']:
        opts['engine'] = opts['engine'].split(PAR_SPLIT, 2)
        comparison = True
    else:
        opts['engine'] = [opts['engine']]*2

    if PAR_SPLIT in opts['count']:
        opts['count'] = [int(k) for k in opts['count'].split(PAR_SPLIT, 2)]
        comparison = True
    else:
        opts['count'] = [int(opts['count'])]*2

    if PAR_SPLIT in opts['cutoff']:
        opts['cutoff'] = [float(k) for k in opts['cutoff'].split(PAR_SPLIT, 2)]
        comparison = True
    else:
        opts['cutoff'] = [float(opts['cutoff'])]*2

    if PAR_SPLIT in opts['res']:
        opts['res'] = [float(k) for k in opts['res'].split(PAR_SPLIT, 2)]
        comparison = True
    else:
        opts['res'] = [float(opts['res'])]*2

    if opts['datafile'] is not None:
        data0 = load_data(os.path.expanduser(opts['datafile']))
        data = data0, data0
    else:
        # Hack around the fact that make_data doesn't take a pair of resolutions
        res = opts['res']
        opts['res'] = res[0]
        data0, _ = make_data(opts)
        if res[0] != res[1]:
            opts['res'] = res[1]
            data1, _ = make_data(opts)
        else:
            data1 = data0
        opts['res'] = res
        data = data0, data1

    base = make_engine(model_info[0], data[0], opts['engine'][0],
                       opts['cutoff'][0], opts['ngauss'][0])
    if comparison:
        comp = make_engine(model_info[1], data[1], opts['engine'][1],
                           opts['cutoff'][1], opts['ngauss'][1])
    else:
        comp = None

    # pylint: disable=bad-whitespace
    # Remember it all
    opts.update({
        'data'      : data,
        'name'      : names,
        'info'      : model_info,
        'engines'   : [base, comp],
        'values'    : values,
    })
    # pylint: enable=bad-whitespace

    # Set the integration parameters to the half sphere
    if opts['sphere'] > 0:
        set_spherical_integration_parameters(opts, opts['sphere'])

    return opts

def set_spherical_integration_parameters(opts, steps):
    # type: (Dict[str, Any], int) -> None
    """
    Set integration parameters for spherical integration over the entire
    surface in theta-phi coordinates.
    """
    # Set the integration parameters to the half sphere
    opts['values'].extend([
        #'theta=90',
        'theta_pd=%g'%(90/np.sqrt(3)),
        'theta_pd_n=%d'%steps,
        'theta_pd_type=rectangle',
        #'phi=0',
        'phi_pd=%g'%(180/np.sqrt(3)),
        'phi_pd_n=%d'%(2*steps),
        'phi_pd_type=rectangle',
        #'background=0',
    ])
    if 'psi' in opts['info'][0].parameters:
        opts['values'].extend([
            #'psi=0',
            'psi_pd=%g'%(180/np.sqrt(3)),
            'psi_pd_n=%d'%(2*steps),
            'psi_pd_type=rectangle',
        ])

def parse_pars(opts, maxdim=None):
    # type: (Dict[str, Any], float) -> Tuple[Dict[str, float], Dict[str, float]]
    """
    Generate parameter sets for base and comparison models.

    Returns a pair of parameter dictionaries.

    The default parameter values come from the model, or a randomized model
    if a seed value is given.  Next, evaluate any parameter expressions,
    constraining the value of the parameter within and between models.

    Note: When generating random parameters, **the seed must already be set**
    with a call to *np.random.seed(opts['seed'])*.

    *opts* controls the parameter generation::

        opts = {
            'info': (model_info 1, model_info 2),
            'seed': -1,         # if seed>=0 then randomize parameters
            'mono': False,      # force monodisperse random parameters
            'magnetic': False,  # force nonmagetic random parameters
            'maxdim': np.inf,   # limit particle size to maxdim for random pars
            'values': ['par=expr', ...],  # override parameter values in model
            'show_pars': False, # Show parameter values
            'is2d': False,      # Show values for orientation parameters
        }

    The values of *par=expr* are evaluated approximately as::

        import numpy as np
        from math import *
        from parameter_set import *

        parameter_set.par = eval(expr)

    That is, you can use arbitrary python math expressions including the
    functions defined in the math library and the numpy library.  You can
    also use the existing parameter values, which will either be the model
    defaults or the randomly generated values if seed is non-negative.

    To compare different values of the same parameter, use *par=expr,expr*.
    The first parameter set will have the values from the first expression
    and the second parameter set will have the values from the second
    expression.  Note that the second expression is evaluated using the
    values from the first expression, which allows things like::

        length=2*radius,length+3

    which will compare length to length+3 when length is set to 2*radius.

    *maxdim* **DEPRECATED** Use *opts['maxdim']* instead.
    """
    # CRUFT: maxdim parameter is deprecated
    if maxdim is not None:
        opts['maxdim'] = maxdim

    model_info, model_info2 = opts['info']

    # Get default parameters from model definition.
    pars = get_pars(model_info)
    pars2 = get_pars(model_info2)
    pars2.update((k, v) for k, v in pars.items() if k in pars2)
    # randomize parameters
    #pars.update(set_pars)  # set value before random to control range
    if opts['seed'] > -1:
        pars = randomize_pars(model_info, pars, maxdim=opts['maxdim'])
        if model_info.id != model_info2.id:
            pars2 = randomize_pars(model_info2, pars2, maxdim=opts['maxdim'])
            # Share values for parameters with the same name
            for k, v in pars.items():
                if k in pars2:
                    pars2[k] = v
        else:
            pars2 = pars.copy()
        constrain_pars(model_info, pars)
        constrain_pars(model_info2, pars2)
        # TODO: randomly contrast match a pair of SLDs with some probability

    # Process -mono and -magnetic command line options.
    if opts['mono']:
        pars = suppress_pd(pars)
        pars2 = suppress_pd(pars2)
    if not opts['magnetic']:
        pars = suppress_magnetism(pars)
        pars2 = suppress_magnetism(pars2)

    # Fill in parameters given on the command line
    presets = {}
    presets2 = {}
    for arg in opts['values']:
        k, v = arg.split('=', 1)
        if k not in pars and k not in pars2:
            # extract base name without polydispersity info
            s = set(p.split('_pd')[0] for p in pars)
            print("%r invalid; parameters are: %s"%(k, ", ".join(sorted(s))))
            return None
        v1, v2 = v.split(PAR_SPLIT, 2) if PAR_SPLIT in v else (v, v)
        if v1 and k in pars:
            presets[k] = float(v1) if isnumber(v1) else v1
        if v2 and k in pars2:
            presets2[k] = float(v2) if isnumber(v2) else v2

    # If pd given on the command line, default pd_n to 35
    for k, v in list(presets.items()):
        if k.endswith('_pd'):
            presets.setdefault(k+'_n', 35.)
    for k, v in list(presets2.items()):
        if k.endswith('_pd'):
            presets2.setdefault(k+'_n', 35.)

    # Evaluate preset parameter expressions
    # Note: need to replace ':' with '_' in parameter names and expressions
    # in order to support math on magnetic parameters.
    context = MATH.copy()
    context['np'] = np
    context.update((k.replace(':', '_'), v) for k, v in pars.items())
    context.update((k, v) for k, v in presets.items() if isinstance(v, float))
    #for k,v in sorted(context.items()): print(k, v)
    for k, v in presets.items():
        if not isinstance(v, float) and not k.endswith('_type'):
            presets[k] = eval(v.replace(':', '_'), context)
    context.update(presets)
    context.update((k.replace(':', '_'), v) for k, v in presets2.items() if isinstance(v, float))
    for k, v in presets2.items():
        if not isinstance(v, float) and not k.endswith('_type'):
            presets2[k] = eval(v.replace(':', '_'), context)

    # update parameters with presets
    pars.update(presets)  # set value after random to control value
    pars2.update(presets2)  # set value after random to control value
    #import pprint; pprint.pprint(model_info)

    # Hack to load user-defined distributions; run through all parameters
    # and make sure any pd_type parameter is a defined distribution.
    if (any(p.endswith('pd_type') and v not in weights.DISTRIBUTIONS
            for p, v in pars.items())
            or any(p.endswith('pd_type') and v not in weights.DISTRIBUTIONS
                   for p, v in pars2.items())):
        weights.load_weights()

    if opts['show_pars']:
        if model_info.name != model_info2.name or pars != pars2:
            print("==== %s ====="%model_info.name)
            print(str(parlist(model_info, pars, opts['is2d'])))
            print("==== %s ====="%model_info2.name)
            print(str(parlist(model_info2, pars2, opts['is2d'])))
        else:
            print(str(parlist(model_info, pars, opts['is2d'])))

    return pars, pars2

def show_docs(opts):
    # type: (Dict[str, Any]) -> None
    """
    show html docs for the model
    """
    from .generate import make_html
    from . import rst2html

    info = opts['info'][0]
    html = make_html(info)
    path = os.path.dirname(info.filename)
    url = "file://" + path.replace("\\", "/")[2:] + "/"
    rst2html.view_html_wxapp(html, url)

def explore(opts):
    # type: (Dict[str, Any]) -> None
    """
    explore the model using the bumps gui.
    """
    import wx  # type: ignore
    from bumps.names import FitProblem  # type: ignore
    from bumps.gui.app_frame import AppFrame  # type: ignore
    from bumps.gui import signal

    is_mac = "cocoa" in wx.version()
    # Create an app if not running embedded
    app = wx.App() if wx.GetApp() is None else None
    model = Explore(opts)
    problem = FitProblem(model)
    frame = AppFrame(parent=None, title="explore", size=(1000, 700))
    if not is_mac:
        frame.Show()
    frame.panel.set_model(model=problem)
    frame.panel.Layout()
    frame.panel.aui.Split(0, wx.TOP)
    def _reset_parameters(event):
        model.revert_values()
        signal.update_parameters(problem)
    frame.Bind(wx.EVT_TOOL, _reset_parameters, frame.ToolBar.GetToolByPos(1))
    if is_mac:
        frame.Show()
    # If running withing an app, start the main loop
    if app:
        app.MainLoop()

class Explore(object):
    """
    Bumps wrapper for a SAS model comparison.

    The resulting object can be used as a Bumps fit problem so that
    parameters can be adjusted in the GUI, with plots updated on the fly.
    """
    def __init__(self, opts):
        # type: (Dict[str, Any]) -> None
        from bumps.cli import config_matplotlib  # type: ignore
        from . import bumps_model
        config_matplotlib()
        self.opts = opts
        opts['pars'] = list(opts['pars'])
        p1, p2 = opts['pars']
        m1, m2 = opts['info']
        self.fix_p2 = m1 != m2 or p1 != p2
        model_info = m1
        pars, pd_types = bumps_model.create_parameters(model_info, **p1)
        # Initialize parameter ranges, fixing the 2D parameters for 1D data.
        if not opts['is2d']:
            for p in model_info.parameters.user_parameters({}, is2d=False):
                for ext in ['', '_pd', '_pd_n', '_pd_nsigma']:
                    k = p.name+ext
                    v = pars.get(k, None)
                    if v is not None:
                        v.range(*parameter_range(k, v.value))
        else:
            for k, v in pars.items():
                v.range(*parameter_range(k, v.value))

        self.pars = pars
        self.starting_values = dict((k, v.value) for k, v in pars.items())
        self.pd_types = pd_types
        self.limits = None

    def revert_values(self):
        # type: () -> None
        """
        Restore starting values of the parameters.
        """
        for k, v in self.starting_values.items():
            self.pars[k].value = v

    def model_update(self):
        # type: () -> None
        """
        Respond to signal that model parameters have been changed.
        """
        pass

    def numpoints(self):
        # type: () -> int
        """
        Return the number of points.
        """
        return len(self.pars) + 1  # so dof is 1

    def parameters(self):
        # type: () -> Any   # Dict/List hierarchy of parameters
        """
        Return a dictionary of parameters.
        """
        return self.pars

    def nllf(self):
        # type: () -> float
        """
        Return cost.
        """
        # pylint: disable=no-self-use
        return 0.  # No nllf

    def plot(self, view='log'):
        # type: (str) -> None
        """
        Plot the data and residuals.
        """
        pars = dict((k, v.value) for k, v in self.pars.items())
        pars.update(self.pd_types)
        self.opts['pars'][0] = pars
        if not self.fix_p2:
            self.opts['pars'][1] = pars
        result = run_models(self.opts)
        limits = plot_models(self.opts, result, limits=self.limits)
        if self.limits is None:
            vmin, vmax = limits
            self.limits = vmax*1e-7, 1.3*vmax
            import pylab
            pylab.clf()
            plot_models(self.opts, result, limits=self.limits)


def main(*argv):
    # type: (*str) -> None
    """
    Main program.
    """
    opts = parse_opts(argv)
    if opts is not None:
        if opts['seed'] > -1:
            print("Randomize using -random=%i"%opts['seed'])
            np.random.seed(opts['seed'])
        if opts['html']:
            show_docs(opts)
        elif opts['explore']:
            opts['pars'] = parse_pars(opts)
            if opts['pars'] is None:
                return
            explore(opts)
        else:
            compare(opts)

if __name__ == "__main__":
    main(*sys.argv[1:])
