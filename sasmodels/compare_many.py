#!/usr/bin/env python
import sys
import traceback

import numpy as np

from . import core
from .kernelcl import environment
from .compare import (MODELS, randomize_model, suppress_pd, eval_sasview,
                      eval_opencl, eval_ctypes, make_data, get_demo_pars,
                      columnize, constrain_pars)

def calc_stats(target, value, index):
    resid = abs(value-target)[index]
    relerr = resid/target[index]
    srel = np.argsort(relerr)
    #p90 = int(len(relerr)*0.90)
    p95 = int(len(relerr)*0.95)
    maxrel = np.max(relerr)
    rel95 = relerr[srel[p95]]
    maxabs = np.max(resid[srel[p95:]])
    maxval = np.max(value[srel[p95:]])
    return maxrel,rel95,maxabs,maxval

def print_column_headers(pars, parts):
    stats = list('Max rel err|95% rel err|Max abs err above 90% rel|Max value above 90% rel'.split('|'))
    groups = ['']
    for p in parts:
        groups.append(p)
        groups.extend(['']*(len(stats)-1))
    groups.append("Parameters")
    columns = ['Seed'] + stats*len(parts) +  list(sorted(pars.keys()))
    print(','.join('"%s"'%c for c in groups))
    print(','.join('"%s"'%c for c in columns))

def compare_instance(name, data, index, N=1, mono=True, cutoff=1e-5,
                     precision='double'):
    model_definition = core.load_model_definition(name)
    pars = get_demo_pars(model_definition)
    header = '\n"Model","%s","Count","%d"'%(name, N)
    if not mono: header += ',"Cutoff",%g'%(cutoff,)
    print(header)

    # Some not very clean macros for evaluating the models and checking the
    # results.  They freely use variables from the current scope, even some
    # which have not been defined yet, complete with abuse of mutable lists
    # to allow them to update values in the current scope since nonlocal
    # declarations are not available in python 2.7.
    def try_model(fn, *args, **kw):
        try:
            result, _ = fn(model_definition, pars_i, data, *args, **kw)
        except KeyboardInterrupt:
            raise
        except:
            traceback.print_exc()
            print("when comparing %s for %d"%(name, seed))
            if hasattr(data, 'qx_data'):
                result = np.NaN*data.data
            else:
                result = np.NaN*data.x
        return result
    def check_model(label, target, value, acceptable):
        stats = calc_stats(target, value, index)
        columns.extend(stats)
        labels.append('GPU single')
        max_diff[0] = max(max_diff[0], stats[0])
        good[0] = good[0] and (stats[0] < acceptable)

    num_good = 0
    first = True
    max_diff = [0]
    for k in range(N):
        print("%s %d"%(name, k))
        pars_i, seed = randomize_model(pars)
        constrain_pars(model_definition, pars_i)
        if mono: suppress_pd(pars_i)

        good = [True]
        labels = []
        columns = []
        target = try_model(eval_sasview)
        #target = try_model(eval_ctypes, dtype='double', cutoff=0.)
        #target = try_model(eval_ctypes, dtype='longdouble', cutoff=0.)
        if precision == 'single':
            value = try_model(eval_opencl, dtype='single', cutoff=cutoff)
            check_model('GPU single', target, value, 5e-5)
            single_value = value  # remember for single/double comparison
        elif precision == 'double':
            if environment().has_type('double'):
                label = 'GPU double'
                value = try_model(eval_opencl, dtype='double', cutoff=cutoff)
            else:
                label = 'CPU double'
                value = try_model(eval_ctypes, dtype='double', cutoff=cutoff)
            check_model(label, target, value, 5e-14)
            double_value = value  # remember for single/double comparison
        elif precision == 'quad':
            value = try_model(eval_opencl, dtype='longdouble', cutoff=cutoff)
            check_model('CPU quad', target, value, 5e-14)
        if 0:
            check_model('single/double', double_value, single_value, 5e-5)

        columns += [v for _,v in sorted(pars_i.items())]
        if first:
            print_column_headers(pars_i, labels)
            first = False
        if good[0]:
            num_good += 1
        else:
            print(("%d,"%seed)+','.join("%g"%v for v in columns))
    print('"good","%d/%d","max diff",%g'%(num_good, N, max_diff[0]))


def print_usage():
    print("usage: compare_many.py MODEL COUNT (1dNQ|2dNQ) (CUTOFF|mono) (single|double|quad)")


def print_models():
    print(columnize(MODELS, indent="  "))


def print_help():
    print_usage()
    print("""\

MODEL is the model name of the model or "all" for all the models
in alphabetical order.

COUNT is the number of randomly generated parameter sets to try. A value
of "10000" is a reasonable check for monodisperse models, or "100" for
polydisperse models.   For a quick check, use "100" and "5" respectively.

NQ is the number of Q values to calculate.  If it starts with "1d", then
it is a 1-dimensional problem, with log spaced Q points from 1e-3 to 1.0.
If it starts with "2d" then it is a 2-dimensional problem, with linearly
spaced points Q points from -1.0 to 1.0 in each dimension. The usual
values are "1d100" for 1-D and "2d32" for 2-D.

CUTOFF is the cutoff value to use for the polydisperse distribution. Weights
below the cutoff will be ignored.  Use "mono" for monodisperse models.  The
choice of polydisperse parameters, and the number of points in the distribution
is set in compare.py defaults for each model.

PRECISION is the floating point precision to use for comparisons.

Available models:
""")
    print_models()

def main():
    if len(sys.argv) != 6:
        print_help()
        sys.exit(1)

    model = sys.argv[1]
    if not (model in MODELS) and (model != "all"):
        print('Bad model %s.  Use "all" or one of:')
        print_models()
        sys.exit(1)
    try:
        count = int(sys.argv[2])
        is2D = sys.argv[3].startswith('2d')
        assert sys.argv[3][1] == 'd'
        Nq = int(sys.argv[3][2:])
        mono = sys.argv[4] == 'mono'
        cutoff = float(sys.argv[4]) if not mono else 0
        precision = sys.argv[5]
    except:
        traceback.print_exc()
        print_usage()
        sys.exit(1)

    data, index = make_data(qmax=1.0, is2D=is2D, Nq=Nq)
    model_list = [model] if model != "all" else MODELS
    for model in model_list:
        compare_instance(model, data, index, N=count, mono=mono,
                         cutoff=cutoff, precision=precision)

if __name__ == "__main__":
    main()
