#!/usr/bin/env python
import sys
import traceback

import numpy as np

from . import core
from .kernelcl import environment
from .compare import (MODELS, randomize_pars, suppress_pd, eval_sasview,
                      eval_opencl, eval_ctypes, make_data, get_demo_pars,
                      columnize, constrain_pars, constrain_new_to_old,
                      make_engine)

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

PRECISION = {
    'fast': 1e-3,
    'half': 1e-3,
    'single': 5e-5,
    'double': 5e-14,
    'single!': 5e-5,
    'double!': 5e-14,
    'quad!': 5e-18,
    'sasview': 5e-14,
}
def compare_instance(name, data, index, N=1, mono=True, cutoff=1e-5,
                     base='sasview', comp='double'):
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
    def try_model(fn, pars):
        try:
            result = fn(**pars)
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
    def check_model(pars):
        base_value = try_model(calc_base, pars)
        comp_value = try_model(calc_comp, pars)
        stats = calc_stats(base_value, comp_value, index)
        max_diff[0] = max(max_diff[0], stats[0])
        good[0] = good[0] and (stats[0] < expected)
        return list(stats)


    calc_base = make_engine(model_definition, data, base, cutoff)
    calc_comp = make_engine(model_definition, data, comp, cutoff)
    expected = max(PRECISION[base], PRECISION[comp])

    num_good = 0
    first = True
    max_diff = [0]
    for k in range(N):
        print("%s %d"%(name, k))
        seed = np.random.randint(1e6)
        pars_i = randomize_pars(pars, seed)
        constrain_pars(model_definition, pars_i)
        constrain_new_to_old(model_definition, pars_i)
        if mono:
            pars_i = suppress_pd(pars_i)

        good = [True]
        columns = check_model(pars_i)
        columns += [v for _,v in sorted(pars_i.items())]
        if first:
            labels = [" vs. ".join((calc_base.engine, calc_comp.engine))]
            print_column_headers(pars_i, labels)
            first = False
        if good[0]:
            num_good += 1
        else:
            print(("%d,"%seed)+','.join("%s"%v for v in columns))
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

PRECISION is the floating point precision to use for comparisons.  If two
precisions are given, then compare one to the other, ignoring sasview.

Available models:
""")
    print_models()

def main():
    if len(sys.argv) not in (6,7):
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
        base = sys.argv[5]
        comp = sys.argv[6] if len(sys.argv) > 6 else "sasview"
    except:
        traceback.print_exc()
        print_usage()
        sys.exit(1)

    data, index = make_data({'qmax':1.0, 'is2d':is2D, 'nq':Nq, 'res':0.,
                              'accuracy': 'Low', 'view':'log'})
    model_list = [model] if model != "all" else MODELS
    for model in model_list:
        compare_instance(model, data, index, N=count, mono=mono,
                         cutoff=cutoff, base=base, comp=comp)

if __name__ == "__main__":
    main()
