#!/usr/bin/env python
import sys
import traceback

import numpy as np

from sasmodels import core
from sasmodels.kernelcl import environment
from compare import (MODELS, randomize_model, suppress_pd, eval_sasview,
                     eval_opencl, eval_ctypes, make_data, get_demo_pars,
                     columnize)

def get_stats(target, value, index):
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

def compare_instance(name, data, index, N=1, mono=True, cutoff=1e-5):
    model_definition = core.load_model_definition(name)
    pars = get_demo_pars(name)
    header = '\n"Model","%s","Count","%d"'%(name, N)
    if not mono: header += ',"Cutoff",%g'%(cutoff,)
    print(header)

    # Stuff the failure flag into a mutable object so we can update it from
    # within the nested function.  Note that the nested function uses "pars"
    # which is dynamically scoped, not lexically scoped in this context.  That
    # is, pars is replaced each time in the loop, so don't assume that it is
    # the default values defined above.
    def trymodel(fn, *args, **kw):
        try:
            result, _ = fn(model_definition, pars, data, *args, **kw)
        except:
            result = np.NaN
            traceback.print_exc()
        return result

    num_good = 0
    first = True
    for _ in range(N):
        pars, seed = randomize_model(name, pars)
        if mono: suppress_pd(pars)

        # Force parameter constraints on a per-model basis.
        if name in ('teubner_strey','broad_peak'):
            pars['scale'] = 1.0
        #if name == 'parallelepiped':
        #    pars['a_side'],pars['b_side'],pars['c_side'] = \
        #        sorted([pars['a_side'],pars['b_side'],pars['c_side']])


        good = True
        labels = []
        columns = []
        if 1:
            sasview_value = trymodel(eval_sasview)
        if 0:
            gpu_single_value = trymodel(eval_opencl, dtype='single', cutoff=cutoff)
            stats = get_stats(sasview_value, gpu_single_value, index)
            columns.extend(stats)
            labels.append('GPU single')
            good = good and (stats[0] < 1e-14)
        if 0 and environment().has_double:
            gpu_double_value = trymodel(eval_opencl, dtype='double', cutoff=cutoff)
            stats = get_stats(sasview_value, gpu_double_value, index)
            columns.extend(stats)
            labels.append('GPU double')
            good = good and (stats[0] < 1e-14)
        if 1:
            cpu_double_value = trymodel(eval_ctypes, dtype='double', cutoff=cutoff)
            stats = get_stats(sasview_value, cpu_double_value, index)
            columns.extend(stats)
            labels.append('CPU double')
            good = good and (stats[0] < 1e-14)
        if 0:
            stats = get_stats(cpu_double_value, gpu_single_value, index)
            columns.extend(stats)
            labels.append('single/double')
            good = good and (stats[0] < 1e-14)

        columns += [v for _,v in sorted(pars.items())]
        if first:
            print_column_headers(pars, labels)
            first = False
        if good:
            num_good += 1
        else:
            print(("%d,"%seed)+','.join("%g"%v for v in columns))
    print '"%d/%d good"'%(num_good, N)


def print_usage():
    print "usage: compare_many.py MODEL COUNT (1dNQ|2dNQ) (CUTOFF|mono)"


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

Available models:
""")
    print_models()

def main():
    if len(sys.argv) == 1:
        print_help()
        sys.exit(1)

    model = sys.argv[1]
    if not (model in MODELS) and (model != "all"):
        print 'Bad model %s.  Use "all" or one of:'
        print_models()
        sys.exit(1)
    try:
        count = int(sys.argv[2])
        is2D = sys.argv[3].startswith('2d')
        assert sys.argv[3][1] == 'd'
        Nq = int(sys.argv[3][2:])
        mono = sys.argv[4] == 'mono'
        cutoff = float(sys.argv[4]) if not mono else 0
    except:
        print_usage()
        sys.exit(1)

    data, index = make_data(qmax=1.0, is2D=is2D, Nq=Nq)
    model_list = [model] if model != "all" else MODELS
    for model in model_list:
        compare_instance(model, data, index, N=count, mono=mono, cutoff=cutoff)

if __name__ == "__main__":
    main()
