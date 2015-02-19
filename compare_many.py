#!/usr/bin/env python

import sys

import numpy as np

from sasmodels.kernelcl import environment
from compare import (MODELS, randomize_model, suppress_pd, eval_sasview,
                     eval_opencl, eval_ctypes, make_data, get_demo_pars)

def get_stats(target, value, index):
    resid = abs(value-target)[index]
    relerr = resid/target[index]
    srel = np.argsort(relerr)
    p90 = int(len(relerr)*0.90)
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
    columns = ['Seed'] + stats*len(parts) +  list(sorted(pars.keys()))
    print(','.join('"%s"'%c for c in groups))
    print(','.join('"%s"'%c for c in columns))

def compare_instance(name, data, index, N=1, mono=True, cutoff=1e-5):
    pars = get_demo_pars(name)
    header = '\n"Model","%s","Count","%d"'%(name, N)
    if not mono: header += ',"Cutoff",%g'%(cutoff,)
    print(header)
    first = True
    for _ in range(N):
        pars, seed = randomize_model(name, pars)
        if mono: suppress_pd(pars)

        target, _ = eval_sasview(name, pars, data)

        env = environment()
        gpu_single_value,_ = eval_opencl(name, pars, data, dtype='single', cutoff=cutoff)
        gpu_single = get_stats(target, gpu_single_value, index)
        if env.has_double:
            gpu_double_value,_ = eval_opencl(name, pars, data, dtype='double', cutoff=cutoff)
            gpu_double = get_stats(target, gpu_double_value, index)
        else:
            gpu_double = [0]*len(gpu_single)
        cpu_double_value,_ =  eval_ctypes(name, pars, data, dtype='double', cutoff=cutoff)
        cpu_double = get_stats(target, cpu_double_value, index)
        single_double = get_stats(cpu_double_value, gpu_single_value, index)

        values = (list(gpu_single) + list(gpu_double) + list(cpu_double)
                  + list(single_double) + [v for _,v in sorted(pars.items())])
        if gpu_single[0] > 5e-5:
            if first:
                print_column_headers(pars,'GPU single|GPU double|CPU double|single/double'.split('|'))
                first = False
            print(("%d,"%seed)+','.join("%g"%v for v in values))

def main():
    try:
        model = sys.argv[1]
        assert (model in MODELS) or (model == "all")
        count = int(sys.argv[2])
        is2D = sys.argv[3].startswith('2d')
        assert sys.argv[3][1] == 'd'
        Nq = int(sys.argv[3][2:])
        mono = sys.argv[4] == 'mono'
        cutoff = float(sys.argv[4]) if not mono else 0
    except:
        import traceback; traceback.print_exc()
        models = "\n    ".join("%-7s: %s"%(k,v.__name__.replace('_',' '))
                               for k,v in sorted(MODELS.items()))
        print("""\
usage: compare_many.py MODEL COUNT (1dNQ|2dNQ) (CUTOFF|mono)

MODEL is the model name of the model, which is one of:
    %s
or "all" for all the models in alphabetical order.

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
"""%(models,))
        sys.exit(1)

    data, index = make_data(qmax=1.0, is2D=is2D, Nq=Nq)
    model_list = [model] if model != "all" else MODELS
    for model in model_list:
        compare_instance(model, data, index, N=count, mono=mono, cutoff=cutoff)

if __name__ == "__main__":
    main()
