#!/usr/bin/env python
"""
Program to compare results from many random parameter sets for a given model.

The result is a comma separated value (CSV) table that can be redirected
from standard output into a file and loaded into a spreadsheet.

The models are compared for each parameter set and if the difference is
greater than expected for that precision, the parameter set is labeled
as bad and written to the output, along with the random seed used to
generate that parameter value.  This seed can be used with :mod:`compare`
to reload and display the details of the model.
"""
from __future__ import print_function

import sys
import traceback

import numpy as np  # type: ignore

from . import core
from .compare import (randomize_pars, suppress_pd, make_data,
                      make_engine, get_pars, columnize,
                      constrain_pars)

MODELS = core.list_models()

def calc_stats(target, value, index):
    """
    Calculate statistics between the target value and the computed value.

    *target* and *value* are the vectors being compared, with the
    difference normalized by *target* to get relative error.  Only
    the elements listed in *index* are used, though index may be
    and empty slice defined by *slice(None, None)*.

    Returns:

        *maxrel* the maximum relative difference

        *rel95* the relative difference with the 5% biggest differences ignored

        *maxabs* the maximum absolute difference for the 5% biggest differences

        *maxval* the maximum value for the 5% biggest differences
    """
    resid = abs(value-target)[index]
    relerr = resid/target[index]
    sorted_rel_index = np.argsort(relerr)
    #p90 = int(len(relerr)*0.90)
    p95 = int(len(relerr)*0.95)
    maxrel = np.max(relerr)
    rel95 = relerr[sorted_rel_index[p95]]
    maxabs = np.max(resid[sorted_rel_index[p95:]])
    maxval = np.max(value[sorted_rel_index[p95:]])
    return maxrel, rel95, maxabs, maxval

def print_column_headers(pars, parts):
    """
    Generate column headers for the differences and for the parameters,
    and print them to standard output.
    """
    stats = list('Max rel err|95% rel err|Max abs err above 90% rel|Max value above 90% rel'.split('|'))
    groups = ['']
    for p in parts:
        groups.append(p)
        groups.extend(['']*(len(stats)-1))
    groups.append("Parameters")
    columns = ['Seed'] + stats*len(parts) +  list(sorted(pars.keys()))
    print(','.join('"%s"'%c for c in groups))
    print(','.join('"%s"'%c for c in columns))

# Target 'good' value for various precision levels.
PRECISION = {
    'fast': 1e-3,
    'half': 1e-3,
    'single': 5e-5,
    'double': 5e-14,
    'single!': 5e-5,
    'double!': 5e-14,
    'quad!': 5e-18,
}
def compare_instance(name, data, index, N=1, mono=True, cutoff=1e-5,
                     base='single', comp='double'):
    r"""
    Compare the model under different calculation engines.

    *name* is the name of the model.

    *data* is the data object giving $q, \Delta q$ calculation points.

    *index* is the active set of points.

    *N* is the number of comparisons to make.

    *cutoff* is the polydispersity weight cutoff to make the calculation
    a little bit faster.

    *base* and *comp* are the names of the calculation engines to compare.
    """

    is_2d = hasattr(data, 'qx_data')
    model_info = core.load_model_info(name)
    pars = get_pars(model_info, use_demo=True)
    header = ('\n"Model","%s","Count","%d","Dimension","%s"'
              % (name, N, "2D" if is_2d else "1D"))
    if not mono:
        header += ',"Cutoff",%g'%(cutoff,)
    print(header)

    if is_2d:
        if not model_info.parameters.has_2d:
            print(',"1-D only"')
            return

    # Some not very clean macros for evaluating the models and checking the
    # results.  They freely use variables from the current scope, even some
    # which have not been defined yet, complete with abuse of mutable lists
    # to allow them to update values in the current scope since nonlocal
    # declarations are not available in python 2.7.
    def try_model(fn, pars):
        """
        Return the model evaluated at *pars*.  If there is an exception,
        print it and return NaN of the right shape.
        """
        try:
            result = fn(**pars)
        except Exception:
            traceback.print_exc()
            print("when comparing %s for %d"%(name, seed))
            if hasattr(data, 'qx_data'):
                result = np.NaN*data.data
            else:
                result = np.NaN*data.x
        return result
    def check_model(pars):
        """
        Run the two calculators against *pars*, returning statistics
        on the differences.  See :func:`calc_stats` for the list of stats.
        """
        base_value = try_model(calc_base, pars)
        comp_value = try_model(calc_comp, pars)
        stats = calc_stats(base_value, comp_value, index)
        max_diff[0] = max(max_diff[0], stats[0])
        good[0] = good[0] and (stats[0] < expected)
        return list(stats)


    try:
        calc_base = make_engine(model_info, data, base, cutoff)
        calc_comp = make_engine(model_info, data, comp, cutoff)
    except Exception as exc:
        #raise
        print('"Error: %s"'%str(exc).replace('"', "'"))
        print('"good","%d of %d","max diff",%g' % (0, N, np.NaN))
        return
    expected = max(PRECISION[base], PRECISION[comp])

    num_good = 0
    first = True
    max_diff = [0]
    for k in range(N):
        print("Model %s %d"%(name, k+1), file=sys.stderr)
        seed = np.random.randint(1e6)
        np.random.seed(seed)
        pars_i = randomize_pars(model_info, pars)
        constrain_pars(model_info, pars_i)
        if mono:
            pars_i = suppress_pd(pars_i)

        good = [True]
        columns = check_model(pars_i)
        columns += [v for _, v in sorted(pars_i.items())]
        if first:
            labels = [" vs. ".join((calc_base.engine, calc_comp.engine))]
            print_column_headers(pars_i, labels)
            first = False
        if good[0]:
            num_good += 1
        else:
            print(("%d,"%seed)+','.join("%s"%v for v in columns))
    print('"good","%d of %d","max diff",%g'%(num_good, N, max_diff[0]))


def print_usage():
    """
    Print the command usage string.
    """
    print("usage: compare_many.py MODEL COUNT (1dNQ|2dNQ) (CUTOFF|mono) (single|double|quad)",
          file=sys.stderr)


def print_models():
    """
    Print the list of available models in columns.
    """
    print(columnize(MODELS, indent="  "))


def print_help():
    """
    Print usage string, the option description and the list of available models.
    """
    print_usage()
    print("""\

MODEL is the model name of the model or one of the model types listed in
sasmodels.core.list_models (all, py, c, double, single, opencl, 1d, 2d,
nonmagnetic, magnetic).  Model types can be combined, such as 2d+single.

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
is set in compare.py defaults for each model.  Polydispersity is given in the
"demo" attribute of each model.

PRECISION is the floating point precision to use for comparisons.  If two
precisions are given, then compare one to the other.  Precision is one of
fast, single, double for GPU or single!, double!, quad! for DLL.  If no
precision is given, then use single and double! respectively.

Available models:
""")
    print_models()

def main(argv):
    """
    Main program.
    """
    if len(argv) not in (3, 4, 5, 6):
        print_help()
        return

    target = argv[0]
    try:
        model_list = [target] if target in MODELS else core.list_models(target)
    except ValueError:
        print('Bad model %s.  Use model type or one of:' % target, file=sys.stderr)
        print_models()
        print('model types: all, py, c, double, single, opencl, 1d, 2d, nonmagnetic, magnetic')
        return
    try:
        count = int(argv[1])
        is2D = argv[2].startswith('2d')
        assert argv[2][1] == 'd'
        Nq = int(argv[2][2:])
        mono = len(argv) <= 3 or argv[3] == 'mono'
        cutoff = float(argv[3]) if not mono else 0
        base = argv[4] if len(argv) > 4 else "single"
        comp = argv[5] if len(argv) > 5 else "double!"
    except Exception:
        traceback.print_exc()
        print_usage()
        return

    data, index = make_data({
        'qmin': 0.001, 'qmax': 1.0, 'is2d': is2D, 'nq': Nq, 'res': 0.,
        'accuracy': 'Low', 'view':'log', 'zero': False
        })
    for model in model_list:
        compare_instance(model, data, index, N=count, mono=mono,
                         cutoff=cutoff, base=base, comp=comp)

if __name__ == "__main__":
    #from .compare import push_seed
    #with push_seed(1): main(sys.argv[1:])
    main(sys.argv[1:])
