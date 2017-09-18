#!/usr/bin/env python
"""
Script to run a batch fit in a series of files and plot the fitted parameters.

Usage syntax::

    python batch_fit.py model.py sample1.dat sample2.dat ... other_sample.dat
    (files named sample1.dat, sample2.dat, ..., other_sample.dat)

or if the file names are numbers (and the extension is .dat)::

    python batch_fit.py model.py 93190 93210
    (files named 093190.dat, 093191.dat, ..., 093210.dat)

or for Grasp-like naming::

    python batch_fit.py model.py 93190 93210 200
    (files named 093190_200.dat, 093191_201.dat, ..., 093210_220.dat)

The script reads a series of files and fits the model defined by model.py.
For example model_ellipsoid_hayter_msa.py fits a model consisting in an
ellipsoid form factor multiplied by a Hayter MSA structure factor.

The file *model.py* must load the data using::

    data = load_data(sys.argv[1])

Include options to bumps (minimizer, steps, etc.) as desired.  For example::

    python batch_fit.py model.py 93190 93210 200 --fit=lm --steps=200 --ftol=1.5e-8 --xtol=1.5e-8

Fit options can come before or after the model and files.

For each file a directory named Fit_filename is created. There the file
model.par contains the fitted parameters.  These are gathered together into
batch_fit.csv in the current directory.

Finally the fitted parameters are plotted for the full series.

Example::

    python batch_fit.py model_ellipsoid_hayter_msa.py 93191 93195 201

Note:

    If sasmodels, sasview or bumps are not in the path, use the PYTHONPATH
    environment variable to set them.
"""
from __future__ import print_function

import sys
import os

import numpy as np
import matplotlib.pyplot as plt
from bumps.dream.views import tile_axes  # make a grid of plots

# GET INPUT AND ENSURE MODEL AND DATA FILES ARE DEFINED

fit_opts = [v for v in sys.argv[1:] if v.startswith('--')]
args = [v for v in sys.argv[1:] if not v.startswith('--')]

nargs = len(args)
if nargs < 2:
    print ("Error in the list of arguments! \n")
    sys.exit()

model_file = args[0]
if not model_file.endswith('.py'):
    print("Expected model.py as the first argument")
    sys.exit(1)

if '.' in args[1]:
    data_files = args[1:]
else:
    first = int(args[1])
    last = int(args[2])
    count = last-first+1
    data_files = []
    if nargs == 3:
        data_files = ['%06d.dat'%(first+i) for i in range(count)]
    elif nargs == 4:
        ext = int(args[3])
        data_files = ['%06d_%d.dat'%(first+i, ext+i) for i in range(count)]
    else:
        print("Unexpected arguments: " + " ".join(args[4:]))
        sys.exit(1)

# CHECK THAT DATA FILES EXIST
missing = [filename for filename in data_files if not os.path.isfile(filename)]
if missing:
    print("Missing data files: %s" % ", ".join(missing))
    sys.exit(1)

# STORE DIRECTORY FOR BUMPS FITS
def fit_dir(filename):
    "Return the store directory name for the given file"
    return "Fit_" + os.path.splitext(filename)[0]

# LOOP OVER FILES AND CALL TO BUMPS FOR EACH OF THEM
bumps_cmd = "python -m bumps.cli --batch"
fit_opts = " ".join(fit_opts)
for data_file in data_files:
    store_opts = "--store=" + fit_dir(data_file)
    cmd = " ".join((bumps_cmd, fit_opts, store_opts, model_file, data_file))
    os.system(cmd)

# GATHER RESULTS
results = {}
par_file = os.path.splitext(model_file)[0] + '.par'
for data_file in data_files:
    with open(os.path.join(fit_dir(data_file), par_file), 'r') as fid:
        for line in fid:
            parameter, value = line.split()
            results.setdefault(parameter, []).append(float(value))

# SAVE RESULTS INTO FILE
with open('batch_fit.csv', 'w') as fid:
    parameters = list(sorted(results.keys()))
    values_by_file = zip(*(v for k, v in sorted(results.items())))
    fid.write(','.join(['filename'] + parameters) + '\n')
    for filename, values in zip(data_files, values_by_file):
        fid.write(','.join([filename] + [str(v) for v in values]) + '\n')

# SHOW FITTED PARAMETERS
nh, nw = tile_axes(len(results))
ticks = np.arange(1, len(data_files)+1)
labels = [os.path.splitext(filename)[0] for filename in data_files]
for k, (parameter, values) in enumerate(sorted(results.items())):
    plt.subplot(nh, nw, k+1)
    plt.plot(ticks, values)
    plt.xlim(ticks[0]-0.5, ticks[-1]+0.5)
    if k%nh == nh-1:
        #plt.xlabel('Dataset #')
        plt.xticks(ticks, labels, rotation=30)
    else:
        plt.xticks(ticks, [' ']*len(labels))
    plt.ylabel(parameter)
plt.suptitle("Fit " + args[0])
plt.show()
