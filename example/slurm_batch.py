#!/usr/bin/env python
"""
Submit a batch fit job to the slurm cluster.

Given a model.py file defining a Bumps problem defined on a single data
file, with the data file specified as a command line argument, run the
bumps fit as a batch over a set of different datafiles independently.
An example model is given in model_ellipsoid_hayter_msa.py, which fits 
the data in 09319*.dat.

To run the fit, use::

    slurm_batch.py [--slurm_opts] model.py *.dat --store=T1 [--bumps_opt ...]

For example::

    slurm_batch.py model_ellipsoid_hayter_msa.py 09319*.dat --store=T1

You may need to run in a particular python environment on the 
compute nodes::

    slurm_batch.py --python=path/to/env/bin/python ...

This creates the T1 subdirectory to hold the fit results and 
prints the real command that is submitted, as well as the job id.

The store directory T1 contains a copy of the model file and 
all the data files.  The fit results for each file will be 
in T1/##/*.  The file T1/files.dat contains the list 
of "subdirectory filename" pairs indicating which ## directory 
contains the resuls for which file.  Check for errors using::

    cat T1/slurm*_1.out

The following slurm options are used::

    --array=1-#files     batch size comes from the file list
    --gres=gpu:1         request a gpu for each fit
    --job-name=model.py  use model file name for job name
    --output=...         log into T1/slurm-job_##.out
    --chdir=...          run fit from store directory
    --time=2             time as number of hours (can override)

To receive an email on job completion or failure, add the following
slurm options before the model file::

    --mail-type=END,FAIL --mail-user=user@mail.domain

Bumps options are described at bumps.readthedocs.org, with the
following set automatically::

    --batch              run in batch mode, without output to .mon
    --view=log           SAS fits want log plots
    --time=2-0.1         slurm time minus 6 minutes for cleanup

The --store and --resume options indicate the parent directory for
the output.  These are modified to store the results in a separate
subdirectory for each file.  Keep in mind that the fit is run from
the store directory, so any files or modules referenced from the
model file will need to use a full path to the original location.

After submitting the job a job id will be printed to the console.
You can check the status of the job using the usual slurm commands
such as::

    squeue

or cancel the job using::

    scancel jobid

The slurm_batch program runs directly from the source tree for sasmodels, 
and requires sasview, bumps and periodictable as sister directories
accessible on the worker nodes.  You can link it into your bin directory
using::

     mkdir ~/bin
     ln -s path/to/slurm_batch.py ~/bin

or if you are a cluster administrator, into /usr/local/bin.
"""

# If called from command line, this submits a job to the slurm queue, with _this_ file
# as the batch script. Before calling it on the worker node, slurm sets the
# SLURM_ARRAY_TASK_ID to the current task so we can tell that we are running
# as a worker and which file we should be working on.

## SBATCH options as comments do not seem to work.  Maybe they neeed to be before 
## the doc string?  For now they are hardcoded in the sbatch call in submit_job.

import sys
import os
import tempfile
import shutil

DEFAULT_TIME_LIMIT = 2

def split_args():
    slurm_opts = []
    bumps_opts = []
    model_file = None
    store = None
    resume = None
    data_files = []
    time_limit = DEFAULT_TIME_LIMIT
    interpreter = sys.executable

    # start with '-' arguments as slurm opts, then after
    # the model file any '-' arguments are bumps opts.
    opts = slurm_opts
    for v in sys.argv[1:]:
        if v.startswith('--store='):
            store = os.path.realpath(os.path.abspath(v[8:]))
        elif v.startswith('--resume='):
            resume = os.path.realpath(os.path.abspath(v[9:]))
        elif v.startswith('--time='):
            time_limit = float(v[7:])
        elif v.startswith('--python='):
            interpreter = v[9:]
        elif v[0] == '-':
            opts.append(v)
        elif model_file is None:
            model_file = v
            opts = bumps_opts
        else:
            data_files.append(v)


    s = time_limit*3600
    slurm_opts.append("--time=%d:%02d:%02d"%(s//3600, (s%3600)//60, s%60))
    bumps_opts.append('--time=%f'%(time_limit - 0.1))  # 6 min to stop cleanly

    return {
        'python': interpreter,
        'slurm': slurm_opts, 
        'model_file': model_file, 
        'data_files': data_files, 
        'store': store, 
        'resume': resume,
        'bumps': bumps_opts,
    }

def dirn(path, n):
    path = os.path.realpath(os.path.abspath(path))
    for _ in range(n):
        path = os.path.dirname(path)
    return path

def submit_job():
    # sbatch --array=1-5 ./slurm_batch.py model_ellipsoid_hayter_msa.py 09*.dat --store=T1 --fit=dream
    opts = split_args()
    store = opts['store']
    model_file = opts['model_file']
    data_files = opts['data_files']
    bumps_opts = opts['bumps']
    slurm_opts = opts['slurm']
    interpreter = opts['python']

    # make sure the store directory exists and save the order of the files, as well
    # as the model and the data files
    if store is not None:
        if not os.path.exists(store):
            os.makedirs(store)

        # save file order
        with open(os.path.join(store, 'files.dat'), 'w') as fid:
            for k, f in enumerate(data_files):
                fid.write("%02d %s\n"%(k+1, f))

        # Copy the model and data files to the root store directory
        # Since bumps changes into the model directory prior to loading
        # the datafiles, strip all leading paths from data and model and
        # set the working directory for the job to the store directory.
        model_copy = os.path.basename(model_file)
        shutil.copy(model_file, os.path.join(store, model_copy))
        data_copy = []
        for f in data_files:
            f_copy = os.path.basename(f)
            shutil.copy(f, os.path.join(store, f_copy))
            data_copy.append(f_copy) 

        model_file = model_copy
        data_files = data_copy


    # build and run the command
    SRC = dirn(__file__, 3) # __file__ is $SRC/sasmodels/example/slurm_batch.py
    parts = [
        "sbatch",
        "--array=1-%d"%len(data_files),
        "--gres=gpu:1",
        "--job-name="+model_file,
        ## since we are setting the current working directory, we don't need
        ## to fiddle the slurm output files
        "--output=%s/slurm-%%A_%%a.out"%store,
        "--chdir=%s"%store,
        ]
    parts.extend(slurm_opts)
    parts.append(__file__)
    # Remember the source root so we can reconstruct the correct python path
    # This is done after the model file so that it doesn't get interpreted
    # as a slurm option.
    parts.append("--python=%s"%interpreter)
    parts.append("--source_root=%s"%SRC)
    parts.append(model_file)
    parts.extend(data_files)
    parts.extend(bumps_opts)
    #if store is not None:
    #    parts.append("--store=" + store)
    command = " ".join(parts)

    print(command)
    os.system(command)

def run_task(task_id):
    opts = split_args()

    # Set environment put compiled sasmodels in user-specific temporary cache
    # We need this because users don't have a home directory on the individual
    # cluster nodes.
    assert opts['slurm'][0].startswith('--source_root=')
    SRC = opts['slurm'][0][14:]
    PACKAGES = ("periodictable", "sasview/src", "bumps", "sasmodels")
    os.environ['PYTHONPATH'] = ":".join(SRC+"/"+v for v in PACKAGES)
    TMP = tempfile.gettempdir()
    cache_path = os.path.join(TMP, os.environ['USER'], '.cache')
    os.environ['SAS_DLL_PATH'] = cache_path
    os.environ['XDG_CACHE_HOME'] = cache_path

    #task_store = "%s/%02d"%(opts['store'], task_id)
    task_store = "%02d"%task_id
    parts = [
       opts['python'], os.path.join(SRC, "bumps", "run.py"), "--batch",
       "--view=log",
       opts['model_file'],
       opts['data_files'][task_id-1],
       ]
    parts.extend(opts['bumps'])
    parts.append('--store='+task_store)
    if opts['resume'] is not None:
        parts.append('--resume='+os.path.join(opts['resume'], task_store))
    command = " ".join(parts)
    print(os.getcwd() + "$ " + command)
    os.system(command)


task_id = int(os.environ.get('SLURM_ARRAY_TASK_ID', -1))
if task_id == -1:
    submit_job()
else:
    run_task(task_id)

