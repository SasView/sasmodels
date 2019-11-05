from __future__ import print_function

import sys
import os
from os.path import basename, dirname, realpath, join as joinpath, exists
import math
import re
import shutil

# CRUFT: python 2.7 backport of makedirs(path, exist_ok=False)
if sys.version_info[0] >= 3:
    from os import makedirs
else:
    def makedirs(path, exist_ok=False):
        try:
            os.makedirs(path)
        except Exception:
            if not exist_ok or not exists(path):
                raise

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, realpath(joinpath(dirname(__file__), '..')))
import sasmodels
from sasmodels import generate, core
from sasmodels.direct_model import DirectModel, call_profile
from sasmodels.data import empty_data1D, empty_data2D

try:
    from typing import Dict, Any
    from matplotlib.axes import Axes
    from sasmodels.kernel import KernelModel
    from sasmodels.modelinfo import ModelInfo
except ImportError:
    pass

def plot_1d(model, opts, ax):
    # type: (KernelModel, Dict[str, Any], Axes) -> None
    """
    Create a 1-D image.
    """
    q_min, q_max, nq = opts['q_min'], opts['q_max'], opts['nq']
    q_min = math.log10(q_min)
    q_max = math.log10(q_max)
    q = np.logspace(q_min, q_max, nq)
    data = empty_data1D(q)
    calculator = DirectModel(data, model)
    Iq1D = calculator()

    ax.plot(q, Iq1D, color='blue', lw=2, label=model.info.name)
    ax.set_xlabel(r'$Q \/(\AA^{-1})$')
    ax.set_ylabel(r'$I(Q) \/(\mathrm{cm}^{-1})$')
    ax.set_xscale(opts['xscale'])
    ax.set_yscale(opts['yscale'])
    #ax.legend(loc='best')

def plot_2d(model, opts, ax):
    # type: (KernelModel, Dict[str, Any], Axes) -> None
    """
    Create a 2-D image.
    """
    qx_max, nq2d = opts['qx_max'], opts['nq2d']
    q = np.linspace(-qx_max, qx_max, nq2d) # type: np.ndarray
    data2d = empty_data2D(q, resolution=0.0)
    calculator = DirectModel(data2d, model)
    Iq2D = calculator() #background=0)
    Iq2D = Iq2D.reshape(nq2d, nq2d)
    if opts['zscale'] == 'log':
        Iq2D = np.log(np.clip(Iq2D, opts['vmin'], np.inf))
    ax.imshow(Iq2D, interpolation='nearest', aspect=1, origin='lower',
              extent=[-qx_max, qx_max, -qx_max, qx_max], cmap=opts['colormap'])
    ax.set_xlabel(r'$Q_x \/(\AA^{-1})$')
    ax.set_ylabel(r'$Q_y \/(\AA^{-1})$')

def plot_profile_inset(model_info, ax):
    # type: (ModelInfo, Axes) -> None
    """
    Plot 1D radial profile as inset plot.
    """
    p = ax.get_position()
    width, height = 0.4*(p.x1-p.x0), 0.4*(p.y1-p.y0)
    left, bottom = p.x1-width, p.y1-height
    inset = plt.gcf().add_axes([left, bottom, width, height])
    x, y, labels = call_profile(model_info)
    inset.plot(x, y, '-')
    inset.locator_params(nbins=4)
    #inset.set_xlabel(labels[0])
    #inset.set_ylabel(labels[1])
    inset.text(0.99, 0.99, "profile",
               horizontalalignment="right",
               verticalalignment="top",
               transform=inset.transAxes)

def figfile(model_info):
    # type: (ModelInfo) -> str
    return model_info.id + '_autogenfig.png'

def make_figure(model_info, opts):
    # type: (ModelInfo, Dict[str, Any]) -> None
    """
    Generate the figure file to include in the docs.
    """
    model = core.build_model(model_info)

    fig_height = 3.0 # in
    fig_left = 0.6 # in
    fig_right = 0.5 # in
    fig_top = 0.6*0.25 # in
    fig_bottom = 0.6*0.75
    if model_info.parameters.has_2d:
        plot_height = fig_height - (fig_top+fig_bottom)
        plot_width = plot_height
        fig_width = 2*(plot_width + fig_left + fig_right)
        aspect = (fig_width, fig_height)
        ratio = aspect[0]/aspect[1]
        ax_left = fig_left/fig_width
        ax_bottom = fig_bottom/fig_height
        ax_height = plot_height/fig_height
        ax_width = ax_height/ratio # square axes
        fig = plt.figure(figsize=aspect)
        ax2d = fig.add_axes([0.5+ax_left, ax_bottom, ax_width, ax_height])
        plot_2d(model, opts, ax2d)
        ax1d = fig.add_axes([ax_left, ax_bottom, ax_width, ax_height])
        plot_1d(model, opts, ax1d)
        #ax.set_aspect('square')
    else:
        plot_height = fig_height - (fig_top+fig_bottom)
        plot_width = (1+np.sqrt(5))/2*fig_height
        fig_width = plot_width + fig_left + fig_right
        ax_left = fig_left/fig_width
        ax_bottom = fig_bottom/fig_height
        ax_width = plot_width/fig_width
        ax_height = plot_height/fig_height
        aspect = (fig_width, fig_height)
        fig = plt.figure(figsize=aspect)
        ax1d = fig.add_axes([ax_left, ax_bottom, ax_width, ax_height])
        plot_1d(model, opts, ax1d)

    if model_info.profile:
        plot_profile_inset(model_info, ax1d)

    # Save image in model/img
    path = joinpath('model', 'img')
    makedirs(path, exist_ok=True)
    path = joinpath('model', 'img', figfile(model_info))
    plt.savefig(path, bbox_inches='tight')
    #print("figure saved in",path)

def copy_if_newer(src, dst):
    # type: (str) -> str
    """
    Copy from *src* to *dst* if *src* is newer or *dst* doesn't exist.
    """
    if not exists(dst):
        path = dirname(dst)
        makedirs(path, exist_ok=True)
        shutil.copy2(src, dst)
    elif os.path.getmtime(src) > os.path.getmtime(dst):
        shutil.copy2(src, dst)

def link_sources(model_info):
    # type: (ModelInfo) -> str
    """
    Add link to model sources from the doc tree.
    """

    # List source files in reverse order of dependency.
    model_file = basename(model_info.filename)
    sources = list(reversed(model_info.source + [model_file]))

    # Copy files to src dir under models directory.  Need to do this
    # because sphinx can't link to an absolute path.
    root = dirname(dirname(realpath(__file__)))
    src = joinpath(root, "sasmodels", "models")
    dst = joinpath(root, "doc", "model", "src")
    for path in sources:
        copy_if_newer(joinpath(src, path), joinpath(dst, path))

    # Link to local copy of the files
    downloads = [":download:`%s <src/%s>`"%(path, path) for path in sources]

    # Could do syntax highlighting on the model files by creating a rst file
    # beside each source file named containing source file with
    #
    #    src/path.rst:
    #
    #    .. {{ path.replace('/','_') }}:
    #
    #    .. literalinclude:: {{ src/path }}
    #        :language: {{ "python" if path.endswith('.py') else "c" }}
    #        :linenos:
    #
    # and link to it using
    #
    #     colors = [":ref:`%s`"%(path.replace('/','_')) for path in sources]
    #
    # Probably need to dump all the rst files into an index.rst to build them.

    # Link to github repo (either the tagged sasmodels version or master)
    #url = "https://github.com/SasView/sasmodels/blob/v%s"%sasmodels.__version__
    #url = "https://github.com/SasView/sasmodels/blob/master"%sasmodels.__version__
    #links = ["`%s <%s/sasmodels/models/%s>`_"%(path, url, path) for path in sources]

    sep = "\n$\\ \\star\\ $ "  # bullet
    body = "\n**Source**\n"
    #body += "\n" + sep.join(links) + "\n\n"
    body += "\n" + sep.join(downloads) + "\n\n"
    return body

def gen_docs(model_info, outfile):
    # type: (ModelInfo, str) -> None
    """
    Generate the doc string with the figure inserted before the references.
    """

    # Load the doc string from the module definition file and store it in rst
    docstr = generate.make_doc(model_info)

    # Auto caption for figure
    captionstr = '\n'
    captionstr += '.. figure:: img/' + figfile(model_info) + '\n'
    captionstr += '\n'
    if model_info.parameters.has_2d:
        captionstr += '    1D and 2D plots corresponding to the default parameters of the model.\n'
    else:
        captionstr += '    1D plot corresponding to the default parameters of the model.\n'
    captionstr += '\n'

    # Add figure reference and caption to documentation (at end, before References)
    pattern = r'\*\*REFERENCE'
    match = re.search(pattern, docstr.upper())

    sources = link_sources(model_info)

    insertion = captionstr + sources

    if match:
        docstr1 = docstr[:match.start()]
        docstr2 = docstr[match.start():]
        docstr = docstr1 + insertion + docstr2
    else:
        print('------------------------------------------------------------------')
        print('References NOT FOUND for model: ', model_info.id)
        print('------------------------------------------------------------------')
        docstr += insertion

    with open(outfile, 'w') as fid:
        fid.write(docstr)

def process_model(infile, outfile):
    # type: (str, str) -> None
    """
    Generate doc file and image file for the given model definition file.
    """
    # Load the model file
    model_name = basename(infile)[:-3]
    model_info = core.load_model_info(model_name)

    # Plotting ranges and options
    opts = {
        'xscale'    : 'log',
        'yscale'    : 'log' if not model_info.structure_factor else 'linear',
        'zscale'    : 'log' if not model_info.structure_factor else 'linear',
        'q_min'     : 0.001,
        'q_max'     : 1.0,
        'nq'        : 1000,
        'nq2d'      : 1000,
        'vmin'      : 1e-3,  # floor for the 2D data results
        'qx_max'    : 0.5,
        #'colormap'  : 'gist_ncar',
        'colormap'  : 'nipy_spectral',
        #'colormap'  : 'jet',
    }

    # Generate the RST file and the figure.  Order doesn't matter.
    gen_docs(model_info, outfile)
    make_figure(model_info, opts)

def main():
    infile, outfile = sys.argv[1], sys.argv[2]
    process_model(infile, outfile)

if __name__ == "__main__":
    main()
