import sys, os, math, re
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath('..'))
from sasmodels import generate, core
from sasmodels.direct_model import DirectModel
from sasmodels.data import empty_data1D, empty_data2D


# Convert ../sasmodels/models/name.py to name
model_name = os.path.basename(sys.argv[1])[:-3]
model_info = core.load_model_info(model_name)
model = core.build_model(model_info)

# Load the doc string from the module definition file and store it in rst
docstr = generate.make_doc(model_info)


# Calculate 1D curve for default parameters
pars = dict((p.name, p.default) for p in model_info['parameters'])

# Plotting ranges and options
opts = {
    'xscale'    : 'log',
    'yscale'    : 'log' if not model_info['structure_factor'] else 'linear',
    'zscale'    : 'log' if not model_info['structure_factor'] else 'linear',
    'q_min'     : 0.001,
    'q_max'     : 1.0,
    'nq'        : 1000,
    'nq2d'      : 400,
    'vmin'      : 1e-3,  # floor for the 2D data results
    'qx_max'    : 0.5,
}


def plot_1d(model, opts, ax):
    q_min, q_max, nq = opts['q_min'], opts['q_max'], opts['nq']
    q_min = math.log10(q_min)
    q_max = math.log10(q_max)
    q = np.logspace(q_min, q_max, nq)
    data = empty_data1D(q)
    calculator = DirectModel(data, model)
    Iq1D = calculator()

    ax.plot(q, Iq1D, color='blue', lw=2, label=model_info['name'])
    ax.set_xlabel(r'$Q \/(\AA^{-1})$')
    ax.set_ylabel(r'$I(Q) \/(\mathrm{cm}^{-1})$')
    ax.set_xscale(opts['xscale'])
    ax.set_yscale(opts['yscale'])
    #ax.legend(loc='best')

def plot_2d(model, opts, ax):
    qx_max, nq2d = opts['qx_max'], opts['nq2d']
    q = np.linspace(-qx_max, qx_max, nq2d)
    data2d = empty_data2D(q, resolution=0.0)
    calculator = DirectModel(data2d, model)
    Iq2D = calculator() #background=0)
    Iq2D = Iq2D.reshape(nq2d, nq2d)
    if opts['zscale'] == 'log':
        Iq2D = np.log(np.clip(Iq2D, opts['vmin'], np.inf))
    h = ax.imshow(Iq2D, interpolation='nearest', aspect=1, origin='upper',
           extent=[-qx_max, qx_max, -qx_max, qx_max], cmap=ice_cm())
    # , vmin=vmin, vmax=vmax)
    ax.set_xlabel(r'$Q_x \/(\AA^{-1})$')
    ax.set_ylabel(r'$Q_y \/(\AA^{-1})$')

def ice_cm():
    from matplotlib._cm import _Blues_data
    from matplotlib import colors
    from matplotlib import rcParams
    def from_white(segments):
        scale = 1.0/segments[0][1]
        return [(k, v*scale, w*scale) for k, v, w in segments]
    ice_data = dict((k,from_white(v)) for k,v in _Blues_data.items())
    ice = colors.LinearSegmentedColormap("ice", ice_data, rcParams['image.lut'])
    return ice


# Generate image (comment IF for 1D/2D for the moment) and generate only 1D
fig_height = 3.0 # in
fig_left = 0.6 # in
fig_right = 0.5 # in
fig_top = 0.6*0.25 # in
fig_bottom = 0.6*0.75
if model_info['has_2d']:
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

# Save image in model/img
figname = model_name + '_autogenfig.png'
filename = os.path.join('model', 'img', figname)
plt.savefig(filename, bbox_inches='tight')
#print "figure saved in",filename

# Auto caption for figure
captionstr = '\n'
captionstr += '.. figure:: img/' + model_info['id'] + '_autogenfig.png\n'
captionstr += '\n'
captionstr += '    1D plot corresponding to the default parameters of the model.\n'
captionstr += '\n'

# Add figure reference and caption to documentation (at end, before References)
pattern = '\*\*REFERENCE'
m = re.search(pattern, docstr.upper())

if m:
    docstr1 = docstr[:m.start()]
    docstr2 = docstr[m.start():]
    docstr = docstr1 + captionstr + docstr2
else:
    print '------------------------------------------------------------------'
    print 'References NOT FOUND for model: ', model_info['id']
    print '------------------------------------------------------------------'
    docstr = docstr + captionstr

open(sys.argv[2],'w').write(docstr)
