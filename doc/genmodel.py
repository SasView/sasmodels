import sys, os, math, re
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath('..'))
from sasmodels import generate, core
from sasmodels.direct_model import DirectModel
from sasmodels.data import empty_data1D, empty_data2D


# Convert ../sasmodels/models/name.py to name
model_name = os.path.basename(sys.argv[1])[:-3]

# Load the doc string from the module definition file and store it in rst
docstr = generate.make_doc(core.load_model_info(model_name))
    
# Generate automatically plot of the model and add it to rst documentation

info = core.load_model_info(model_name)

# Calculate 1D curve for default parameters
pars = dict((p[0], p[2]) for p in info['parameters'])

# Plotting ranges and options
opts = {
        'xscale'    : 'log',
        'yscale'    : 'log' if not info['structure_factor'] else 'linear',
        'qmin'      : 0.001,
        'qmax'      : 1.0,
        'nq'        : 1000,
        'nq2d'      : 100,
}

qmin, qmax, nq = opts['qmin'], opts['qmax'], opts['nq']
qmin = math.log10(qmin)
qmax = math.log10(qmax)
q = np.logspace(qmin, qmax, nq)
data = empty_data1D(q)
model = core.load_model(model_name)
calculator = DirectModel(data, model)
Iq1D = calculator()

# TO DO: Generation of 2D plots 
# Problem in sasmodels.direct_model._calc_theory
# There self._kernel.q_input.nq  gets a value of 0 in the 2D case
# and returns a 0 numpy array (it does not call the C code)

# If 2D model, compute 2D image
#if info['has_2d'] != []:
#    qmax, nq2d = opts['qmax'], opts['nq2d']
#    data2d = empty_data2D(np.linspace(-qmax, qmax, nq2d), resolution=0.0) 
#    #model = core.load_model(model_name)
#    calculator = DirectModel(data2d, model)
#    Iq2D = calculator()

# Generate image (comment IF for 1D/2D for the moment) and generate only 1D
#if info['has_2d'] == []:
#    fig = plt.figure()
#    ax = fig.add_subplot(1,1,1)
#    ax.plot(q, Iq1D, color='blue', lw=2, label=model_name)
#    ax.set_xlabel(r'$Q \/(\AA^{-1})$')
#    ax.set_xscale(opts['xscale'])   
#    ax.set_ylabel(r'$I(Q) \/(\mathrm{cm}^{-1})$')
#    ax.set_yscale(opts['yscale'])  
#    ax.legend()
#else:
#    # need figure with 1D + 2D
#    pass
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(q, Iq1D, color='blue', lw=2, label=info['name'])
ax.set_xlabel(r'$Q \/(\AA^{-1})$')
ax.set_xscale(opts['xscale'])   
ax.set_ylabel(r'$I(Q) \/(\mathrm{cm}^{-1})$')
ax.set_yscale(opts['yscale'])  
ax.legend()
 

# Save image in model/img
figname = model_name + '_autogenfig.png'
filename = os.path.join('model', 'img', figname)
plt.savefig(filename)

# Auto caption for figure
captionstr = '\n'
captionstr += '.. figure:: img/' + model_name + '_autogenfig.png\n'
captionstr += '\n'
#if info['has_2d'] == []:
#    captionstr += '    1D plot corresponding to the default parameters of the model.\n'
#else:
#    captionstr += '    1D and 2D plots corresponding to the default parameters of the model.\n'
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
    print 'References NOT FOUND for model: ', model_name
    print '------------------------------------------------------------------'
    docstr = docstr + captionstr

open(sys.argv[2],'w').write(docstr)
