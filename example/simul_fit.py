from bumps.names import *
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasmodels.data import load_data, plot_data

# latex data, same sample usans and sans
# particles radius ~2300, uniform dispersity
datasets = load_data('latex_smeared.xml', index='all')
#[print(data) for data in datasets]

# A single sphere model to share between the datasets.  We will use
# FreeVariables below to set the parameters that are independent between
# the datasets.
kernel = load_model('sphere')
pars = dict(scale=0.01, background=0.0, sld=5.0, sld_solvent=0.0, radius=1500.,
            #radius_pd=0.1, radius_pd_n=35,
            )
model = Model(kernel, **pars)

# radius and polydispersity (if any) are shared
model.radius.range(0, inf)
#model.radius_pd.range(0, 1)

# Contrast and dilution are the same for both measurements, but are not
# separable with a single measurement (i.e., I(q) ~ F(q) contrast^2 Vf),
# so fit one of scale, sld or solvent sld.  With absolute scaling from
# data reduction, can use the same parameter for both datasets.
model.scale.range(0, inf)
#model.sld.range(-inf, inf)
#model.sld_solvent.range(-inf, inf)

# Background is different for sans and usans so set it as a free variable
# in the model.
free = FreeVariables(
    names=[data.run[0] for data in datasets],
    background=model.background,
    )
free.background.range(-inf, inf)

# Note: can access the parameters for the individual models using
# free.background[0] and free.background[1], setting constraints or
# ranges as appropriate.

# For more complex systems where different datasets require independent models,
# separate models can be defined, with parameters tied together using
# constraint expressions.  For example, the following could be used to fit
# data set 1 to spheres and data set 2 to cylinders of the same volume:
#    model1 = Model(load_model('sphere'))
#    model2 = Model(load_model('cylinder'))
#    model1.sld = model2.sld
#    model1.sld_solvent = model2.sld_solvent
#    model1.scale = model2.scale
#    # set cylinders and spheres to the same volume
#    model1.radius = (3/4*model2.radius**2*model2.length)**(1/3)
#    model1.background.range(0, 2)
#    model2.background.range(0, 2)

# Setup the experiments, sharing the same model across all datasets.
M = [Experiment(data=data, model=model, name=data.run[0]) for data in datasets]

problem = FitProblem(M, freevars=free)