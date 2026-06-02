from bumps.names import FitProblem, FreeVariables

from sasdata import data_path

from sasmodels.bumps_model import Experiment, Model
from sasmodels.core import load_model
from sasmodels.data import load_data

# latex data, same sample usans and sans
# particles radius ~2300, uniform dispersity
all_data = load_data(str(data_path / '1d_data' / 'latex_smeared.xml'), index='all')
#[print(data) for data in datasets]
datasets = all_data # Both SANS and USANS
#datasets = all_data[0:1] # SANS
#datasets = all_data[1:2] # USANS

# A single model to share between the datasets.  We will use
# FreeVariables below to set the parameters that are independent between
# the datasets.
model_name = "ellipsoid"
if model_name == "sphere":
    kernel = load_model("sphere")
    pars = dict(
        scale=0.01, background=0.0, sld=5.0, sld_solvent=0.0, radius=1500.,
        #radius_pd=0.1, radius_pd_n=35,
        )
elif model_name == "ellipsoid":
    kernel = load_model("ellipsoid")
    pars = dict(
        scale=0.01, background=0.0, sld=5.0, sld_solvent=0.0,
        radius_polar=1500., radius_equatorial=1500.,
        )
elif model_name == "triaxial_ellipsoid":
    kernel = load_model("triaxial_ellipsoid", dtype="double!")
    pars = dict(
        scale=0.01, background=0.0, sld=5.0, sld_solvent=0.0,
        radius_polar=1500., radius_equat_major=1500., radius_equat_minor=1500.,
        )
else:
    raise ValueError("Use model_name in sphere, ellipsoid, triaxial_ellipsoid")

model = Model(kernel, **pars)

# radius and polydispersity (if any) are shared
if model_name == "sphere":
    model.radius.range(0, 3000.)
    # model.radius_pd.range(0, 1)
elif model_name == "ellipsoid":
    model.radius_equatorial.range(0, 3000.)
    model.radius_polar.range(0, 3000.)
elif model_name == "triaxial_ellipsoid":
    model.radius_equat_major.range(0, 3000.)
    model.radius_equat_minor.range(0, 3000.)
    model.radius_polar.range(0, 3000.)


# Contrast and dilution are the same for both measurements, but are not
# separable with a single measurement (i.e., I(q) ~ F(q) contrast^2 Vf),
# so fit one of scale, sld or solvent sld.  With absolute scaling from
# data reduction, can use the same parameter for both datasets.
model.scale.range(0, 1.0)
#model.sld.range(-inf, inf)
#model.sld_solvent.range(-inf, inf)

# Background is different for sans and usans so set it as a free variable
# in the model.
free = FreeVariables(
    names=[data.run[0] for data in datasets],
    background=model.background,
    )
free.background.range(-1, 1)

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

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    problem.plot()
    plt.show()
