import logging

from bumps.names import *
from sasmodels import core, bumps_model, sesans
from sasdata.dataloader.loader import Loader

def get_bumps_model(model_name):
    kernel = core.load_model(model_name)
    model = bumps_model.Model(kernel)
    return model

def sesans_fit(file, model, initial_vals={}, custom_params={}, param_range=[],
               acceptance_angle=None):
    """

    @param file: SESANS file location
    @param model: Bumps model object or model name - can be model, model_1 * model_2, and/or model_1 + model_2
    @param initial_vals: dictionary of {param_name : initial_value}
    @param custom_params: dictionary of {custom_parameter_name : Parameter() object}
    @param param_range: dictionary of {parameter_name : [minimum, maximum]}
    @param constraints: dictionary of {parameter_name : constraint}
    @return: FitProblem for Bumps usage
    """
    logging.basicConfig()

    initial_vals['background'] = 0.0
    try:
        loader = Loader()
        data = loader.load(file)[0]
        if data is None:
            raise IOError("Could not load file %r"%(file))

    except Exception:
        raise
        # If no loadable data file, generate random data
        SElength = np.linspace(0, 2400, 61) # [A]
        data = np.ones_like(SElength)
        err_data = np.ones_like(SElength)*0.03

        class Sample:
            zacceptance = 0.1 # [A^-1]
            thickness = 0.2 # [cm]

        class SESANSData1D:
            #q_zmax = 0.23 # [A^-1]
            lam = 0.2 # [nm]
            x = SElength
            y = data
            dy = err_data
            sample = Sample()
            has_no_finite_acceptance = acceptance_angle is not None
        data = SESANSData1D()
        data.acceptance_angle = acceptance_angle

    data.has_no_finite_acceptance = acceptance_angle is not None
    radius = initial_vals.get("radius", 1000)
    data.Rmax = 30*radius # [A]

    if isinstance(model, str):
        model = get_bumps_model(model)

    # Load custom parameters, initial values and parameter ranges
    pars = model.parameters().copy()
    pars.update(custom_params)
    for k, v in custom_params.items():
        if k not in pars:
            par_error("Can't set parameter %r in model"%k, pars)
        setattr(model, k, v)
        model._parameter_names.append(k)
    for k, v in initial_vals.items():
        if k not in pars:
            par_error("Can't get parameter %r from model"%k, pars)
        param = pars[k]
        param.value = v
    for k, v in param_range.items():
        if k not in pars:
            par_error("Can't set range on parameter %r in model"%k, pars)
        param = pars[k]
        param.range(*v)

    if False: # for future implementation
        M_sesans = bumps_model.Experiment(data=data, model=model)
        M_sans = bumps_model.Experiment(data=sans_data, model=model)
        problem = FitProblem([M_sesans, M_sans])
    else:
        M_sesans = bumps_model.Experiment(data=data, model=model)
        problem = FitProblem(M_sesans)
    return problem


def par_error(msg, pars):
    raise ValueError(msg+"\nAvailable parameters: %s"%", ".join(sorted(pars.keys())))
