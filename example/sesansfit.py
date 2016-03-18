from bumps.names import *
from sasmodels import core, bumps_model, sesans

HAS_CONVERTER = True
try:
    from sas.sascalc.data_util.nxsunit import Converter
except ImportError:
    HAS_CONVERTER = False

def get_bumps_model(model_name):
    kernel = core.load_model(model_name)
    model = bumps_model.Model(kernel)
    return model

def sesans_fit(file, model, initial_vals={}, custom_params={}, param_range=[]):
    """
    @param file: SESANS file location
    @param model: Bumps model object or model name - can be model, model_1 * model_2, and/or model_1 + model_2
    @param initial_vals: dictionary of {param_name : initial_value}
    @param custom_params: dictionary of {custom_parameter_name : Parameter() object}
    @param param_range: dictionary of {parameter_name : [minimum, maximum]}
    @param constraints: dictionary of {parameter_name : constraint}
    @return: FitProblem for Bumps usage
    """
    try:
        from sas.sascalc.dataloader.loader import Loader
        loader = Loader()
        data = loader.load(file)
        if data is None: raise IOError("Could not load file %r"%(file))
        if HAS_CONVERTER == True:
            default_unit = "A"
            data_conv_q = Converter(data._xunit)
            for x in data.x:
                print x
            data.x = data_conv_q(data.x, units=default_unit)
            for x in data.x:
                print x
            data._xunit = default_unit

    except:
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
        data = SESANSData1D()

    if "radius" in initial_vals:
        radius = initial_vals.get("radius")
    else:
        radius = 1000
    data.Rmax = 3*radius # [A]

    if isinstance(model, basestring):
        model = get_bumps_model(model)

    # Load custom parameters, initial values and parameter ranges
    for k, v in custom_params.items():
        setattr(model, k, v)
        model._parameter_names.append(k)
    for k, v in initial_vals.items():
        param = model.parameters().get(k)
        setattr(param, "value", v)
    for k, v in param_range.items():
        param = model.parameters().get(k)
        if param is not None:
            setattr(param.bounds, "limits", v)

    if False: # for future implementation
        M_sesans = bumps_model.Experiment(data=data, model=model)
        M_sans = bumps_model.Experiment(data=sans_data, model=model)
        problem = FitProblem([M_sesans, M_sans])
    else:
        M_sesans = bumps_model.Experiment(data=data, model=model)
        problem = FitProblem(M_sesans)
    return problem