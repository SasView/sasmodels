"""
Product model
-------------

The product model multiplies the structure factor by the form factor,
modulated by the effective radius of the form.  The resulting model
has a attributes of both the model description (with parameters, etc.)
and the module evaluator (with call, release, etc.).

To use it, first load form factor P and structure factor S, then create
*ProductModel(P, S)*.
"""
import numpy as np

from .details import dispersion_mesh
from .modelinfo import suffix_parameter, ParameterTable, ModelInfo
from .kernel import KernelModel, Kernel

try:
    from typing import Tuple
    from .modelinfo import ParameterSet
    from .details import CallDetails
except ImportError:
    pass

# TODO: make estimates available to constraints
#ESTIMATED_PARAMETERS = [
#    ["est_effect_radius", "A", 0.0, [0, np.inf], "", "Estimated effective radius"],
#    ["est_volume_ratio", "", 1.0, [0, np.inf], "", "Estimated volume ratio"],
#]

# TODO: core_shell_sphere model has suppressed the volume ratio calculation
# revert it after making VR and ER available at run time as constraints.
def make_product_info(p_info, s_info):
    # type: (ModelInfo, ModelInfo) -> ModelInfo
    """
    Create info block for product model.
    """
    p_id, p_name, p_pars = p_info.id, p_info.name, p_info.parameters
    s_id, s_name, s_pars = s_info.id, s_info.name, s_info.parameters
    p_set = set(p.id for p in p_pars.call_parameters)
    s_set = set(p.id for p in s_pars.call_parameters)

    if p_set & s_set:
        # there is some overlap between the parameter names; tag the
        # overlapping S parameters with name_S
        s_list = [(suffix_parameter(par, "_S") if par.id in p_set else par)
                  for par in s_pars.kernel_parameters]
        combined_pars = p_pars.kernel_parameters + s_list
    else:
        combined_pars = p_pars.kernel_parameters + s_pars.kernel_parameters
    parameters = ParameterTable(combined_pars)

    model_info = ModelInfo()
    model_info.id = '*'.join((p_id, s_id))
    model_info.name = ' X '.join((p_name, s_name))
    model_info.filename = None
    model_info.title = 'Product of %s and %s'%(p_name, s_name)
    model_info.description = model_info.title
    model_info.docs = model_info.title
    model_info.category = "custom"
    model_info.parameters = parameters
    #model_info.single = p_info.single and s_info.single
    model_info.structure_factor = False
    model_info.variant_info = None
    #model_info.tests = []
    #model_info.source = []
    # Iq, Iqxy, form_volume, ER, VR and sesans
    model_info.composition = ('product', [p_info, s_info])
    return model_info

class ProductModel(KernelModel):
    def __init__(self, model_info, P, S):
        # type: (ModelInfo, KernelModel, KernelModel) -> None
        self.info = model_info
        self.P = P
        self.S = S

    def __call__(self, q_vectors):
        # type: (List[np.ndarray]) -> Kernel
        # Note: may be sending the q_vectors to the GPU twice even though they
        # are only needed once.  It would mess up modularity quite a bit to
        # handle this optimally, especially since there are many cases where
        # separate q vectors are needed (e.g., form in python and structure
        # in opencl; or both in opencl, but one in single precision and the
        # other in double precision).
        p_kernel = self.P.make_kernel(q_vectors)
        s_kernel = self.S.make_kernel(q_vectors)
        return ProductKernel(self.info, p_kernel, s_kernel)

    def release(self):
        # type: (None) -> None
        """
        Free resources associated with the model.
        """
        self.P.release()
        self.S.release()


class ProductKernel(Kernel):
    def __init__(self, model_info, p_kernel, s_kernel):
        # type: (ModelInfo, Kernel, Kernel) -> None
        self.info = model_info
        self.p_kernel = p_kernel
        self.s_kernel = s_kernel

    def __call__(self, details, weights, values, cutoff):
        # type: (CallDetails, np.ndarray, np.ndarray, float) -> np.ndarray
        effect_radius, vol_ratio = call_ER_VR(self.p_kernel.info, vol_pars)

        p_fixed[SCALE] = s_volfraction
        p_fixed[BACKGROUND] = 0.0
        s_fixed[SCALE] = scale
        s_fixed[BACKGROUND] = 0.0
        s_fixed[2] = s_volfraction/vol_ratio
        s_pd[0] = [effect_radius], [1.0]

        p_res = self.p_kernel(p_details, p_weights, p_values, cutoff)
        s_res = self.s_kernel(s_details, s_weights, s_values, cutoff)
        #print s_fixed, s_pd, p_fixed, p_pd

        return p_res*s_res + background

    def release(self):
        # type: () -> None
        self.p_kernel.release()
        self.s_kernel.release()

def call_ER_VR(model_info, pars):
    """
    Return effect radius and volume ratio for the model.
    """
    if model_info.ER is None and model_info.VR is None:
        return 1.0, 1.0

    value, weight = _vol_pars(model_info, pars)

    if model_info.ER is not None:
        individual_radii = model_info.ER(*value)
        effect_radius = np.sum(weight*individual_radii) / np.sum(weight)
    else:
        effect_radius = 1.0

    if model_info.VR is not None:
        whole, part = model_info.VR(*value)
        volume_ratio = np.sum(weight*part)/np.sum(weight*whole)
    else:
        volume_ratio = 1.0

    return effect_radius, volume_ratio

def _vol_pars(model_info, pars):
    # type: (ModelInfo, ParameterSet) -> Tuple[np.ndarray, np.ndarray]
    vol_pars = [get_weights(p, pars)
                for p in model_info.parameters.call_parameters
                if p.type == 'volume']
    value, weight = dispersion_mesh(model_info, vol_pars)
    return value, weight

