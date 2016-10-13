"""
Product model
-------------

The product model multiplies the structure factor by the form factor,
modulated by the effective radius of the form.  The resulting model
has a attributes of both the model description (with parameters, etc.)
and the module evaluator (with call, release, etc.).

To use it, first load form factor P and structure factor S, then create
*make_product_info(P, S)*.
"""
from __future__ import print_function, division

from copy import copy
import numpy as np  # type: ignore

from .modelinfo import Parameter, ParameterTable, ModelInfo
from .kernel import KernelModel, Kernel
from .details import make_details, dispersion_mesh

try:
    from typing import Tuple
except ImportError:
    pass

# TODO: make estimates available to constraints
#ESTIMATED_PARAMETERS = [
#    ["est_radius_effective", "A", 0.0, [0, np.inf], "", "Estimated effective radius"],
#    ["est_volume_ratio", "", 1.0, [0, np.inf], "", "Estimated volume ratio"],
#]

ER_ID = "radius_effective"
VF_ID = "volfraction"

# TODO: core_shell_sphere model has suppressed the volume ratio calculation
# revert it after making VR and ER available at run time as constraints.
def make_product_info(p_info, s_info):
    # type: (ModelInfo, ModelInfo) -> ModelInfo
    """
    Create info block for product model.
    """
    # Make sure effective radius is the first parameter and
    # make sure volume fraction is the second parameter of the
    # structure factor calculator.  Structure factors should not
    # have any magnetic parameters
    assert(s_info.parameters.kernel_parameters[0].id == ER_ID)
    assert(s_info.parameters.kernel_parameters[1].id == VF_ID)
    assert(s_info.parameters.magnetism_index == [])
    p_id, p_name, p_pars = p_info.id, p_info.name, p_info.parameters
    s_id, s_name, s_pars = s_info.id, s_info.name, s_info.parameters
    p_set = set(p.id for p in p_pars.call_parameters)
    s_set = set(p.id for p in s_pars.call_parameters)

    if p_set & s_set:
        # there is some overlap between the parameter names; tag the
        # overlapping S parameters with name_S.
        # Skip the first parameter of s, which is effective radius
        s_list = [(suffix_parameter(par) if par.id in p_set else par)
                  for par in s_pars.kernel_parameters[1:]]
    else:
        # Skip the first parameter of s, which is effective radius
        s_list = s_pars.kernel_parameters[1:]
    translate_name = dict((old.id, new.id) for old, new
                          in zip(s_pars.kernel_parameters[1:], s_list))
    demo = {}
    demo.update((k, v) for k, v in p_info.demo.items()
                if k not in ("background", "scale"))
    demo.update((translate_name[k], v) for k, v in s_info.demo.items()
                if k not in ("background", "scale") and not k.startswith(ER_ID))
    combined_pars = p_pars.kernel_parameters + s_list
    parameters = ParameterTable(combined_pars)
    parameters.max_pd = p_pars.max_pd + s_pars.max_pd

    model_info = ModelInfo()
    model_info.id = '*'.join((p_id, s_id))
    model_info.name = '*'.join((p_name, s_name))
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
    # Remember the component info blocks so we can build the model
    model_info.composition = ('product', [p_info, s_info])
    model_info.demo = {}
    return model_info

def suffix_parameter(par, suffix):
    par = copy(par)
    par.name = par.name + " S"
    par.id = par.id + "_S"
    return par

class ProductModel(KernelModel):
    def __init__(self, model_info, P, S):
        # type: (ModelInfo, KernelModel, KernelModel) -> None
        self.info = model_info
        self.P = P
        self.S = S

    def make_kernel(self, q_vectors):
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
        self.dtype = p_kernel.dtype
        self.results = []  # type: List[np.ndarray]

    def __call__(self, call_details, values, cutoff, magnetic):
        # type: (CallDetails, np.ndarray, float, bool) -> np.ndarray
        p_info, s_info = self.info.composition[1]

        # if there are magnetic parameters, they will only be on the
        # form factor P, not the structure factor S.
        nmagnetic = len(self.info.parameters.magnetism_index)
        if nmagnetic:
            spin_index = self.info.parameters.npars + 2
            magnetism = values[spin_index: spin_index+3+3*nmagnetic]
        else:
            magnetism = []
        nvalues = self.info.parameters.nvalues
        nweights = call_details.num_weights
        weights = values[nvalues:nvalues + 2*nweights]

        # Construct the calling parameters for P.
        p_npars = p_info.parameters.npars
        p_length = call_details.length[:p_npars]
        p_offset = call_details.offset[:p_npars]
        p_details = make_details(p_info, p_length, p_offset, nweights)
        # Set p scale to the volume fraction in s, which is the first of the
        # 'S' parameters in the parameter list, or 2+np in 0-origin.
        volfrac = values[2+p_npars]
        p_values = [[volfrac, 0.0], values[2:2+p_npars], magnetism, weights]
        spacer = (32 - sum(len(v) for v in p_values)%32)%32
        p_values.append([0.]*spacer)
        p_values = np.hstack(p_values).astype(self.p_kernel.dtype)

        # Call ER and VR for P since these are needed for S.
        p_er, p_vr = calc_er_vr(p_info, p_details, p_values)
        s_vr = (volfrac/p_vr if p_vr != 0. else volfrac)
        #print("volfrac:%g p_er:%g p_vr:%g s_vr:%g"%(volfrac,p_er,p_vr,s_vr))

        # Construct the calling parameters for S.
        # The  effective radius is not in the combined parameter list, so
        # the number of 'S' parameters is one less than expected.  The
        # computed effective radius needs to be added into the weights
        # vector, especially since it is a polydisperse parameter in the
        # stand-alone structure factor models.  We will added it at the
        # end so the remaining offsets don't need to change.
        s_npars = s_info.parameters.npars-1
        s_length = call_details.length[p_npars:p_npars+s_npars]
        s_offset = call_details.offset[p_npars:p_npars+s_npars]
        s_length = np.hstack((1, s_length))
        s_offset = np.hstack((nweights, s_offset))
        s_details = make_details(s_info, s_length, s_offset, nweights+1)
        v, w = weights[:nweights], weights[nweights:]
        s_values = [
            # scale=1, background=0, radius_effective=p_er, volfraction=s_vr
            [1., 0., p_er, s_vr],
            # structure factor parameters start after scale, background and
            # all the form factor parameters.  Skip the volfraction parameter
            # as well, since it is computed elsewhere, and go to the end of the
            # parameter list.
            values[2+p_npars+1:2+p_npars+s_npars],
            # no magnetism parameters to include for S
            # add er into the (value, weights) pairs
            v, [p_er], w, [1.0]
        ]
        spacer = (32 - sum(len(v) for v in s_values)%32)%32
        s_values.append([0.]*spacer)
        s_values = np.hstack(s_values).astype(self.s_kernel.dtype)

        # Call the kernels
        p_result = self.p_kernel(p_details, p_values, cutoff, magnetic)
        s_result = self.s_kernel(s_details, s_values, cutoff, False)

        #print("p_npars",p_npars,s_npars,p_er,s_vr,values[2+p_npars+1:2+p_npars+s_npars])
        #call_details.show(values)
        #print("values", values)
        #p_details.show(p_values)
        #print("=>", p_result)
        #s_details.show(s_values)
        #print("=>", s_result)

        # remember the parts for plotting later
        self.results = [p_result, s_result]

        #import pylab as plt
        #plt.subplot(211); plt.loglog(self.p_kernel.q_input.q, p_result, '-')
        #plt.subplot(212); plt.loglog(self.s_kernel.q_input.q, s_result, '-')
        #plt.figure()

        return values[0]*(p_result*s_result) + values[1]

    def release(self):
        # type: () -> None
        self.p_kernel.release()
        self.s_kernel.release()


def calc_er_vr(model_info, call_details, values):
    # type: (ModelInfo, ParameterSet) -> Tuple[float, float]

    if model_info.ER is None and model_info.VR is None:
        return 1.0, 1.0

    nvalues = model_info.parameters.nvalues
    value = values[nvalues:nvalues + call_details.num_weights]
    weight = values[nvalues + call_details.num_weights: nvalues + 2*call_details.num_weights]
    npars = model_info.parameters.npars
    pairs = [(value[offset:offset+length], weight[offset:offset+length])
             for p, offset, length
             in zip(model_info.parameters.call_parameters[2:2+npars],
                    call_details.offset,
                    call_details.length)
             if p.type == 'volume']
    value, weight = dispersion_mesh(model_info, pairs)

    if model_info.ER is not None:
        individual_radii = model_info.ER(*value)
        radius_effective = np.sum(weight*individual_radii) / np.sum(weight)
    else:
        radius_effective = 1.0

    if model_info.VR is not None:
        whole, part = model_info.VR(*value)
        volume_ratio = np.sum(weight*part)/np.sum(weight*whole)
    else:
        volume_ratio = 1.0

    return radius_effective, volume_ratio
