"""
Mixture model
-------------

The product model multiplies the structure factor by the form factor,
modulated by the effective radius of the form.  The resulting model
has a attributes of both the model description (with parameters, etc.)
and the module evaluator (with call, release, etc.).

To use it, first load form factor P and structure factor S, then create
*ProductModel(P, S)*.
"""
from copy import copy
import numpy as np

from .modelinfo import Parameter, ParameterTable, ModelInfo
from .kernel import KernelModel, Kernel

try:
    from typing import List
    from .details import CallDetails
except ImportError:
    pass

def make_mixture_info(parts):
    # type: (List[ModelInfo]) -> ModelInfo
    """
    Create info block for product model.
    """
    flatten = []
    for part in parts:
        if part.composition and part.composition[0] == 'mixture':
            flatten.extend(part.composition[1])
        else:
            flatten.append(part)
    parts = flatten

    # Build new parameter list
    combined_pars = []
    for k, part in enumerate(parts):
        # Parameter prefix per model, A_, B_, ...
        # Note that prefix must also be applied to id and length_control
        # to support vector parameters
        prefix = chr(ord('A')+k) + '_'
        combined_pars.append(Parameter(prefix+'scale'))
        for p in part.parameters.kernel_parameters:
            p = copy(p)
            p.name = prefix + p.name
            p.id = prefix + p.id
            if p.length_control is not None:
                p.length_control = prefix + p.length_control
            combined_pars.append(p)
    parameters = ParameterTable(combined_pars)

    model_info = ModelInfo()
    model_info.id = '+'.join(part.id for part in parts)
    model_info.name = ' + '.join(part.name for part in parts)
    model_info.filename = None
    model_info.title = 'Mixture model with ' + model_info.name
    model_info.description = model_info.title
    model_info.docs = model_info.title
    model_info.category = "custom"
    model_info.parameters = parameters
    #model_info.single = any(part['single'] for part in parts)
    model_info.structure_factor = False
    model_info.variant_info = None
    #model_info.tests = []
    #model_info.source = []
    # Iq, Iqxy, form_volume, ER, VR and sesans
    # Remember the component info blocks so we can build the model
    model_info.composition = ('mixture', parts)


class MixtureModel(KernelModel):
    def __init__(self, model_info, parts):
        # type: (ModelInfo, List[KernelModel]) -> None
        self.info = model_info
        self.parts = parts

    def __call__(self, q_vectors):
        # type: (List[np.ndarray]) -> MixtureKernel
        # Note: may be sending the q_vectors to the n times even though they
        # are only needed once.  It would mess up modularity quite a bit to
        # handle this optimally, especially since there are many cases where
        # separate q vectors are needed (e.g., form in python and structure
        # in opencl; or both in opencl, but one in single precision and the
        # other in double precision).
        kernels = [part.make_kernel(q_vectors) for part in self.parts]
        return MixtureKernel(self.info, kernels)

    def release(self):
        # type: () -> None
        """
        Free resources associated with the model.
        """
        for part in self.parts:
            part.release()


class MixtureKernel(Kernel):
    def __init__(self, model_info, kernels):
        # type: (ModelInfo, List[Kernel]) -> None
        self.dim = kernels[0].dim
        self.info =  model_info
        self.kernels = kernels

    def __call__(self, call_details, value, weight, cutoff):
        # type: (CallDetails, np.ndarray, np.ndarry, float) -> np.ndarray
        scale, background = value[0:2]
        total = 0.0
        # remember the parts for plotting later
        self.results = []
        for kernel, kernel_details in zip(self.kernels, call_details.parts):
            part_result = kernel(kernel_details, value, weight, cutoff)
            total += part_result
            self.results.append(part_result)

        return scale*total + background

    def release(self):
        # type: () -> None
        for k in self.kernels:
            k.release()

