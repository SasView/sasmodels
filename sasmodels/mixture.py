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
from __future__ import print_function

from copy import copy
import numpy as np  # type: ignore

from .modelinfo import Parameter, ParameterTable, ModelInfo
from .kernel import KernelModel, Kernel
from .details import make_details

try:
    from typing import List
except ImportError:
    pass

def make_mixture_info(parts):
    # type: (List[ModelInfo]) -> ModelInfo
    """
    Create info block for mixture model.
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
    demo = {}
    for k, part in enumerate(parts):
        # Parameter prefix per model, A_, B_, ...
        # Note that prefix must also be applied to id and length_control
        # to support vector parameters
        prefix = chr(ord('A')+k) + '_'
        scale =  Parameter(prefix+'scale', default=1.0,
                           description="model intensity for " + part.name)
        combined_pars.append(scale)
        for p in part.parameters.kernel_parameters:
            p = copy(p)
            p.name = prefix + p.name
            p.id = prefix + p.id
            if p.length_control is not None:
                p.length_control = prefix + p.length_control
            combined_pars.append(p)
        demo.update((prefix+k, v) for k, v in part.demo.items()
                    if k != "background")
    #print("pars",combined_pars)
    parameters = ParameterTable(combined_pars)
    parameters.max_pd = sum(part.parameters.max_pd for part in parts)

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
    model_info.demo = demo
    return model_info


class MixtureModel(KernelModel):
    def __init__(self, model_info, parts):
        # type: (ModelInfo, List[KernelModel]) -> None
        self.info = model_info
        self.parts = parts

    def make_kernel(self, q_vectors):
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
        self.dtype = self.kernels[0].dtype
        self.results = []  # type: List[np.ndarray]

    def __call__(self, call_details, values, cutoff, magnetic):
        # type: (CallDetails, np.ndarray, np.ndarry, float, bool) -> np.ndarray
        scale, background = values[0:2]
        total = 0.0
        # remember the parts for plotting later
        self.results = []  # type: List[np.ndarray]
        offset = 2 # skip scale & background
        parts = MixtureParts(self.info, self.kernels, call_details, values)
        for kernel, kernel_details, kernel_values in parts:
            #print("calling kernel", kernel.info.name)
            result = kernel(kernel_details, kernel_values, cutoff, magnetic)
            #print(kernel.info.name, result)
            total += result
            self.results.append(result)

        return scale*total + background

    def release(self):
        # type: () -> None
        for k in self.kernels:
            k.release()


class MixtureParts(object):
    def __init__(self, model_info, kernels, call_details, values):
        # type: (ModelInfo, List[Kernel], CallDetails, np.ndarray) -> None
        self.model_info = model_info
        self.parts = model_info.composition[1]
        self.kernels = kernels
        self.call_details = call_details
        self.values = values
        self.spin_index = model_info.parameters.npars + 2
        #call_details.show(values)

    def __iter__(self):
        # type: () -> PartIterable
        self.part_num = 0
        self.par_index = 2
        self.mag_index = self.spin_index + 3
        return self

    def next(self):
        # type: () -> Tuple[List[Callable], CallDetails, np.ndarray]
        if self.part_num >= len(self.parts):
            raise StopIteration()
        info = self.parts[self.part_num]
        kernel = self.kernels[self.part_num]
        call_details = self._part_details(info, self.par_index)
        values = self._part_values(info, self.par_index, self.mag_index)
        values = values.astype(kernel.dtype)
        #call_details.show(values)

        self.part_num += 1
        self.par_index += info.parameters.npars + 1
        self.mag_index += 3 * len(info.parameters.magnetism_index)

        return kernel, call_details, values

    def _part_details(self, info, par_index):
        # type: (ModelInfo, int) -> CallDetails
        full = self.call_details
        # par_index is index into values array of the current parameter,
        # which includes the initial scale and background parameters.
        # We want the index into the weight length/offset for each parameter.
        # Exclude the initial scale and background, so subtract two, but each
        # component has its own scale factor which we need to skip when
        # constructing the details for the kernel, so add one, giving a
        # net subtract one.
        index = slice(par_index - 1, par_index - 1 + info.parameters.npars)
        length = full.length[index]
        offset = full.offset[index]
        # The complete weight vector is being sent to each part so that
        # offsets don't need to be adjusted.
        part = make_details(info, length, offset, full.num_weights)
        return part

    def _part_values(self, info, par_index, mag_index):
        # type: (ModelInfo, int, int) -> np.ndarray
        #print(info.name, par_index, self.values[par_index:par_index + info.parameters.npars + 1])
        scale = self.values[par_index]
        pars = self.values[par_index + 1:par_index + info.parameters.npars + 1]
        nmagnetic = len(info.parameters.magnetism_index)
        if nmagnetic:
            spin_state = self.values[self.spin_index:self.spin_index + 3]
            mag_index = self.values[mag_index:mag_index + 3 * nmagnetic]
        else:
            spin_state = []
            mag_index = []
        nvalues = self.model_info.parameters.nvalues
        nweights = self.call_details.num_weights
        weights = self.values[nvalues:nvalues+2*nweights]
        zero = self.values.dtype.type(0.)
        values = [[scale, zero], pars, spin_state, mag_index, weights]
        # Pad value array to a 32 value boundary
        spacer = (32 - sum(len(v) for v in values)%32)%32
        values.append([zero]*spacer)
        values = np.hstack(values).astype(self.kernels[0].dtype)
        return values
