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
from collections import OrderedDict

import numpy as np  # type: ignore

from .modelinfo import Parameter, ParameterTable, ModelInfo
from .modelinfo import NUM_MAGFIELD_PARS, NUM_MAGNETIC_PARS, NUM_COMMON_PARS
from .kernel import KernelModel, Kernel
from .details import make_details

# pylint: disable=unused-import
try:
    from typing import List, Tuple, Optional, Callable, Any
    from .details import CallDetails
except ImportError:
    pass
# pylint: enable=unused-import

def make_mixture_info(parts, operation='+'):
    # type: (List[ModelInfo], str) -> ModelInfo
    """
    Create info block for mixture model.
    """
    # Build new parameter list
    combined_pars = []

    # When creating a mixture model that is a sum of product models (ie (1*2)+(3*4))
    # the parameters for models 1 & 2 will be prefixed with A & B respectively,
    # but so will the parameters for models 3 & 4. We need to rename models 3 & 4
    # so that they are prefixed with C & D to avoid overlap of parameter names.
    used_prefixes = []
    for part in parts:
        if part.composition and part.composition[0] == 'mixture':
            i = 0
            for submodel in part.composition[1]:
                npars = len(submodel.parameters.kernel_parameters)
                # List of params of one of the constituent models of part
                submodel_pars = part.parameters.kernel_parameters[i:i+npars]
                # Prefix of the constituent model
                prefix = submodel_pars[0].name[0]
                if prefix not in used_prefixes: # Haven't seen this prefix so far
                    used_prefixes.append(prefix)
                    i += npars
                    continue
                # TODO: don't modify submodel --- it may be used elsewhere
                # Existing code probably doesn't keep a handle on the model
                # parts so its probably okay, but it's possible that a mix
                # on user defined mixture models models will change the
                # parameters used for the parts in the GUI. Even worse if the
                # same plugin is used twice. For example, twosphere.py
                # contains sphere+sphere and you create twosphere+twosphere.
                while prefix in used_prefixes:
                    # This prefix has been already used, so change it to the
                    # next letter that hasn't been used
                    prefix = chr(ord(prefix) + 1)
                used_prefixes.append(prefix)
                prefix += "_"
                # Update the parameters of this constituent model to use the
                # new prefix
                for par in submodel_pars:
                    # Strip {prefix}_ using par.name[2:], etc.
                    # TODO: fails for AB_scale
                    par.id = prefix + par.id[2:]
                    par.name = prefix + par.name[2:]
                    if par.length_control is not None:
                        par.length_control = prefix + par.length_control[2:]
                i += npars

    for part in parts:
        # Parameter prefix per model, A_, B_, ...
        # Note that prefix must also be applied to id and length_control
        # to support vector parameters
        prefix = ''
        if not part.composition or part.composition[0] == 'product':
            # Model isn't a composition model, so its parameters don't have a
            # a prefix. Add the next available prefix
            prefix = chr(ord('A')+len(used_prefixes))
            used_prefixes.append(prefix)
            prefix += '_'

        if operation == '+':
            # If model is a sum model, each constituent model gets its own scale parameter
            scale_prefix = prefix
            if prefix == '' and getattr(part, "operation", '') == '*':
                # `part` is a composition product model. Find the prefixes of
                # its parameters to form a new prefix for the scale.
                # For example, a model with A*B*C will have ABC_scale.
                sub_prefixes = []
                for param in part.parameters.kernel_parameters:
                    # Prefix of constituent model
                    sub_prefix = param.id.split('_')[0]
                    if sub_prefix not in sub_prefixes:
                        sub_prefixes.append(sub_prefix)
                # Concatenate sub_prefixes to form prefix for the scale
                scale_prefix = ''.join(sub_prefixes) + '_'
            scale = Parameter(scale_prefix + 'scale', default=1.0,
                              description="model intensity for " + part.name)
            combined_pars.append(scale)
        for p in part.parameters.kernel_parameters:
            p = copy(p)
            p.name = prefix + p.name
            p.id = prefix + p.id
            if p.length_control is not None:
                p.length_control = prefix + p.length_control
            combined_pars.append(p)
    parameters = ParameterTable(combined_pars)
    # Allow for the scenario in which each component has all its PD parameters
    # active simultaneously.  details.make_details() will throw an error if
    # too many are used from any one component.
    parameters.max_pd = sum(part.parameters.max_pd for part in parts)

    def random():
        """Random set of model parameters for mixture model"""
        combined_pars = {}
        for k, part in enumerate(parts):
            prefix = chr(ord('A')+k) + '_'
            pars = part.random()
            combined_pars.update((prefix+k, v) for k, v in pars.items())
        return combined_pars

    model_info = ModelInfo()
    model_info.id = operation.join(part.id for part in parts)
    model_info.operation = operation
    model_info.name = '(' + operation.join(part.name for part in parts) + ')'
    model_info.filename = None
    model_info.title = 'Mixture model with ' + model_info.name
    model_info.description = model_info.title
    model_info.docs = model_info.title
    model_info.category = "custom"
    model_info.parameters = parameters
    model_info.random = random
    #model_info.single = any(part['single'] for part in parts)
    model_info.structure_factor = False
    #model_info.tests = []
    #model_info.source = []
    # Remember the component info blocks so we can build the model
    model_info.composition = ('mixture', parts)
    return model_info


class MixtureModel(KernelModel):
    """
    Model definition for mixture of models.
    """
    def __init__(self, model_info, parts):
        # type: (ModelInfo, List[KernelModel]) -> None
        self.info = model_info
        self.parts = parts
        self.dtype = parts[0].dtype

    def make_kernel(self, q_vectors):
        # type: (List[np.ndarray]) -> MixtureKernel
        # Note: may be sending the q_vectors to the n times even though they
        # are only needed once.  It would mess up modularity quite a bit to
        # handle this optimally, especially since there are many cases where
        # separate q vectors are needed (e.g., form in python and structure
        # in opencl; or both in opencl, but one in single precision and the
        # other in double precision).
        kernels = [part.make_kernel(q_vectors) for part in self.parts]
        return MixtureKernel(self.info, kernels, q_vectors)
    make_kernel.__doc__ = KernelModel.make_kernel.__doc__

    def release(self):
        # type: () -> None
        """Free resources associated with the model."""
        for part in self.parts:
            part.release()
    release.__doc__ = KernelModel.release.__doc__


def _intermediates(q, results):
    # type: (np.ndarray, List[Tuple[Kernel, np.ndarray, Optional[Callable]]]) -> OrderedDict[str, Any]
    """
    Returns intermediate results for mixture model.
    """
    parts = OrderedDict()
    for part_num, (kernel, Iq, intermediates) in enumerate(results):
        part_id = chr(ord('A') + part_num) + "_" + kernel.info.name + "(Q)"
        parts[part_id] = (q, Iq)
        if intermediates is not None:
            parts[part_id + " parts"] = intermediates()
    return parts


class MixtureKernel(Kernel):
    """
    Instantiated kernel for mixture of models.
    """
    def __init__(self, model_info, kernels, q):
        # type: (ModelInfo, List[Kernel], Tuple[np.ndarray]) -> None
        self.dim = kernels[0].dim
        self.info = model_info
        self.q = q
        self.kernels = kernels
        self.dtype = self.kernels[0].dtype
        self.operation = model_info.operation
        self.results = None  # type: Callable[[], OrderedDict]

    def Iq(self, call_details, values, cutoff, magnetic):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool) -> np.ndarray
        scale, background = values[0:2]
        total = 0.0
        # remember the parts for plotting later
        results = []
        parts = _MixtureParts(self.info, self.kernels, call_details, values)
        for kernel, kernel_details, kernel_values in parts:
            #print("calling kernel", kernel.info.name)
            result = kernel(kernel_details, kernel_values, cutoff, magnetic)
            result = np.array(result).astype(kernel.dtype)
            # print(kernel.info.name, result)
            if self.operation == '+':
                total += result
            elif self.operation == '*':
                if np.all(total) == 0.0:
                    total = result
                else:
                    total *= result
            results.append((kernel, result, getattr(kernel, 'results', None)))

        self.results = lambda: _intermediates(self.q, results)

        return scale*total + background

    Iq.__doc__ = Kernel.Iq.__doc__
    __call__ = Iq

    def release(self):
        # type: () -> None
        """Free resources associated with the kernel."""
        for k in self.kernels:
            k.release()


# Note: _MixtureParts doesn't implement iteration correctly, and only allows
# a single iterator to be active at once.  It doesn't matter in this case
# since _MixtureParts is only used in one place, but it is not clean style.
class _MixtureParts(object):
    """
    Mixture component iterator.
    """
    def __init__(self, model_info, kernels, call_details, values):
        # type: (ModelInfo, List[Kernel], CallDetails, np.ndarray) -> None
        self.model_info = model_info
        self.parts = model_info.composition[1]
        self.kernels = kernels
        self.call_details = call_details
        self.values = values
        self.spin_index = model_info.parameters.npars + NUM_COMMON_PARS
        # The following are redefined by __iter__, but set them here so that
        # lint complains a little less.
        self.part_num = -1
        self.par_index = -1
        self.mag_index = -1
        #call_details.show(values)

    def __iter__(self):
        # type: () -> _MixtureParts
        self.part_num = 0
        self.par_index = NUM_COMMON_PARS
        self.mag_index = self.spin_index + NUM_MAGFIELD_PARS
        return self

    def __next__(self):
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
        self.par_index += info.parameters.npars
        if self.model_info.operation == '+':
            self.par_index += 1 # Account for each constituent model's scale param
        self.mag_index += NUM_MAGNETIC_PARS*len(info.parameters.magnetism_index)

        return kernel, call_details, values

    # CRUFT: py2 support
    next = __next__

    def _part_details(self, info, par_index):
        # type: (ModelInfo, int) -> CallDetails
        full = self.call_details
        # par_index is index into values array of the current parameter,
        # which includes the initial scale and background parameters.
        # We want the index into the weight length/offset for each parameter.
        # Exclude the initial scale and background, so subtract two. If we're
        # building an addition model, each component has its own scale factor
        # which we need to skip when constructing the details for the kernel, so
        # add one, giving a net subtract one.
        diff = NUM_COMMON_PARS-1 if self.model_info.operation == '+' else NUM_COMMON_PARS
        index = slice(par_index - diff, par_index - diff + info.parameters.npars)
        length = full.length[index]
        offset = full.offset[index]
        # The complete weight vector is being sent to each part so that
        # offsets don't need to be adjusted.
        part = make_details(info, length, offset, full.num_weights)
        return part

    def _part_values(self, info, par_index, mag_index):
        # type: (ModelInfo, int, int) -> np.ndarray
        # Set each constituent model's scale to 1 if this is a multiplication model
        scale = self.values[par_index] if self.model_info.operation == '+' else 1.0
        diff = 1 if self.model_info.operation == '+' else 0 # Skip scale if addition model
        pars = self.values[par_index + diff:par_index + info.parameters.npars + diff]
        nmagnetic = len(info.parameters.magnetism_index)
        if nmagnetic:
            spin_state = self.values[self.spin_index:self.spin_index + NUM_MAGFIELD_PARS]
            mag_index = self.values[mag_index:mag_index + NUM_MAGNETIC_PARS*nmagnetic]
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
