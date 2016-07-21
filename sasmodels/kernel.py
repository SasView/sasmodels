"""
Execution kernel interface
==========================

:class:`KernelModel` defines the interface to all kernel models.
In particular, each model should provide a :meth:`KernelModel.make_kernel`
call which returns an executable kernel, :class:`Kernel`, that operates
on the given set of *q_vector* inputs.  On completion of the computation,
the kernel should be released, which also releases the inputs.
"""

from __future__ import division, print_function

import numpy as np
from .details import mono_details, poly_details

try:
    from typing import List
except ImportError:
    pass
else:
    from .details import CallDetails
    from .modelinfo import ModelInfo
    import numpy as np  # type: ignore

class KernelModel(object):
    info = None  # type: ModelInfo
    dtype = None # type: np.dtype
    def make_kernel(self, q_vectors):
        # type: (List[np.ndarray]) -> "Kernel"
        raise NotImplementedError("need to implement make_kernel")

    def release(self):
        # type: () -> None
        pass

class Kernel(object):
    #: kernel dimension, either "1d" or "2d"
    dim = None  # type: str
    info = None  # type: ModelInfo
    results = None # type: List[np.ndarray]

    def __call__(self, call_details, values, cutoff, magnetic):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool) -> np.ndarray
        raise NotImplementedError("need to implement __call__")

    def release(self):
        # type: () -> None
        pass

try:
    np.meshgrid([])
    meshgrid = np.meshgrid
except ValueError:
    # CRUFT: np.meshgrid requires multiple vectors
    def meshgrid(*args):
        if len(args) > 1:
            return np.meshgrid(*args)
        else:
            return [np.asarray(v) for v in args]

def dispersion_mesh(model_info, pars):
    """
    Create a mesh grid of dispersion parameters and weights.

    Returns [p1,p2,...],w where pj is a vector of values for parameter j
    and w is a vector containing the products for weights for each
    parameter set in the vector.
    """
    value, weight = zip(*pars)
    weight = [w if w else [1.] for w in weight]
    weight = np.vstack([v.flatten() for v in meshgrid(*weight)])
    weight = np.prod(weight, axis=0)
    value = [v.flatten() for v in meshgrid(*value)]
    lengths = [par.length for par in model_info.parameters.kernel_parameters
               if par.type == 'volume']
    if any(n > 1 for n in lengths):
        pars = []
        offset = 0
        for n in lengths:
            pars.append(np.vstack(value[offset:offset+n]) if n > 1 else value[offset])
            offset += n
        value = pars
    return value, weight



def build_details(kernel, pairs):
    # type: (Kernel, Tuple[List[np.ndarray], List[np.ndarray]]) -> Tuple[CallDetails, np.ndarray, np.ndarray]
    """
    Construct the kernel call details object for calling the particular kernel.
    """
    values, weights = zip(*pairs)
    scalars = [v[0] for v in values]
    if all(len(w)==1 for w in weights):
        call_details = mono_details(kernel.info)
        data = np.array(scalars, dtype=kernel.dtype)
    else:
        call_details = poly_details(kernel.info, weights)
        data = np.hstack(scalars+list(values)+list(weights)).astype(kernel.dtype)
    return call_details, data
