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

# pylint: disable=unused-import
try:
    from typing import List
except ImportError:
    pass
else:
    import numpy as np
    from .details import CallDetails
    from .modelinfo import ModelInfo
# pylint: enable=unused-import

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
    dtype = None  # type: np.dtype

    def __call__(self, call_details, values, cutoff, magnetic):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool) -> np.ndarray
        raise NotImplementedError("need to implement __call__")

    def release(self):
        # type: () -> None
        pass
