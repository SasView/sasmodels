"""
Execution kernel interface
==========================

:class:`KernelModel` defines the interface to all kernel models.
In particular, each model should provide a :meth:`KernelModel.make_kernel`
call which returns an executable kernel, :class:`Kernel`, that operates
on the given set of *q_vector* inputs.  On completion of the computation,
the kernel should be released, which also releases the inputs.
"""

try:
    from typing import List
    from .details import CallDetails
    from .modelinfo import ModelInfo
    import numpy as np  # type: ignore
except ImportError:
    pass

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

    def __call__(self, call_details, weights, values, cutoff):
        # type: (CallDetails, np.ndarray, np.ndarray, float) -> np.ndarray
        raise NotImplementedError("need to implement __call__")

    def release(self):
        # type: () -> None
        pass
