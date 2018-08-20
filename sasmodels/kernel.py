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

    def Iq(self, call_details, values, cutoff, magnetic):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool) -> np.ndarray
        Pq, Reff = self.Pq_Reff(call_details, values, cutoff, magnetic, effective_radius_type=0)
        scale = values[0]
        background = values[1]
        return scale*Pq + background
    __call__ = Iq

    def Pq_Reff(self, call_details, values, cutoff, magnetic, effective_radius_type):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool, int) -> np.ndarray
        self._call_kernel(call_details, values, cutoff, magnetic, effective_radius_type)
        #print("returned",self.q_input.q, self.result)
        nout = 2 if self.info.have_Fq else 1
        total_weight = self.result[nout*self.q_input.nq + 0]
        if total_weight == 0.:
            total_weight = 1.
        weighted_volume = self.result[nout*self.q_input.nq + 1]
        weighted_radius = self.result[nout*self.q_input.nq + 2]
        effective_radius = weighted_radius/total_weight
        # compute I = scale*P + background
        #           = scale*(sum(w*F^2)/sum w)/(sum (w*V)/sum w) + background
        #           = scale/sum (w*V) * sum(w*F^2) + background
        scale = values[0]/(weighted_volume if weighted_volume != 0.0 else 1.0)
        background = values[1]
        #print("scale",scale,background)
        return self.result[0:nout*self.q_input.nq:nout], effective_radius

    def beta(self, call_details, values, cutoff, magnetic, effective_radius_type):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool, int) -> np.ndarray
        if self.dim == '2d':
            raise NotImplementedError("beta not yet supported for 2D")
        if not self.info.have_Fq:
            raise NotImplementedError("beta not yet supported for "+self.info.id)
        self._call_kernel(call_details, values, cutoff, magnetic, effective_radius_type)
        total_weight = self.result[2*self.q_input.nq + 0]
        if total_weight == 0.:
            total_weight = 1.
        weighted_volume = self.result[2*self.q_input.nq + 1]
        weighted_radius = self.result[2*self.q_input.nq + 2]
        volume_average = weighted_volume/total_weight
        effective_radius = weighted_radius/total_weight
        F2 = self.result[0:2*self.q_input.nq:2]/total_weight
        F1 = self.result[1:2*self.q_input.nq:2]/total_weight
        return F1, F2, volume_average, effective_radius

    def release(self):
        # type: () -> None
        pass

    def _call_kernel(self, call_details, values, cutoff, magnetic, effective_radius_type):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool, int) -> np.ndarray
        raise NotImpmentedError()