"""
Python driver for python kernels

Calls the kernel with a vector of $q$ values for a single parameter set.
Polydispersity is supported by looping over different parameter sets and
summing the results.  The interface to :class:`PyModel` matches those for
:class:`kernelcl.GpuModel` and :class:`kerneldll.DllModel`.
"""
from __future__ import division, print_function

import logging

import numpy as np  # type: ignore

from .generate import F64
from .kernel import KernelModel, Kernel

# pylint: disable=unused-import
try:
    from typing import Union, Callable
except ImportError:
    pass
else:
    from . import details
    DType = Union[None, str, np.dtype]
# pylint: enable=unused-import

logger = logging.getLogger(__name__)

class PyModel(KernelModel):
    """
    Wrapper for pure python models.
    """
    def __init__(self, model_info):
        # Make sure Iq is available and vectorized
        _create_default_functions(model_info)
        self.info = model_info
        self.dtype = np.dtype('d')
        logger.info("make python model " + self.info.name)

    def make_kernel(self, q_vectors):
        q_input = PyInput(q_vectors, dtype=F64)
        return PyKernel(self.info, q_input)

    def release(self):
        """
        Free resources associated with the model.
        """
        pass

class PyInput(object):
    """
    Make q data available to the gpu.

    *q_vectors* is a list of q vectors, which will be *[q]* for 1-D data,
    and *[qx, qy]* for 2-D data.  Internally, the vectors will be reallocated
    to get the best performance on OpenCL, which may involve shifting and
    stretching the array to better match the memory architecture.  Additional
    points will be evaluated with *q=1e-3*.

    *dtype* is the data type for the q vectors. The data type should be
    set to match that of the kernel, which is an attribute of
    :class:`GpuProgram`.  Note that not all kernels support double
    precision, so even if the program was created for double precision,
    the *GpuProgram.dtype* may be single precision.

    Call :meth:`release` when complete.  Even if not called directly, the
    buffer will be released when the data object is freed.
    """
    def __init__(self, q_vectors, dtype):
        self.nq = q_vectors[0].size
        self.dtype = dtype
        self.is_2d = (len(q_vectors) == 2)
        if self.is_2d:
            self.q = np.empty((self.nq, 2), dtype=dtype)
            self.q[:, 0] = q_vectors[0]
            self.q[:, 1] = q_vectors[1]
        else:
            self.q = np.empty(self.nq, dtype=dtype)
            self.q[:self.nq] = q_vectors[0]

    def release(self):
        """
        Free resources associated with the model inputs.
        """
        self.q = None

class PyKernel(Kernel):
    """
    Callable SAS kernel.

    *kernel* is the kernel object to call.

    *model_info* is the module information

    *q_input* is the DllInput q vectors at which the kernel should be
    evaluated.

    The resulting call method takes the *pars*, a list of values for
    the fixed parameters to the kernel, and *pd_pars*, a list of (value,weight)
    vectors for the polydisperse parameters.  *cutoff* determines the
    integration limits: any points with combined weight less than *cutoff*
    will not be calculated.

    Call :meth:`release` when done with the kernel instance.
    """
    def __init__(self, model_info, q_input):
        # type: (callable, ModelInfo, List[np.ndarray]) -> None
        self.dtype = np.dtype('d')
        self.info = model_info
        self.q_input = q_input
        self.res = np.empty(q_input.nq, q_input.dtype)
        self.dim = '2d' if q_input.is_2d else '1d'

        partable = model_info.parameters
        #kernel_parameters = (partable.iqxy_parameters if q_input.is_2d
        #                     else partable.iq_parameters)
        kernel_parameters = partable.iq_parameters
        volume_parameters = partable.form_volume_parameters

        # Create an array to hold the parameter values.  There will be a
        # single array whose values are updated as the calculator goes
        # through the loop.  Arguments to the kernel and volume functions
        # will use views into this vector, relying on the fact that a
        # an array of no dimensions acts like a scalar.
        parameter_vector = np.empty(len(partable.call_parameters)-2, 'd')

        # Create views into the array to hold the arguments
        offset = 0
        kernel_args, volume_args = [], []
        for p in partable.kernel_parameters:
            if p.length == 1:
                # Scalar values are length 1 vectors with no dimensions.
                v = parameter_vector[offset:offset+1].reshape(())
            else:
                # Vector values are simple views.
                v = parameter_vector[offset:offset+p.length]
            offset += p.length
            if p in kernel_parameters:
                kernel_args.append(v)
            if p in volume_parameters:
                volume_args.append(v)

        # Hold on to the parameter vector so we can use it to call kernel later.
        # This may also be required to preserve the views into the vector.
        self._parameter_vector = parameter_vector

        # Generate a closure which calls the kernel with the views into the
        # parameter array.
        if q_input.is_2d:
            form = model_info.Iqxy
            qx, qy = q_input.q[:, 0], q_input.q[:, 1]
            self._form = lambda: form(qx, qy, *kernel_args)
        else:
            form = model_info.Iq
            q = q_input.q
            self._form = lambda: form(q, *kernel_args)

        # Generate a closure which calls the form_volume if it exists.
        form_volume = model_info.form_volume
        self._volume = ((lambda: form_volume(*volume_args)) if form_volume else
                        (lambda: 1.0))

    def __call__(self, call_details, values, cutoff, magnetic):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool) -> np.ndarray
        if magnetic:
            raise NotImplementedError("Magnetism not implemented for pure python models")
        #print("Calling python kernel")
        #call_details.show(values)
        res = _loops(self._parameter_vector, self._form, self._volume,
                     self.q_input.nq, call_details, values, cutoff)
        return res

    def release(self):
        # type: () -> None
        """
        Free resources associated with the kernel.
        """
        self.q_input.release()
        self.q_input = None

def _loops(parameters,    # type: np.ndarray
           form,          # type: Callable[[], np.ndarray]
           form_volume,   # type: Callable[[], float]
           nq,            # type: int
           call_details,  # type: details.CallDetails
           values,        # type: np.ndarray
           cutoff         # type: float
          ):
    # type: (...) -> None
    ################################################################
    #                                                              #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #   !!                                                    !!   #
    #   !!  KEEP THIS CODE CONSISTENT WITH KERNEL_TEMPLATE.C  !!   #
    #   !!                                                    !!   #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #                                                              #
    ################################################################
    n_pars = len(parameters)
    parameters[:] = values[2:n_pars+2]
    if call_details.num_active == 0:
        pd_norm = float(form_volume())
        scale = values[0]/(pd_norm if pd_norm != 0.0 else 1.0)
        background = values[1]
        return scale*form() + background

    pd_value = values[2+n_pars:2+n_pars + call_details.num_weights]
    pd_weight = values[2+n_pars + call_details.num_weights:]

    pd_norm = 0.0
    partial_weight = np.NaN
    weight = np.NaN

    p0_par = call_details.pd_par[0]
    p0_length = call_details.pd_length[0]
    p0_index = p0_length
    p0_offset = call_details.pd_offset[0]

    pd_par = call_details.pd_par[:call_details.num_active]
    pd_offset = call_details.pd_offset[:call_details.num_active]
    pd_stride = call_details.pd_stride[:call_details.num_active]
    pd_length = call_details.pd_length[:call_details.num_active]

    total = np.zeros(nq, 'd')
    for loop_index in range(call_details.num_eval):
        # update polydispersity parameter values
        if p0_index == p0_length:
            pd_index = (loop_index//pd_stride)%pd_length
            parameters[pd_par] = pd_value[pd_offset+pd_index]
            partial_weight = np.prod(pd_weight[pd_offset+pd_index][1:])
            p0_index = loop_index%p0_length

        weight = partial_weight * pd_weight[p0_offset + p0_index]
        parameters[p0_par] = pd_value[p0_offset + p0_index]
        p0_index += 1
        if weight > cutoff:
            # Call the scattering function
            # Assume that NaNs are only generated if the parameters are bad;
            # exclude all q for that NaN.  Even better would be to have an
            # INVALID expression like the C models, but that is too expensive.
            Iq = np.asarray(form(), 'd')
            if np.isnan(Iq).any():
                continue

            # update value and norm
            total += weight * Iq
            pd_norm += weight * form_volume()

    scale = values[0]/(pd_norm if pd_norm != 0.0 else 1.0)
    background = values[1]
    return scale*total + background


def _create_default_functions(model_info):
    """
    Autogenerate missing functions, such as Iqxy from Iq.

    This only works for Iqxy when Iq is written in python. :func:`make_source`
    performs a similar role for Iq written in C.  This also vectorizes
    any functions that are not already marked as vectorized.
    """
    # Note: must call create_vector_Iq before create_vector_Iqxy
    _create_vector_Iq(model_info)
    _create_vector_Iqxy(model_info)


def _create_vector_Iq(model_info):
    """
    Define Iq as a vector function if it exists.
    """
    Iq = model_info.Iq
    if callable(Iq) and not getattr(Iq, 'vectorized', False):
        #print("vectorizing Iq")
        def vector_Iq(q, *args):
            """
            Vectorized 1D kernel.
            """
            return np.array([Iq(qi, *args) for qi in q])
        vector_Iq.vectorized = True
        model_info.Iq = vector_Iq


def _create_vector_Iqxy(model_info):
    """
    Define Iqxy as a vector function if it exists, or default it from Iq().
    """
    Iqxy = getattr(model_info, 'Iqxy', None)
    if callable(Iqxy):
        if not getattr(Iqxy, 'vectorized', False):
            #print("vectorizing Iqxy")
            def vector_Iqxy(qx, qy, *args):
                """
                Vectorized 2D kernel.
                """
                return np.array([Iqxy(qxi, qyi, *args) for qxi, qyi in zip(qx, qy)])
            vector_Iqxy.vectorized = True
            model_info.Iqxy = vector_Iqxy
    else:
        #print("defaulting Iqxy")
        # Iq is vectorized because create_vector_Iq was already called.
        Iq = model_info.Iq
        def default_Iqxy(qx, qy, *args):
            """
            Default 2D kernel.
            """
            return Iq(np.sqrt(qx**2 + qy**2), *args)
        default_Iqxy.vectorized = True
        model_info.Iqxy = default_Iqxy
