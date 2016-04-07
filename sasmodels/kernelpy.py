"""
Python driver for python kernels

Calls the kernel with a vector of $q$ values for a single parameter set.
Polydispersity is supported by looping over different parameter sets and
summing the results.  The interface to :class:`PyModel` matches those for
:class:`kernelcl.GpuModel` and :class:`kerneldll.DllModel`.
"""
import numpy as np
from numpy import pi, cos

from .generate import F64

class PyModel(object):
    """
    Wrapper for pure python models.
    """
    def __init__(self, model_info):
        self.info = model_info

    def make_kernel(self, q_vectors):
        q_input = PyInput(q_vectors, dtype=F64)
        kernel = self.info['Iqxy'] if q_input.is_2d else self.info['Iq']
        return PyKernel(kernel, self.info, q_input)

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

class PyKernel(object):
    """
    Callable SAS kernel.

    *kernel* is the DllKernel object to call.

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
    def __init__(self, kernel, model_info, q_input):
        self.dtype = np.dtype('d')
        self.info = model_info
        self.q_input = q_input
        self.res = np.empty(q_input.nq, q_input.dtype)
        self.kernel = kernel
        self.dim = '2d' if q_input.is_2d else '1d'

        partable = model_info['parameters']
        kernel_parameters = (partable.iqxy_parameters if q_input.is_2d
                             else partable.iq_parameters)
        volume_parameters = partable.form_volume_parameters

        # Create an array to hold the parameter values.  There will be a
        # single array whose values are updated as the calculator goes
        # through the loop.  Arguments to the kernel and volume functions
        # will use views into this vector, relying on the fact that a
        # an array of no dimensions acts like a scalar.
        parameter_vector = np.empty(len(partable.call_parameters), 'd')

        # Create views into the array to hold the arguments
        offset = 2
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
                volume_args.append(p)

        # Hold on to the parameter vector so we can use it to call kernel later.
        # This may also be required to preserve the views into the vector.
        self._parameter_vector = parameter_vector

        # Generate a closure which calls the kernel with the views into the
        # parameter array.
        if q_input.is_2d:
            form = model_info['Iqxy']
            qx, qy = q_input.q[:,0], q_input.q[:,1]
            self._form = lambda: form(qx, qy, *kernel_args)
        else:
            form = model_info['Iq']
            q = q_input.q
            self._form = lambda: form(q, *kernel_args)

        # Generate a closure which calls the form_volume if it exists.
        form_volume = model_info['form_volume']
        self._volume = ((lambda: form_volume(*volume_args)) if form_volume
                        else (lambda: 1.0))

    def __call__(self, details, weights, values, cutoff):
        # type: (.generate.CoordinationDetails, np.ndarray, np.ndarray, float) -> np.ndarray
        res = _loops(self._parameter_vector, self._form, self._volume,
                     self.q_input.nq, details, weights, values, cutoff)
        return res

    def release(self):
        """
        Free resources associated with the kernel.
        """
        self.q_input = None

def _loops(parameters,  # type: np.ndarray
           form,        # type: Callable[[], np.ndarray]
           form_volume, # type: Callable[[], float]
           nq,          # type: int
           details,     # type: .generate.CoordinationDetails
           weights,     # type: np.ndarray
           values,      # type: np.ndarray
           cutoff,      # type: float
           ):           # type: (...) -> None
    ################################################################
    #                                                              #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #   !!                                                    !!   #
    #   !!  KEEP THIS CODE CONSISTENT WITH KERNEL_TEMPLATE.C  !!   #
    #   !!                                                    !!   #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #                                                              #
    ################################################################
    parameters[2:] = values[details.par_offset]
    scale, background = values[0], values[1]
    if details.num_active == 0:
        norm = float(form_volume())
        if norm > 0.0:
            return (scale/norm)*form() + background
        else:
            return np.ones(nq, 'd')*background

    partial_weight = np.NaN
    spherical_correction = 1.0
    pd_stride = details.pd_stride[:details.num_active]
    pd_length = details.pd_length[:details.num_active]
    pd_offset = details.pd_offset[:details.num_active]
    pd_index = np.empty_like(pd_offset)
    offset = np.empty_like(details.par_offset)
    theta = details.theta_par
    fast_length = pd_length[0]
    pd_index[0] = fast_length
    total = np.zeros(nq, 'd')
    norm = 0.0
    for loop_index in range(details.total_pd):
        # update polydispersity parameter values
        if pd_index[0] == fast_length:
            pd_index[:] = (loop_index/pd_stride)%pd_length
            partial_weight = np.prod(weights[pd_offset+pd_index][1:])
            for k in range(details.num_coord):
                par = details.par_coord[k]
                coord = details.pd_coord[k]
                this_offset = details.par_offset[par]
                block_size = 1
                for bit in xrange(32):
                    if coord&1:
                        this_offset += block_size * pd_index[bit]
                        block_size *= pd_length[bit]
                    coord >>= 1
                    if coord == 0: break
                offset[par] = this_offset
                parameters[par] = values[this_offset]
                if par == theta and not (details.par_coord[k]&1):
                    spherical_correction = max(abs(cos(pi/180 * parameters[theta])), 1e-6)
        for k in range(details.num_coord):
            if details.pd_coord[k]&1:
                par = details.par_coord[k]
                parameters[par] = values[offset[par]]
                offset[par] += 1
                if par == theta:
                    spherical_correction = max(abs(cos(pi/180 * parameters[theta])), 1e-6)

        weight = partial_weight * weights[pd_offset[0] + pd_index[0]]
        pd_index[0] += 1
        if weight > cutoff:
            # Call the scattering function
            # Assume that NaNs are only generated if the parameters are bad;
            # exclude all q for that NaN.  Even better would be to have an
            # INVALID expression like the C models, but that is too expensive.
            I = form()
            if np.isnan(I).any(): continue

            # update value and norm
            weight *= spherical_correction
            total += weight * I
            norm += weight * form_volume()

    if norm > 0.0:
        return (scale/norm)*total + background
    else:
        return np.ones(nq, 'd')*background
