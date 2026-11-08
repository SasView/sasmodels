"""
Python driver for python kernels

Calls the kernel with a vector of $q$ values for a single parameter set.
Polydispersity is supported by looping over different parameter sets and
summing the results.  The interface to :class:`PyModel` matches those for
:class:`.kernelcl.GpuModel` and :class:`.kerneldll.DllModel`.
"""

import logging

import numpy as np  # type: ignore
from numpy import pi

try:
    from numpy import cbrt
except ImportError:
    def cbrt(x):
        """Return cubed root of x."""
        return x ** (1.0/3.0)

from .generate import F64
from .jitter import Rx, Ry, Rz
from .kernel import Kernel, KernelModel

# pylint: disable=unused-import
try:
    from typing import Callable

    from .details import CallDetails
    from .modelinfo import ModelInfo
except ImportError:
    pass
# pylint: enable=unused-import

logger = logging.getLogger(__name__)


class PyModel(KernelModel):
    """
    Wrapper for pure python models.
    """
    def __init__(self, model_info):
        # Make sure Iq is available and vectorized.
        _create_default_functions(model_info)
        self.info = model_info
        self.dtype = np.dtype('d')
        logger.info("make python model %s", self.info.name)

    def make_kernel(self, q_vectors, wavelength=None):
        """Instantiate the python kernel with input *q_vectors*"""
        q_input = PyInput(q_vectors, dtype=F64, wavelength=wavelength)
        return PyKernel(self.info, q_input)

    def release(self):
        """
        Free resources associated with the model.
        """
        pass


class PyInput:
    """
    Make q data available to the gpu.

    *q_vectors* is a list of q vectors, which will be *[q]* for 1-D data,
    and *[qx, qy]* for 2-D data.  Internally, the vectors will be reallocated
    to get the best performance on OpenCL, which may involve shifting and
    stretching the array to better match the memory architecture.  Additional
    points will be evaluated with *q=1e-3*.

    *dtype* is the data type for the q vectors. The data type should be
    set to match that of the kernel, which is an attribute of
    :class:`PyModel`.  Note that not all kernels support double
    precision, so even if the program was created for double precision,
    the *GpuProgram.dtype* may be single precision.

    Call :meth:`release` when complete.  Even if not called directly, the
    buffer will be released when the data object is freed.
    """
    def __init__(self, q_vectors, dtype, wavelength=None):
        self.nq = q_vectors[0].size
        self.dtype = dtype
        self.is_2d = 2 <= len(q_vectors) <= 3
        if self.is_2d:
            qx, qy = q_vectors[:2]
            self.q = np.empty((self.nq, 2), dtype=dtype)
            self.q[:, 0] = qx
            self.q[:, 1] = qy
            # TODO: should be putting this in self.q
            if len(q_vectors) == 3:
                self.qz = q_vectors[3]
            else:
                self.qz = deduce_qz(qx, qy, wavelength)
        else:
            self.q = np.empty(self.nq, dtype=dtype)
            self.q[:self.nq] = q_vectors[0]

    def release(self):
        """
        Free resources associated with the model inputs.
        """
        self.q = None


def deduce_qz(qx, qy, wavelength):
    if wavelength is not None:
        # we can deduce the value of qz from qx and qy
        # if we have the wavelength
        # this is how you convert qx, qy, and wavelength to qz
        k0 = 2.0 * np.pi / wavelength
        qz = k0 - np.sqrt(k0**2 - qx**2 - qy**2)
    else:
        qz = np.zeros_like(qx)

    return qz


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
        # type: (ModelInfo, list[np.ndarray]) -> None
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

        # Create views into the array to hold the arguments.
        offset = 0
        kernel_args, volume_args = [], []
        for p in partable.kernel_parameters:
            if p.type == "orientation":
                continue
            if p.length == 1:
                # Scalar values are length 1 vectors with no dimensions.
                # The following expression captures the location of the
                # parameter in the vector so that it can be updated according
                # to the index into the polydispersity mesh within the _loops
                # method below.
                v = parameter_vector[offset:offset+1].reshape(())
            else:
                # Vector values are simple views.
                v = parameter_vector[offset:offset+p.length]
            offset += p.length
            if p in kernel_parameters:
                kernel_args.append(v)
            if p in volume_parameters:
                volume_args.append(v)

        # Hold on to the parameter vector so it doesn't get collected
        # This may also be required to preserve the views into the vector.
        self._parameter_vector = parameter_vector
        # Similarly for the orientation matrix

        # Generate a closure which calls the kernel with the views into the
        # parameter array.
        if q_input.is_2d:
            qx, qy, qz = q_input.q[:, 0], q_input.q[:, 1], q_input.qz
            assert (qz == 0.).all()
            if model_info.Iqabc:
                def _form():
                    dtheta, dphi, dpsi = parameter_vector[-3:]
                    #jitter = Rx(dphi)@Ry(dtheta)@Rz(dpsi)
                    jitter = Rz(-dpsi)@Ry(-dtheta)@Rx(-dphi)
                    q = self._rotate@jitter@np.vstack((qx, qy, qz))
                    # prism with n_sides=4 is rotated 45 degrees relative to parallelepiped
                    q = Rz(-45)@q
                    return model_info.Iqabc(q[0], q[1], q[2], *kernel_args)
                self._form = _form
            elif model_info.Iqac:
                def _form():
                    dtheta, dphi = parameter_vector[-2:]
                    jitter = Ry(-dtheta)@Rz(-dphi)
                    q = self._rotate@jitter@np.vstack((qx, qy, qz))
                    return model_info.Iqac(q[0], q[2], *kernel_args)
                self._form = _form
            else:
                form = model_info.Iqxy
                self._form = lambda: form(qx, qy, *kernel_args)
            self._nq_copies = 1
        else:
            end = -3 if model_info.Iqabc else -2 if model_info.Iqac else len(parameter_vector)
            q = q_input.q
            if model_info.have_Fq:
                def _form():
                    return _pack_fq(*model_info.Fq(q, *kernel_args))
                self._nq_copies = 2
            else:
                form = model_info.Iq
                def _form():
                    return model_info.Iq(q, *kernel_args)
                self._form = _form
                self._nq_copies = 1
            self._form = _form

        # Generate a closure which calls the form_volume if it exists.
        self._volume_args = volume_args
        volume = model_info.form_volume
        shell = model_info.shell_volume
        radius = model_info.radius_effective
        self._volume = (
            (lambda: (shell(*volume_args), volume(*volume_args))) if shell and volume
            else (lambda: [volume(*volume_args)]*2) if volume
            else (lambda: (1.0, 1.0)))
        self._radius = (
            (lambda mode: radius(mode, *volume_args)) if radius
            else (lambda mode: cbrt(0.75/pi*volume(*volume_args))) if volume
            else (lambda mode: 1.0))

    def _call_kernel(self, call_details, values, cutoff, magnetic, radius_effective_mode):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool) -> None
        if magnetic:
            raise NotImplementedError("Magnetism not implemented for pure python models")
        #print("Calling python kernel")
        #call_details.show(values)
        form_radius = (
            (lambda: 0.0)
            if radius_effective_mode == 0
            else (lambda: self._radius(radius_effective_mode))
        )
        self.result = self._loops(form_radius, call_details, values, cutoff)

    def _loops(
        self,
        form_radius: Callable[[], float],
        call_details: CallDetails,
        values: np.ndarray,
        cutoff: float
    ) -> None:
        ################################################################
        #                                                              #
        #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
        #   !!                                                    !!   #
        #   !!  KEEP THIS CODE CONSISTENT WITH KERNEL_TEMPLATE.C  !!   #
        #   !!                                                    !!   #
        #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
        #                                                              #
        ################################################################

        # WARNING: Trickery ahead
        # The parameters[] vector is embedded in the closures for form(),
        # form_volume() and form_radius().  We set the initial vector from
        # the values for the model parameters. As we loop through the polydispesity
        # mesh, we update the components with the polydispersity values before
        # calling the respective functions.
        parameters = self._parameter_vector
        n_pars = len(parameters)
        # Skip scale and background
        parameters[:] = values[2:n_pars+2]
        nq = self.q_input.nq*self._nq_copies

        if self.info.Iqabc:
            theta, phi, psi = parameters[-3:]
            self._rotate = Rz(-psi)@Ry(-theta)@Rz(-phi)
            parameters[-3:] = 0. # Default is no jitter

            _ = """
            sin_theta, cos_theta = np.sin(np.degrees(theta)), np.cos(np.degrees(theta))
            sin_phi, cos_phi = np.sin(np.degrees(phi)), np.cos(np.degrees(phi))
            sin_psi, cos_psi = np.sin(np.degrees(psi)), np.cos(np.degrees(psi))
            V11 = -sin_phi*sin_psi + cos_phi*cos_psi*cos_theta
            V12 = sin_phi*cos_psi*cos_theta + sin_psi*cos_phi
            V21 = -sin_phi*cos_psi - sin_psi*cos_phi*cos_theta
            V22 = -sin_phi*sin_psi*cos_theta + cos_phi*cos_psi
            V31 = sin_theta*cos_phi
            V32 = sin_phi*sin_theta
            rot = np.array([[V11, V12], [V21, V22], [V31, V32]])
            if (rot != self._rotate[:2]).any():
                print(rot)
                print(self._rotate[:2])
                raise Hell
            """

        elif self.info.Iqac:
            theta, phi = parameters[-2:]
            self._rotate = Ry(-theta)@Rz(-phi)
            parameters[-2:] = 0. # Default is no jitter
        else:
            # Unused. self._rotate is used within form, which
            # is specialized for the Iqabc and Iqac cases.
            self._rotate = np.eye(3)

        if call_details.num_active == 0:
            total = self._form()
            weight_norm = 1.0
            weighted_shell, weighted_form = self._volume()
            weighted_radius = form_radius()

        else:
            pd_value = values[2+n_pars:2+n_pars + call_details.num_weights]
            pd_weight = values[2+n_pars + call_details.num_weights:]

            weight_norm = 0.0
            weighted_form = 0.0
            weighted_shell = 0.0
            weighted_radius = 0.0
            partial_weight = np.nan
            weight = np.nan

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
                # Update polydispersity parameter values.
                if p0_index == p0_length:
                    pd_index = (loop_index//pd_stride)%pd_length
                    parameters[pd_par] = pd_value[pd_offset+pd_index]
                    partial_weight = np.prod(pd_weight[pd_offset+pd_index][1:])
                    p0_index = loop_index%p0_length

                weight = partial_weight * pd_weight[p0_offset + p0_index]
                parameters[p0_par] = pd_value[p0_offset + p0_index]
                p0_index += 1
                if weight > cutoff:
                    # Call the scattering function.
                    # Assume that NaNs are only generated if the parameters are bad;
                    # exclude all q for that NaN.  Even better would be to have an
                    # INVALID expression like the C models, but that is expensive.
                    Iq = np.asarray(self._form(), 'd')
                    if np.isnan(Iq).any():
                        continue

                    # Update value and norm.
                    total += weight * Iq
                    weight_norm += weight
                    unweighted_shell, unweighted_form = self._volume()
                    weighted_shell += weight * unweighted_shell
                    weighted_form += weight * unweighted_form
                    weighted_radius += weight * form_radius()

        result = np.hstack((total, weight_norm, weighted_form, weighted_shell, weighted_radius))
        return result

    def release(self):
        # type: () -> None
        """
        Free resources associated with the kernel.
        """
        self.q_input.release()
        self.q_input = None


def _pack_fq(fq, fqsq):
    result = np.empty(2*fq.shape[0], dtype=fq.dtype)
    # PAK: For reasons unknown, the original author (me) decided to pack the
    # values into the returned array as (f^2, f) pairs.
    result[0::2] = fqsq
    result[1::2] = fq
    return result


def _create_default_functions(model_info):
    """
    Autogenerate missing functions, such as Iqxy from Iq.

    This only works for Iqxy when Iq is written in python. :func:`make_source`
    performs a similar role for Iq written in C.  This also vectorizes
    any functions that are not already marked as vectorized.
    """
    # Note: Must call create_vector_Iq before create_vector_Iqxy.
    _create_vector_Iq(model_info)
    _create_vector_Fq(model_info)
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

def _create_vector_Fq(model_info):
    """
    Define Fq as a vector function if it exists.
    """
    # Note that this is doing slightly too much work since we are composing
    # two vector results which are then interlaced via _pack_fq above.
    Fq = model_info.Fq
    if callable(Fq) and not getattr(Fq, 'vectorized', False):
        def vector_Fq(q, *args):
            """
            Vectorized 1D kernel returning Fq, Fq**2
            """
            fq, fqsq = zip(*(Fq(qi, *args) for qi in q))
            return np.array(fq), np.array(fqsq)
        vector_Fq.vectorized = True
        model_info.Fq = vector_Fq

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
