"""
Kernel Call Details
===================

When calling sas computational kernels with polydispersity there are a
number of details that need to be sent to the caller.  This includes the
list of polydisperse parameters, the number of points in the polydispersity
weight distribution, and which parameter is the "theta" parameter for
polar coordinate integration.  The :class:`CallDetails` object maintains
this data.  Use :func:`build_details` to build a *details* object which
can be passed to one of the computational kernels.
"""

from __future__ import print_function

import numpy as np  # type: ignore
from numpy import pi, cos, sin

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

try:
    from typing import List
except ImportError:
    pass
else:
    from .modelinfo import ModelInfo


class CallDetails(object):
    """
    Manage the polydispersity information for the kernel call.

    Conceptually, a polydispersity calculation is an integral over a mesh
    in n-D space where n is the number of polydisperse parameters.  In order
    to keep the program responsive, and not crash the GPU, only a portion
    of the mesh is computed at a time.  Meshes with a large number of points
    will therefore require many calls to the polydispersity loop.  Restarting
    a nested loop in the middle requires that the indices of the individual
    mesh dimensions can be computed for the current loop location.  This
    is handled by the *pd_stride* vector, with n//stride giving the loop
    index and n%stride giving the position in the sub loops.

    One of the parameters may be the latitude.  When integrating in polar
    coordinates, the total circumference decreases as latitude varies from
    pi r^2 at the equator to 0 at the pole, and the weight associated
    with a range of phi values needs to be scaled by this circumference.
    This scale factor needs to be updated each time the theta value
    changes.  *theta_par* indicates which of the values in the parameter
    vector is the latitude parameter, or -1 if there is no latitude
    parameter in the model.  In practice, the normalization term cancels
    if the latitude is not a polydisperse parameter.
    """
    parts = None  # type: List["CallDetails"]
    def __init__(self, model_info):
        # type: (ModelInfo) -> None
        parameters = model_info.parameters
        max_pd = parameters.max_pd

        # Structure of the call details buffer:
        #   pd_par[max_pd]     pd params in order of length
        #   pd_length[max_pd]  length of each pd param
        #   pd_offset[max_pd]  offset of pd values in parameter array
        #   pd_stride[max_pd]  index of pd value in loop = n//stride[k]
        #   pd_prod            total length of pd loop
        #   pd_sum             total length of the weight vector
        #   num_active         number of pd params
        #   theta_par          parameter number for theta parameter
        self.buffer = np.zeros(4*max_pd + 4, 'i4')

        # generate views on different parts of the array
        self._pd_par = self.buffer[0 * max_pd:1 * max_pd]
        self._pd_length = self.buffer[1 * max_pd:2 * max_pd]
        self._pd_offset = self.buffer[2 * max_pd:3 * max_pd]
        self._pd_stride = self.buffer[3 * max_pd:4 * max_pd]

        # theta_par is fixed
        self.theta_par = parameters.theta_offset

    @property
    def pd_par(self):
        """List of polydisperse parameters"""
        return self._pd_par

    @property
    def pd_length(self):
        """Number of weights for each polydisperse parameter"""
        return self._pd_length

    @property
    def pd_offset(self):
        """Offsets for the individual weight vectors in the set of weights"""
        return self._pd_offset

    @property
    def pd_stride(self):
        """Stride in the pd mesh for each pd dimension"""
        return self._pd_stride

    @property
    def pd_prod(self):
        """Total size of the pd mesh"""
        return self.buffer[-4]

    @pd_prod.setter
    def pd_prod(self, v):
        """Total size of the pd mesh"""
        self.buffer[-4] = v

    @property
    def pd_sum(self):
        """Total length of all the weight vectors"""
        return self.buffer[-3]

    @pd_sum.setter
    def pd_sum(self, v):
        """Total length of all the weight vectors"""
        self.buffer[-3] = v

    @property
    def num_active(self):
        """Number of active polydispersity loops"""
        return self.buffer[-2]

    @num_active.setter
    def num_active(self, v):
        """Number of active polydispersity loops"""
        self.buffer[-2] = v

    @property
    def theta_par(self):
        """Location of the theta parameter in the parameter vector"""
        return self.buffer[-1]

    @theta_par.setter
    def theta_par(self, v):
        """Location of the theta parameter in the parameter vector"""
        self.buffer[-1] = v

    def show(self):
        """Print the polydispersity call details to the console"""
        print("num_active", self.num_active)
        print("pd_prod", self.pd_prod)
        print("pd_sum", self.pd_sum)
        print("theta par", self.theta_par)
        print("pd_par", self.pd_par)
        print("pd_length", self.pd_length)
        print("pd_offset", self.pd_offset)
        print("pd_stride", self.pd_stride)


def mono_details(model_info):
    # type: (ModelInfo) -> CallDetails
    """
    Return a :class:`CallDetails` object for a monodisperse calculation
    of the model defined by *model_info*.
    """
    call_details = CallDetails(model_info)
    call_details.pd_prod = 1
    call_details.pd_sum = model_info.parameters.nvalues
    call_details.pd_par[:] = np.arange(0, model_info.parameters.max_pd)
    call_details.pd_length[:] = 1
    call_details.pd_offset[:] = np.arange(0, model_info.parameters.max_pd)
    call_details.pd_stride[:] = 1
    return call_details


def poly_details(model_info, weights):
    """
    Return a :class:`CallDetails` object for a polydisperse calculation
    of the model defined by *model_info* for the given set of *weights*.
    *weights* is a list of vectors, one for each parameter.  Monodisperse
    parameters should provide a weight vector containing [1.0].
    """
    # type: (ModelInfo) -> CallDetails
    #print("weights",weights)
    #weights = weights[2:] # Skip scale and background

    # Decreasing list of polydispersity lengths
    #print([p.id for p in model_info.parameters.call_parameters])
    pd_length = np.array([len(w)
                          for w in weights[2:2+model_info.parameters.npars]])
    num_active = np.sum(pd_length > 1)
    max_pd = model_info.parameters.max_pd
    if num_active > max_pd:
        raise ValueError("Too many polydisperse parameters")

    pd_offset = np.cumsum(np.hstack((0, pd_length)))
    #print(", ".join(str(i)+"-"+p.id for i,p in enumerate(model_info.parameters.call_parameters)))
    #print("len:",pd_length)
    #print("off:",pd_offset)
    # Note: the reversing view, x[::-1], does not require a copy
    idx = np.argsort(pd_length)[::-1][:max_pd]
    pd_stride = np.cumprod(np.hstack((1, pd_length[idx])))

    call_details = CallDetails(model_info)
    call_details.pd_par[:max_pd] = idx
    call_details.pd_length[:max_pd] = pd_length[idx]
    call_details.pd_offset[:max_pd] = pd_offset[idx]
    call_details.pd_stride[:max_pd] = pd_stride[:-1]
    call_details.pd_prod = pd_stride[-1]
    call_details.pd_sum = sum(len(w) for w in weights)
    call_details.num_active = num_active
    #call_details.show()
    return call_details


def dispersion_mesh(model_info, pars):
    # type: (ModelInfo) -> Tuple[List[np.ndarray], List[np.ndarray]]
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
            pars.append(np.vstack(value[offset:offset+n])
                        if n > 1 else value[offset])
            offset += n
        value = pars
    return value, weight



def build_details(kernel, pairs):
    # type: (Kernel, Tuple[List[np.ndarray], List[np.ndarray]]) -> Tuple[CallDetails, np.ndarray, bool]
    """
    Converts (value, weight) pairs into parameters for the kernel call.

    Returns a CallDetails object indicating the polydispersity, a data object
    containing the different values, and the magnetic flag indicating whether
    any magnetic magnitudes are non-zero. Magnetic vectors (M0, phi, theta) are
    converted to rectangular coordinates (mx, my, mz).
    """
    values, weights = zip(*pairs)
    scalars = [v[0] for v in values]
    if all(len(w) == 1 for w in weights):
        call_details = mono_details(kernel.info)
        # Pad value array to a 32 value boundary
        data_len = 3*len(scalars)
        extra = ((data_len+31)//32)*32 - data_len
        data = np.array(scalars+scalars+[1.]*len(scalars)+[0.]*extra,
                        dtype=kernel.dtype)
    else:
        call_details = poly_details(kernel.info, weights)
        # Pad value array to a 32 value boundary
        data_len = len(scalars) + 2*sum(len(v) for v in values)
        extra = ((data_len+31)//32)*32 - data_len
        data = np.hstack(scalars+list(values)+list(weights)+[0.]*extra)
        data = data.astype(kernel.dtype)
    is_magnetic = convert_magnetism(kernel.info.parameters, data)
    #call_details.show()
    return call_details, data, is_magnetic

def convert_magnetism(parameters, values):
    """
    Convert magnetism values from polar to rectangular coordinates.

    Returns True if any magnetism is present.
    """
    mag = values[parameters.nvalues-3*parameters.nmagnetic:parameters.nvalues]
    mag = mag.reshape(-1, 3)
    scale = mag[:,0]
    if np.any(scale):
        theta, phi = mag[:, 1]*pi/180., mag[:, 2]*pi/180.
        cos_theta = cos(theta)
        mag[:, 0] = scale*cos_theta*cos(phi)  # mx
        mag[:, 1] = scale*sin(theta)  # my
        mag[:, 2] = -scale*cos_theta*sin(phi)  # mz
        return True
    else:
        return False
