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
from numpy import cos, sin, radians

try:
    np.meshgrid([])
    meshgrid = np.meshgrid
except Exception:
    # CRUFT: np.meshgrid requires multiple vectors
    def meshgrid(*args):
        """See docs from a recent version of numpy"""
        if len(args) > 1:
            return np.meshgrid(*args)
        else:
            return [np.asarray(v) for v in args]

# pylint: disable=unused-import
try:
    from typing import List, Tuple, Sequence
except ImportError:
    pass
else:
    from .modelinfo import ModelInfo
    from .modelinfo import ParameterTable
# pylint: enable=unused-import


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
    with a range of latitude values needs to be scaled by this circumference.
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
        #   num_eval           total length of pd loop
        #   num_weights        total length of the weight vector
        #   num_active         number of pd params
        #   theta_par          parameter number for theta parameter
        self.buffer = np.empty(4*max_pd + 4, 'i4')

        # generate views on different parts of the array
        self._pd_par = self.buffer[0 * max_pd:1 * max_pd]
        self._pd_length = self.buffer[1 * max_pd:2 * max_pd]
        self._pd_offset = self.buffer[2 * max_pd:3 * max_pd]
        self._pd_stride = self.buffer[3 * max_pd:4 * max_pd]

        # theta_par is fixed
        self.theta_par = parameters.theta_offset

        # offset and length are for all parameters, not just pd parameters
        # They are not sent to the kernel function, though they could be.
        # They are used by the composite models (sum and product) to
        # figure out offsets into the combined value list.
        self.offset = None  # type: np.ndarray
        self.length = None  # type: np.ndarray

        # keep hold of ifno show() so we can break a values vector
        # into the individual components
        self.info = model_info

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
    def num_eval(self):
        """Total size of the pd mesh"""
        return self.buffer[-4]

    @num_eval.setter
    def num_eval(self, v):
        """Total size of the pd mesh"""
        self.buffer[-4] = v

    @property
    def num_weights(self):
        """Total length of all the weight vectors"""
        return self.buffer[-3]

    @num_weights.setter
    def num_weights(self, v):
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

    def show(self, values=None):
        """Print the polydispersity call details to the console"""
        print("===== %s details ===="%self.info.name)
        print("num_active:%d  num_eval:%d  num_weights:%d  theta=%d"
              % (self.num_active, self.num_eval, self.num_weights, self.theta_par))
        if self.pd_par.size:
            print("pd_par", self.pd_par)
            print("pd_length", self.pd_length)
            print("pd_offset", self.pd_offset)
            print("pd_stride", self.pd_stride)
        if values is not None:
            nvalues = self.info.parameters.nvalues
            print("scale, background", values[:2])
            print("val", values[2:nvalues])
            print("pd", values[nvalues:nvalues+self.num_weights])
            print("wt", values[nvalues+self.num_weights:nvalues+2*self.num_weights])
            print("offsets", self.offset)


def make_details(model_info, length, offset, num_weights):
    # type: (ModelInfo, np.ndarray, np.ndarray, int) -> CallDetails
    """
    Return a :class:`CallDetails` object for a polydisperse calculation
    of the model defined by *model_info*.  Polydispersity is defined by
    the *length* of the polydispersity distribution for each parameter
    and the *offset* of the distribution in the polydispersity array.
    Monodisperse parameters should use a polydispersity length of one
    with weight 1.0. *num_weights* is the total length of the polydispersity
    array.
    """
    #pars = model_info.parameters.call_parameters[2:model_info.parameters.npars+2]
    #print(", ".join(str(i)+"-"+p.id for i,p in enumerate(pars)))
    #print("len:",length)
    #print("off:",offset)

    # Check that we aren't using too many polydispersity loops
    num_active = np.sum(length > 1)
    max_pd = model_info.parameters.max_pd
    if num_active > max_pd:
        raise ValueError("Too many polydisperse parameters")

    # Decreasing list of polydpersity lengths
    # Note: the reversing view, x[::-1], does not require a copy
    idx = np.argsort(length)[::-1][:max_pd]
    pd_stride = np.cumprod(np.hstack((1, length[idx])))

    call_details = CallDetails(model_info)
    call_details.pd_par[:max_pd] = idx
    call_details.pd_length[:max_pd] = length[idx]
    call_details.pd_offset[:max_pd] = offset[idx]
    call_details.pd_stride[:max_pd] = pd_stride[:-1]
    call_details.num_eval = pd_stride[-1]
    call_details.num_weights = num_weights
    call_details.num_active = num_active
    call_details.length = length
    call_details.offset = offset
    #call_details.show()
    return call_details


ZEROS = tuple([0.]*31)
def make_kernel_args(kernel, # type: Kernel
                     mesh    # type: Tuple[List[np.ndarray], List[np.ndarray]]
                    ):
    # type: (...) -> Tuple[CallDetails, np.ndarray, bool]
    """
    Converts (value, dispersity, weight) for each parameter into kernel pars.

    Returns a CallDetails object indicating the polydispersity, a data object
    containing the different values, and the magnetic flag indicating whether
    any magnetic magnitudes are non-zero. Magnetic vectors (M0, phi, theta) are
    converted to rectangular coordinates (mx, my, mz).
    """
    npars = kernel.info.parameters.npars
    nvalues = kernel.info.parameters.nvalues
    scalars = [value for value, _dispersity, _weight in mesh]
    # skipping scale and background when building values and weights
    _values, dispersity, weights = zip(*mesh[2:npars+2]) if npars else ((), (), ())
    #weights = correct_theta_weights(kernel.info.parameters, dispersity, weights)
    length = np.array([len(w) for w in weights])
    offset = np.cumsum(np.hstack((0, length)))
    call_details = make_details(kernel.info, length, offset[:-1], offset[-1])
    # Pad value array to a 32 value boundary
    data_len = nvalues + 2*sum(len(v) for v in dispersity)
    extra = (32 - data_len%32)%32
    data = np.hstack((scalars,) + dispersity + weights + ZEROS[:extra])
    data = data.astype(kernel.dtype)
    is_magnetic = convert_magnetism(kernel.info.parameters, data)
    #call_details.show()
    #print("data", data)
    return call_details, data, is_magnetic

def correct_theta_weights(parameters, # type: ParameterTable
                          dispersity, # type: Sequence[np.ndarray]
                          weights     # type: Sequence[np.ndarray]
                         ):
    # type: (...) -> Sequence[np.ndarray]
    """
    **Deprecated** Theta weights will be computed in the kernel wrapper if
    they are needed.

    If there is a theta parameter, update the weights of that parameter so that
    the cosine weighting required for polar integration is preserved.

    Avoid evaluation strictly at the pole, which would otherwise send the
    weight to zero.  This is probably not a problem in practice (if dispersity
    is +/- 90, then you probably should be using a 1-D model of the circular
    average).

    Note: scale and background parameters are not include in the tuples for
    dispersity and weights, so index is parameters.theta_offset, not
    parameters.theta_offset+2

    Returns updated weights vectors
    """
    # Apparently the parameters.theta_offset similarly skips scale and
    # and background, so the indexing works out, but they are still shipped
    # to the kernel, so we need to add two there.
    if parameters.theta_offset >= 0:
        index = parameters.theta_offset
        theta = dispersity[index]
        theta_weight = abs(cos(radians(theta)))
        weights = tuple(theta_weight*w if k == index else w
                        for k, w in enumerate(weights))
    return weights


def convert_magnetism(parameters, values):
    # type: (ParameterTable, Sequence[np.ndarray]) -> bool
    """
    Convert magnetism values from polar to rectangular coordinates.

    Returns True if any magnetism is present.
    """
    mag = values[parameters.nvalues-3*parameters.nmagnetic:parameters.nvalues]
    mag = mag.reshape(-1, 3)
    if np.any(mag[:, 0] != 0.0):
        M0 = mag[:, 0].copy()
        theta, phi = radians(mag[:, 1]), radians(mag[:, 2])
        mag[:, 0] = +M0*cos(theta)*cos(phi)  # mx
        mag[:, 1] = +M0*sin(theta) # my
        mag[:, 2] = -M0*cos(theta)*sin(phi)  # mz
        return True
    else:
        return False


def dispersion_mesh(model_info, mesh):
    # type: (ModelInfo) -> Tuple[List[np.ndarray], List[np.ndarray]]
    """
    Create a mesh grid of dispersion parameters and weights.

    *mesh* is a list of (value, dispersity, weights), where the values
    are the individual parameter values, and (dispersity, weights) is
    the distribution of parameter values.

    Only the volume parameters should be included in this list.  Orientation
    parameters do not affect the calculation of effective radius or volume
    ratio.  This is convenient since it avoids the distinction between
    value and dispersity that is present in orientation parameters but not
    shape parameters.

    Returns [p1,p2,...],w where pj is a vector of values for parameter j
    and w is a vector containing the products for weights for each
    parameter set in the vector.
    """
    _, dispersity, weight = zip(*mesh)
    #weight = [w if len(w)>0 else [1.] for w in weight]
    weight = np.vstack([v.flatten() for v in meshgrid(*weight)])
    weight = np.prod(weight, axis=0)
    dispersity = [v.flatten() for v in meshgrid(*dispersity)]
    lengths = [par.length for par in model_info.parameters.kernel_parameters
               if par.type == 'volume']
    if any(n > 1 for n in lengths):
        pars = []
        offset = 0
        for n in lengths:
            pars.append(np.vstack(dispersity[offset:offset+n])
                        if n > 1 else dispersity[offset])
            offset += n
        dispersity = pars
    return dispersity, weight
