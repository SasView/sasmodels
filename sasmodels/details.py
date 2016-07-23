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
        self._pd_par     = self.buffer[0 * max_pd:1 * max_pd]
        self._pd_length  = self.buffer[1 * max_pd:2 * max_pd]
        self._pd_offset  = self.buffer[2 * max_pd:3 * max_pd]
        self._pd_stride  = self.buffer[3 * max_pd:4 * max_pd]

        # theta_par is fixed
        self.theta_par = parameters.theta_offset

    @property
    def pd_par(self): return self._pd_par

    @property
    def pd_length(self): return self._pd_length

    @property
    def pd_offset(self): return self._pd_offset

    @property
    def pd_stride(self): return self._pd_stride

    @property
    def pd_prod(self): return self.buffer[-4]
    @pd_prod.setter
    def pd_prod(self, v): self.buffer[-4] = v

    @property
    def pd_sum(self): return self.buffer[-3]
    @pd_sum.setter
    def pd_sum(self, v): self.buffer[-3] = v

    @property
    def num_active(self): return self.buffer[-2]
    @num_active.setter
    def num_active(self, v): self.buffer[-2] = v

    @property
    def theta_par(self): return self.buffer[-1]
    @theta_par.setter
    def theta_par(self, v): self.buffer[-1] = v

    def show(self):
        print("num_active", self.num_active)
        print("pd_prod", self.pd_prod)
        print("pd_sum", self.pd_sum)
        print("theta par", self.theta_par)
        print("pd_par", self.pd_par)
        print("pd_length", self.pd_length)
        print("pd_offset", self.pd_offset)
        print("pd_stride", self.pd_stride)


def mono_details(model_info):
    call_details = CallDetails(model_info)
    call_details.pd_prod = 1
    call_details.pd_sum = model_info.parameters.nvalues
    call_details.pd_par[:] = np.arange(0, model_info.parameters.max_pd)
    call_details.pd_length[:] = 1
    call_details.pd_offset[:] = np.arange(0, model_info.parameters.max_pd)
    call_details.pd_stride[:] = 1
    return call_details


def poly_details(model_info, weights):
    #print("weights",weights)
    #weights = weights[2:] # Skip scale and background

    # Decreasing list of polydispersity lengths
    #print([p.id for p in model_info.parameters.call_parameters])
    pd_length = np.array([len(w) for w in weights[2:2+model_info.parameters.npars]])
    num_active = np.sum(pd_length>1)
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
    if all(len(w)==1 for w in weights):
        call_details = mono_details(kernel.info)
        data = np.array(scalars+scalars+[1]*len(scalars), dtype=kernel.dtype)
    else:
        call_details = poly_details(kernel.info, weights)
        data = np.hstack(scalars+list(values)+list(weights)).astype(kernel.dtype)
    is_magnetic = convert_magnetism(kernel.info.parameters, data)
    #call_details.show()
    return call_details, data, is_magnetic

def convert_magnetism(parameters, values):
    """
    Convert magnetism in value table from polar to rectangular coordinates.

    Returns True if any magnetism is present.
    """
    mag = values[parameters.nvalues-3*parameters.nmagnetic:parameters.nvalues]
    mag = mag.reshape(-1, 3)
    M0 = mag[:,0]
    if np.any(M0):
        theta, phi = mag[:,1]*pi/180., mag[:,2]*pi/180.
        cos_theta = cos(theta)
        mx = M0*cos_theta*cos(phi)
        my = M0*sin(theta)
        mz = -M0*cos_theta*sin(phi)
        mag[:,0], mag[:,1], mag[:,2] = mx, my, mz
        return True
    else:
        return False
