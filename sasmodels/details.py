import numpy as np  # type: ignore

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
        npars = parameters.npars
        par_offset = 4*max_pd
        self.buffer = np.zeros(par_offset + 3 * npars + 4, 'i4')

        # generate views on different parts of the array
        self._pd_par     = self.buffer[0 * max_pd:1 * max_pd]
        self._pd_length  = self.buffer[1 * max_pd:2 * max_pd]
        self._pd_offset  = self.buffer[2 * max_pd:3 * max_pd]
        self._pd_stride  = self.buffer[3 * max_pd:4 * max_pd]
        self._par_offset = self.buffer[par_offset + 0 * npars:par_offset + 1 * npars]
        self._par_coord  = self.buffer[par_offset + 1 * npars:par_offset + 2 * npars]
        self._pd_coord   = self.buffer[par_offset + 2 * npars:par_offset + 3 * npars]

        # theta_par is fixed
        self.buffer[-1] = parameters.theta_offset

    @property
    def pd_par(self): return self._pd_par

    @property
    def pd_length(self): return self._pd_length

    @property
    def pd_offset(self): return self._pd_offset

    @property
    def pd_stride(self): return self._pd_stride

    @property
    def pd_coord(self): return self._pd_coord

    @property
    def par_coord(self): return self._par_coord

    @property
    def par_offset(self): return self._par_offset

    @property
    def num_active(self): return self.buffer[-4]
    @num_active.setter
    def num_active(self, v): self.buffer[-4] = v

    @property
    def total_pd(self): return self.buffer[-3]
    @total_pd.setter
    def total_pd(self, v): self.buffer[-3] = v

    @property
    def num_coord(self): return self.buffer[-2]
    @num_coord.setter
    def num_coord(self, v): self.buffer[-2] = v

    @property
    def theta_par(self): return self.buffer[-1]

    def show(self):
        print("total_pd", self.total_pd)
        print("num_active", self.num_active)
        print("pd_par", self.pd_par)
        print("pd_length", self.pd_length)
        print("pd_offset", self.pd_offset)
        print("pd_stride", self.pd_stride)
        print("par_offsets", self.par_offset)
        print("num_coord", self.num_coord)
        print("par_coord", self.par_coord)
        print("pd_coord", self.pd_coord)
        print("theta par", self.buffer[-1])

def mono_details(model_info):
    call_details = CallDetails(model_info)
    # The zero defaults for monodisperse systems are mostly fine
    call_details.par_offset[:] = np.arange(2, len(call_details.par_offset)+2)
    return call_details

def poly_details(model_info, weights):
    #print("weights",weights)
    weights = weights[2:] # Skip scale and background

    # Decreasing list of polydispersity lengths
    # Note: the reversing view, x[::-1], does not require a copy
    pd_length = np.array([len(w) for w in weights])
    num_active = np.sum(pd_length>1)
    if num_active > model_info.parameters.max_pd:
        raise ValueError("Too many polydisperse parameters")

    pd_offset = np.cumsum(np.hstack((0, pd_length)))
    idx = np.argsort(pd_length)[::-1][:num_active]
    par_length = np.array([max(len(w),1) for w in weights])
    pd_stride = np.cumprod(np.hstack((1, par_length[idx])))
    par_offsets = np.cumsum(np.hstack((2, par_length)))

    call_details = CallDetails(model_info)
    call_details.pd_par[:num_active] = idx
    call_details.pd_length[:num_active] = pd_length[idx]
    call_details.pd_offset[:num_active] = pd_offset[idx]
    call_details.pd_stride[:num_active] = pd_stride[:-1]
    call_details.par_offset[:] = par_offsets[:-1]
    call_details.total_pd = pd_stride[-1]
    call_details.num_active = num_active
    # Without constraints coordinated parameters are just the pd parameters
    call_details.par_coord[:num_active] = idx
    call_details.pd_coord[:num_active] = 2**np.arange(num_active)
    call_details.num_coord = num_active
    #call_details.show()
    return call_details

def constrained_poly_details(model_info, weights, constraints):
    # Need to find the independently varying pars and sort them
    # Need to build a coordination list for the dependent variables
    # Need to generate a constraints function which takes values
    # and weights, returning par blocks
    raise NotImplementedError("Can't handle constraints yet")

