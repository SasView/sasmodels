from __future__ import print_function

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
    return call_details

def poly_details(model_info, weights):
    #print("weights",weights)
    #weights = weights[2:] # Skip scale and background

    # Decreasing list of polydispersity lengths
    pd_length = np.array([len(w) for w in weights])
    num_active = np.sum(pd_length>1)
    if num_active > model_info.parameters.max_pd:
        raise ValueError("Too many polydisperse parameters")

    pd_offset = np.cumsum(np.hstack((0, pd_length)))
    #print(", ".join(str(i)+"-"+p.id for i,p in enumerate(model_info.parameters.call_parameters)))
    #print("len:",pd_length)
    #print("off:",pd_offset)
    # Note: the reversing view, x[::-1], does not require a copy
    idx = np.argsort(pd_length)[::-1][:num_active]
    par_length = np.array([len(w) for w in weights])
    pd_stride = np.cumprod(np.hstack((1, par_length[idx])))

    call_details = CallDetails(model_info)
    call_details.pd_par[:num_active] = idx - 2  # skip background & scale
    call_details.pd_length[:num_active] = pd_length[idx]
    call_details.pd_offset[:num_active] = pd_offset[idx]
    call_details.pd_stride[:num_active] = pd_stride[:-1]
    call_details.pd_prod = pd_stride[-1]
    call_details.pd_sum = np.sum(par_length)
    call_details.num_active = num_active
    #call_details.show()
    return call_details
