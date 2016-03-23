"""
Mixture model
-------------

The product model multiplies the structure factor by the form factor,
modulated by the effective radius of the form.  The resulting model
has a attributes of both the model description (with parameters, etc.)
and the module evaluator (with call, release, etc.).

To use it, first load form factor P and structure factor S, then create
*ProductModel(P, S)*.
"""
from copy import copy
import numpy as np

from .modelinfo import Parameter, ParameterTable

SCALE=0
BACKGROUND=1
EFFECT_RADIUS=2
VOLFRACTION=3

def make_mixture_info(parts):
    """
    Create info block for product model.
    """
    flatten = []
    for part in parts:
        if part['composition'] and part['composition'][0] == 'mixture':
            flatten.extend(part['compostion'][1])
        else:
            flatten.append(part)
    parts = flatten

    # Build new parameter list
    pars = []
    for k, part in enumerate(parts):
        # Parameter prefix per model, A_, B_, ...
        # Note that prefix must also be applied to id and length_control
        # to support vector parameters
        prefix = chr(ord('A')+k) + '_'
        pars.append(Parameter(prefix+'scale'))
        for p in part['parameters'].kernel_pars:
            p = copy(p)
            p.name = prefix+p.name
            p.id = prefix+p.id
            if p.length_control is not None:
                p.length_control = prefix+p.length_control
            pars.append(p)
    partable = ParameterTable(pars)

    model_info = {}
    model_info['id'] = '+'.join(part['id'])
    model_info['name'] = ' + '.join(part['name'])
    model_info['filename'] = None
    model_info['title'] = 'Mixture model with ' + model_info['name']
    model_info['description'] = model_info['title']
    model_info['docs'] = model_info['title']
    model_info['category'] = "custom"
    model_info['parameters'] = partable
    #model_info['single'] = any(part['single'] for part in parts)
    model_info['structure_factor'] = False
    model_info['variant_info'] = None
    #model_info['tests'] = []
    #model_info['source'] = []
    # Iq, Iqxy, form_volume, ER, VR and sesans
    # Remember the component info blocks so we can build the model
    model_info['composition'] = ('mixture', parts)


class MixtureModel(object):
    def __init__(self, model_info, parts):
        self.info = model_info
        self.parts = parts

    def __call__(self, q_vectors):
        # Note: may be sending the q_vectors to the n times even though they
        # are only needed once.  It would mess up modularity quite a bit to
        # handle this optimally, especially since there are many cases where
        # separate q vectors are needed (e.g., form in python and structure
        # in opencl; or both in opencl, but one in single precision and the
        # other in double precision).
        kernels = [part(q_vectors) for part in self.parts]
        return MixtureKernel(self.info, kernels)

    def release(self):
        """
        Free resources associated with the model.
        """
        for part in self.parts:
            part.release()


class MixtureKernel(object):
    def __init__(self, model_info, kernels):
        dim = '2d' if kernels[0].q_input.is_2d else '1d'

        # fixed offsets starts at 2 for scale and background
        fixed_pars, pd_pars = [], []
        offsets = [[2, 0]]
        #vol_index = []
        def accumulate(fixed, pd, volume):
            # subtract 1 from fixed since we are removing background
            fixed_offset, pd_offset = offsets[-1]
            #vol_index.extend(k+pd_offset for k,v in pd if v in volume)
            offsets.append([fixed_offset + len(fixed) - 1, pd_offset + len(pd)])
            pd_pars.append(pd)
        if dim == '2d':
            for p in kernels:
                partype = p.info['partype']
                accumulate(partype['fixed-2d'], partype['pd-2d'], partype['volume'])
        else:
            for p in kernels:
                partype = p.info['partype']
                accumulate(partype['fixed-1d'], partype['pd-1d'], partype['volume'])

        #self.vol_index = vol_index
        self.offsets = offsets
        self.fixed_pars = fixed_pars
        self.pd_pars = pd_pars
        self.info = model_info
        self.kernels = kernels
        self.results = None

    def __call__(self, fixed_pars, pd_pars, cutoff=1e-5):
        scale, background = fixed_pars[0:2]
        total = 0.0
        self.results = []  # remember the parts for plotting later
        for k in range(len(self.offsets)-1):
            start_fixed, start_pd = self.offsets[k]
            end_fixed, end_pd = self.offsets[k+1]
            part_fixed = [fixed_pars[start_fixed], 0.0] + fixed_pars[start_fixed+1:end_fixed]
            part_pd = [pd_pars[start_pd], 0.0] + pd_pars[start_pd+1:end_pd]
            part_result = self.kernels[k](part_fixed, part_pd)
            total += part_result
            self.results.append(scale*sum+background)

        return scale*total + background

    def release(self):
        self.p_kernel.release()
        self.q_kernel.release()

