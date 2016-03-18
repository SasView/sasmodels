"""
Product model
-------------

The product model multiplies the structure factor by the form factor,
modulated by the effective radius of the form.  The resulting model
has a attributes of both the model description (with parameters, etc.)
and the module evaluator (with call, release, etc.).

To use it, first load form factor P and structure factor S, then create
*ProductModel(P, S)*.
"""
import numpy as np

from .core import call_ER_VR
from .generate import process_parameters

SCALE=0
BACKGROUND=1
RADIUS_EFFECTIVE=2
VOLFRACTION=3

def make_product_info(p_info, s_info):
    """
    Create info block for product model.
    """
    p_id, p_name, p_pars = p_info['id'], p_info['name'], p_info['parameters']
    s_id, s_name, s_pars = s_info['id'], s_info['name'], s_info['parameters']
    # We require models to start with scale and background
    assert s_pars[SCALE].name == 'scale'
    assert s_pars[BACKGROUND].name == 'background'
    # We require structure factors to start with effect radius and volfraction
    assert s_pars[RADIUS_EFFECTIVE].name == 'radius_effective'
    assert s_pars[VOLFRACTION].name == 'volfraction'
    # Combine the parameter sets.  We are skipping the first three
    # parameters of S since scale, background are defined in P and
    # effect_radius is set from P.ER().
    pars = p_pars + s_pars[3:]
    # check for duplicates; can't use assertion since they may never be checked
    if len(set(p.name for p in pars)) != len(pars):
        raise ValueError("Duplicate parameters in %s and %s"%(p_id))
    # For comparison with sasview, determine the old parameters.
    oldname = [p_info['oldname'], s_info['oldname']]
    oldpars = {'scale':'scale_factor'}
    oldpars.update(p_info['oldpars'])
    oldpars.update(s_info['oldpars'])

    model_info = {}
    model_info['id'] = '*'.join((p_id, s_id))
    model_info['name'] = ' X '.join((p_name, s_name))
    model_info['filename'] = None
    model_info['title'] = 'Product of %s and structure factor %s'%(p_name, s_name)
    model_info['description'] = model_info['title']
    model_info['docs'] = model_info['title']
    model_info['category'] = "custom"
    model_info['parameters'] = pars
    #model_info['single'] = p_info['single'] and s_info['single']
    model_info['structure_factor'] = False
    model_info['variant_info'] = None
    #model_info['tests'] = []
    #model_info['source'] = []
    # Iq, Iqxy, form_volume, ER, VR and sesans
    model_info['oldname'] = oldname
    model_info['oldpars'] = oldpars
    model_info['composition'] = ('product', [p_info, s_info])
    process_parameters(model_info)
    return model_info

class ProductModel(object):
    def __init__(self, P, S):
        self.P = P
        self.S = S
        self.info = make_product_info(P.info, S.info)

    def __call__(self, q_vectors):
        # Note: may be sending the q_vectors to the GPU twice even though they
        # are only needed once.  It would mess up modularity quite a bit to
        # handle this optimally, especially since there are many cases where
        # separate q vectors are needed (e.g., form in python and structure
        # in opencl; or both in opencl, but one in single precision and the
        # other in double precision).
        p_kernel = self.P(q_vectors)
        s_kernel = self.S(q_vectors)
        return ProductKernel(self.info, p_kernel, s_kernel)

    def release(self):
        """
        Free resources associated with the model.
        """
        self.P.release()
        self.S.release()


class ProductKernel(object):
    def __init__(self, model_info, p_kernel, s_kernel):
        dim = '2d' if p_kernel.q_input.is_2d else '1d'

        # Need to know if we want 2D and magnetic parameters when constructing
        # a parameter map.
        par_map = {}
        p_info = p_kernel.info['partype']
        s_info = s_kernel.info['partype']
        vol_pars = set(p_info['volume'])
        if dim == '2d':
            num_p_fixed = len(p_info['fixed-2d'])
            num_p_pd = len(p_info['pd-2d'])
            num_s_fixed = len(s_info['fixed-2d'])
            num_s_pd = len(s_info['pd-2d']) - 1 # exclude effect_radius
            # volume parameters are amongst the pd pars for P, not S
            vol_par_idx = [k for k,v in enumerate(p_info['pd-2d'])
                           if v in vol_pars]
        else:
            num_p_fixed = len(p_info['fixed-1d'])
            num_p_pd = len(p_info['pd-1d'])
            num_s_fixed = len(s_info['fixed-1d'])
            num_s_pd = len(s_info['pd-1d']) - 1  # exclude effect_radius
            # volume parameters are amongst the pd pars for P, not S
            vol_par_idx = [k for k,v in enumerate(p_info['pd-1d'])
                           if v in vol_pars]

        start = 0
        par_map['p_fixed'] = np.arange(start, start+num_p_fixed)
        # User doesn't set scale, background or effect_radius for S in P*S,
        # so borrow values from end of p_fixed.  This makes volfraction the
        # first S parameter.
        start += num_p_fixed
        par_map['s_fixed'] = np.hstack(([start,start],
                                        np.arange(start, start+num_s_fixed-2)))
        par_map['volfraction'] = num_p_fixed
        start += num_s_fixed-2
        # vol pars offset from the start of pd pars
        par_map['vol_pars'] = [start+k for k in vol_par_idx]
        par_map['p_pd'] = np.arange(start, start+num_p_pd)
        start += num_p_pd-1
        par_map['s_pd'] = np.hstack((start,
                                     np.arange(start, start+num_s_pd-1)))

        self.fixed_pars = model_info['partype']['fixed-' + dim]
        self.pd_pars = model_info['partype']['pd-' + dim]
        self.info = model_info
        self.p_kernel = p_kernel
        self.s_kernel = s_kernel
        self.par_map = par_map

    def __call__(self, fixed_pars, pd_pars, cutoff=1e-5):
        pars = fixed_pars + pd_pars
        scale = pars[SCALE]
        background = pars[BACKGROUND]
        s_volfraction = pars[self.par_map['volfraction']]
        p_fixed = [pars[k] for k in self.par_map['p_fixed']]
        s_fixed = [pars[k] for k in self.par_map['s_fixed']]
        p_pd = [pars[k] for k in self.par_map['p_pd']]
        s_pd = [pars[k] for k in self.par_map['s_pd']]
        vol_pars = [pars[k] for k in self.par_map['vol_pars']]

        effect_radius, vol_ratio = call_ER_VR(self.p_kernel.info, vol_pars)

        p_fixed[SCALE] = s_volfraction
        p_fixed[BACKGROUND] = 0.0
        s_fixed[SCALE] = scale
        s_fixed[BACKGROUND] = 0.0
        s_fixed[2] = s_volfraction/vol_ratio
        s_pd[0] = [effect_radius], [1.0]

        p_res = self.p_kernel(p_fixed, p_pd)
        s_res = self.s_kernel(s_fixed, s_pd)
        #print s_fixed, s_pd, p_fixed, p_pd

        return p_res*s_res + background

    def release(self):
        self.p_kernel.release()
        self.q_kernel.release()

