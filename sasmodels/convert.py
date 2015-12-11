"""
Convert models to and from sasview.
"""
import warnings

PD_DOT = [
    ("", ""),
    ("_pd", ".width"),
    ("_pd_n", ".npts"),
    ("_pd_nsigma", ".nsigmas"),
    ("_pd_type", ".type"),
    ]
def _convert_pars(pars, mapping):
    """
    Rename the parameters and any associated polydispersity attributes.
    """
    newpars = pars.copy()
    for new, old in mapping.items():
        if old == new: continue
        for pd, dot in PD_DOT:
            if old+dot in newpars:
                if new is not None:
                    newpars[new+pd] = pars[old+dot]
                del newpars[old+dot]
    return newpars

def _rescale_sld(pars):
    """
    rescale all sld parameters in the new model definition by 1e6 so the
    numbers are nicer.  Relies on the fact that all sld parameters in the
    new model definition end with sld.
    """
    return dict((p, (v*1e6 if p.endswith('sld') else v))
                for p,v in pars.items())

def convert_model(name, pars):
    """
    Convert model from old style parameter names to new style.
    """
    _,_ = name,pars # lint
    raise NotImplementedError
    # need to load all new models in order to determine old=>new
    # model name mapping

def _unscale_sld(pars):
    """
    rescale all sld parameters in the new model definition by 1e6 so the
    numbers are nicer.  Relies on the fact that all sld parameters in the
    new model definition end with sld.
    """
    return dict((p, (v*1e-6 if p.endswith('sld') else v))
                for p,v in pars.items())

def _remove_pd(pars, key, name):
    """
    Remove polydispersity from the parameter list.

    Note: operates in place
    """
    # Bumps style parameter names
    pd = pars.pop(key+".width", 0.0)
    pd_n = pars.pop(key+".npts", 0)
    if pd != 0.0 and pd_n != 0:
        warnings.warn("parameter %s not polydisperse in sasview %s"%(key, name))
    pars.pop(key+".nsigmas", None)
    pars.pop(key+".type", None)
    return pars

def _revert_pars(pars, mapping):
    """
    Rename the parameters and any associated polydispersity attributes.
    """
    newpars = pars.copy()

    for new, old in mapping.items():
        for pd, dot in PD_DOT:
            if old+pd == new+dot:
                continue
            if new+pd in newpars:
                if old is not None:
                    newpars[old+dot] = pars[new+pd]
                del newpars[new+pd]
    for k in list(newpars.keys()):
        for pd, dot in PD_DOT[1:]:  # skip "" => ""
            if k.endswith(pd):
                newpars[k[:-len(pd)]+dot] = newpars[k]
                del newpars[k]
    return newpars

def revert_model(model_definition, pars):
    """
    Convert model from new style parameter names to old style.
    """
    mapping = model_definition.oldpars
    oldname = model_definition.oldname
    oldpars = _revert_pars(_unscale_sld(pars), mapping)

    # Note: update compare.constrain_pars to match
    name = model_definition.name
    if name in ('teubner_strey', 'broad_peak'):
        if oldpars.pop('scale', 1.0) != 1.0:
            warnings.warn("parameter scale not used in sasview %s"%name)
    elif name in ('guinier',):
        if oldpars.pop('background', 0.0) != 0.0:
            warnings.warn("parameter background not used in sasview %s"%name)
    elif getattr(model_definition, 'category', None) == 'structure-factor':
        if oldpars.pop('scale', 1.0) != 1.0:
            warnings.warn("parameter scale not used in sasview %s"%name)
        if oldpars.pop('background', 0.0) != 0.0:
            warnings.warn("parameter background not used in sasview %s"%name)
    elif name == 'pearl_necklace':
        _remove_pd(oldpars, 'num_pearls', name)
        _remove_pd(oldpars, 'thick_string', name)

    return oldname, oldpars

def constrain_new_to_old(model_definition, pars):
    """
    Restrict parameter values to those that will match sasview.
    """
    # Note: update convert.revert_model to match
    name = model_definition.name
    if name in ('teubner_strey','broad_peak'):
        pars['scale'] = 1
    elif name in ('guinier',):
        pars['background'] = 0
    elif name == 'pearl_necklace':
        pars['string_thickness_pd_n'] = 0
        pars['number_of_pearls_pd_n'] = 0
    elif getattr(model_definition, 'category', None) == 'structure-factor':
        pars['scale'], pars['background'] = 1, 0

