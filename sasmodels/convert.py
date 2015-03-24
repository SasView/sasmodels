"""
Convert models to and from sasview.
"""

def _rename_pars(pars, mapping):
    """
    Rename the parameters and any associated polydispersity attributes.
    """
    newpars = pars.copy()
    for old,new in mapping.items():
        if old == new: continue
        # Bumps style parameter names
        for variant in ("", "_pd", "_pd_n", "_pd_nsigma", "_pd_type"):
            if old+variant in newpars:
                if new is not None:
                    newpars[new+variant] = pars[old+variant]
                del newpars[old+variant]
        # Sasview style parameter names
        for variant in (".width", ".nsigmas", ".type", ".npts"):
            if old+variant in newpars:
                if new is not None:
                    newpars[new+variant] = pars[old+variant]
                del newpars[old+variant]
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

def revert_model(model_definition, pars):
    """
    Convert model from new style parameter names to old style.
    """
    mapping = model_definition.oldpars
    oldname = model_definition.oldname
    oldpars = _rename_pars(_unscale_sld(pars), mapping)
    return oldname, oldpars

