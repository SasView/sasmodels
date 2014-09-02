"""
Convert models to and from sasview.
"""
PARAMETER_MAP = {
    'CylinderModel': dict(
        name='cylinder',
        cyl_theta='theta', cyl_phi='phi',
        sldCyl='sld', sldSolv='solvent_sld',
        ),
    'EllipsoidModel': dict(
        name='ellipsoid',
        axis_theta='theta', axis_phi='phi',
        sldEll='sld', sldSolv='solvent_sld',
        radius_a='rpolar', radius_b='requatorial',
        ),
    'CoreShellCylinderModel': dict(
        name='core_shell_cylinder',
        axis_theta='theta', axis_phi='phi',
        ),
    'TriaxialEllipsoidModel': dict(
        name='triaxial_ellipsoid',
        axis_theta='theta', axis_phi='phi', axis_psi='psi',
        sldEll='sld', sldSolv='solvent_sld',
        semi_axisA='req_minor', semi_axisB='req_major', semi_axisC='rpolar',
        ),
    'LamellarModel': dict(
        name='lamellar',
        sld_bi='sld', sld_sol='solvent_sld',
        bi_thick='thickness',
        ),
    'CappedCylinderModel': dict(
        name='capped_cylinder',
        sld_capcyl='sld', sld_solv='solvent_sld',
        len_cyl='length', rad_cyl='radius', rad_cap='cap_radius',
        ),
    'SphereModel': dict(
        name='sphere',
        sldSph='sld', sldSolv='solvent_sld',
        radius='radius',  # listing identical parameters is optional
        ),
    }

def _reverse_map():
    retval = {}
    for old_name,old_map in PARAMETER_MAP.items():
        new_name = old_map['name']
        new_map = dict((v,k) for k,v in old_map.items() if k != 'name')
        new_map['name'] = old_name
        retval[new_name] = new_map
    return retval
REVERSE_MAP = _reverse_map()
del _reverse_map


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
                newpars[new+variant] = pars[old+variant]
                del newpars[old+variant]
        # Sasview style parameter names
        for variant in (".width", ".nsigmas", ".type", ".npts"):
            if old+variant in newpars:
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
    mapping = PARAMETER_MAP[name]
    newname = mapping['name']
    newpars = _rescale_sld(_rename_pars(pars, mapping))
    return newname, newpars

def _unscale_sld(pars):
    """
    rescale all sld parameters in the new model definition by 1e6 so the
    numbers are nicer.  Relies on the fact that all sld parameters in the
    new model definition end with sld.
    """
    return dict((p, (v*1e-6 if p.endswith('sld') else v))
                for p,v in pars.items())

def revert_model(name, pars):
    """
    Convert model from new style parameter names to old style.
    """
    mapping = REVERSE_MAP[name]
    oldname = mapping['name']
    oldpars = _rename_pars(_unscale_sld(pars), mapping)
    return oldname, oldpars

