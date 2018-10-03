"""
Convert models to and from sasview.
"""
from __future__ import print_function, division

import math
import warnings

import numpy as np

from .conversion_table import CONVERSION_TABLE
from .core import load_model_info

# List of models which SasView versions don't contain the explicit 'scale' argument.
# When converting such a model, please update this list.
MODELS_WITHOUT_SCALE = [
    'teubner_strey',
    'broad_peak',
    'two_lorentzian',
    "two_power_law",
    'gauss_lorentz_gel',
    'be_polyelectrolyte',
    'correlation_length',
    'fractal_core_shell',
    'binary_hard_sphere',
    'raspberry'
]

# List of models which SasView versions don't contain the explicit 'background' argument.
# When converting such a model, please update this list.
MODELS_WITHOUT_BACKGROUND = [
    'guinier',
]

MODELS_WITHOUT_VOLFRACTION = [
    'fractal',
    'vesicle',
    'multilayer_vesicle',
]

MAGNETIC_SASVIEW_MODELS = [
    'core_shell',
    'core_multi_shell',
    'cylinder',
    'parallelepiped',
    'sphere',
]


# Convert new style names for polydispersity info to old style names
PD_DOT = [
    ("_pd", ".width"),
    ("_pd_n", ".npts"),
    ("_pd_nsigma", ".nsigmas"),
    ("_pd_type", ".type"),
    (".lower", ".lower"),
    (".upper", ".upper"),
    (".fittable", ".fittable"),
    (".std", ".std"),
    (".units", ".units"),
    ("", "")
    ]

def _rescale(par, scale):
    return [pk*scale for pk in par] if isinstance(par, list) else par*scale

def _is_sld(model_info, par):
    """
    Return True if parameter is a magnetic magnitude or SLD parameter.
    """
    if par.startswith('M0:'):
        return True
    if '_pd' in par or '.' in par:
        return False
    for p in model_info.parameters.call_parameters:
        if p.id == par:
            return p.type == 'sld'
    # check through kernel parameters in case it is a named as a vector
    for p in model_info.parameters.kernel_parameters:
        if p.id == par:
            return p.type == 'sld'
    return False

def _rescale_sld(model_info, pars, scale):
    """
    rescale all sld parameters in the new model definition by *scale* so the
    numbers are nicer.  Relies on the fact that all sld parameters in the
    new model definition end with sld.  For backward conversion use
    *scale=1e-6*.  For forward conversion use *scale=1e6*.
    """
    return dict((par, (_rescale(v, scale) if _is_sld(model_info, par) else v))
                for par, v in pars.items())


def _get_translation_table(model_info, version=(3, 1, 2)):
    conv_param = CONVERSION_TABLE.get(version, {}).get(model_info.id, [None, {}])
    translation = conv_param[1].copy()
    for p in model_info.parameters.kernel_parameters:
        if p.length > 1:
            newid = p.id
            oldid = translation.get(p.id, p.id)
            translation.pop(newid, None)
            for k in range(1, p.length+1):
                if newid+str(k) not in translation:
                    translation[newid+str(k)] = oldid+str(k)
    # Remove control parameter from the result
    if model_info.control:
        translation[model_info.control] = "CONTROL"
    return translation

# ========= FORWARD CONVERSION sasview 3.x => sasmodels ===========
def _dot_pd_to_underscore_pd(par):
    if par.endswith(".width"):
        return par[:-6]+"_pd"
    elif par.endswith(".type"):
        return par[:-5]+"_pd_type"
    elif par.endswith(".nsigmas"):
        return par[:-8]+"_pd_nsigma"
    elif par.endswith(".npts"):
        return par[:-5]+"_pd_n"
    else:
        return par

def _pd_to_underscores(pars):
    return dict((_dot_pd_to_underscore_pd(k), v) for k, v in pars.items())

def _convert_pars(pars, mapping):
    """
    Rename the parameters and any associated polydispersity attributes.
    """
    newpars = pars.copy()
    for new, old in mapping.items():
        if old == new:
            continue
        if old is None:
            continue
        for _, dot in PD_DOT:
            source = old+dot
            if source in newpars:
                if new is not None:
                    target = new+dot
                else:
                    target = None
                if source != target:
                    if target:
                        newpars[target] = pars[old+dot]
                    del newpars[source]
    return newpars

def _conversion_target(model_name, version=(3, 1, 2)):
    """
    Find the sasmodel name which translates into the sasview name.

    Note: *CoreShellEllipsoidModel* translates into *core_shell_ellipsoid:1*.
    This is necessary since there is only one variant in sasmodels for the
    two variants in sasview.
    """
    for sasmodels_name, sasview_dict in \
            CONVERSION_TABLE.get(version, {}).items():
        if sasview_dict[0] == model_name:
            return sasmodels_name
    return None

def _hand_convert(name, oldpars, version=(3, 1, 2)):
    if version == (3, 1, 2):
        oldpars = _hand_convert_3_1_2_to_4_1(name, oldpars)
    if version < (4, 2, 0):
        oldpars = _rename_magnetic_pars(oldpars)
    return oldpars

def _rename_magnetic_pars(pars):
    """
    Change from M0:par to par_M0, etc.
    """
    keys = list(pars.items())
    for k in keys:
        if k.startswith('M0:'):
            pars[k[3:]+'_M0'] = pars.pop(k)
        elif k.startswith('mtheta:'):
            pars[k[7:]+'_mtheta'] = pars.pop(k)
        elif k.startswith('mphi:'):
            pars[k[5:]+'_mphi'] = pars.pop(k)
        elif k.startswith('up:'):
            pars['up_'+k[3:]] = pars.pop(k)
    return pars

def _hand_convert_3_1_2_to_4_1(name, oldpars):
    if name == 'core_shell_parallelepiped':
        # Make sure pd on rim parameters defaults to zero
        # ... probably not necessary.
        oldpars['rimA.width'] = 0.0
        oldpars['rimB.width'] = 0.0
        oldpars['rimC.width'] = 0.0
    elif name == 'core_shell_ellipsoid:1':
        # Reverse translation (from new to old), from core_shell_ellipsoid.c
        #    equat_shell = equat_core + thick_shell
        #    polar_core = equat_core * x_core
        #    polar_shell = equat_core * x_core + thick_shell*x_polar_shell
        # Forward translation (from old to new), inverting reverse translation:
        #    thick_shell = equat_shell - equat_core
        #    x_core = polar_core / equat_core
        #    x_polar_shell = (polar_shell - polar_core)/(equat_shell - equat_core)
        # Auto translation (old <=> new) happens after hand_convert
        #    equat_shell <=> thick_shell
        #    polar_core <=> x_core
        #    polar_shell <=> x_polar_shell
        # So...
        equat_core, equat_shell = oldpars['equat_core'], oldpars['equat_shell']
        polar_core, polar_shell = oldpars['polar_core'], oldpars['polar_shell']
        oldpars['equat_shell'] = equat_shell - equat_core
        oldpars['polar_core'] = polar_core / equat_core
        oldpars['polar_shell'] = (polar_shell-polar_core)/(equat_shell-equat_core)
    elif name == 'hollow_cylinder':
        # now uses radius and thickness
        thickness = oldpars['radius'] - oldpars['core_radius']
        oldpars['radius'] = thickness
        if 'radius.width' in oldpars:
            pd = oldpars['radius.width']*oldpars['radius']/thickness
            oldpars['radius.width'] = pd
    elif name == 'multilayer_vesicle':
        if 'scale' in oldpars:
            oldpars['volfraction'] = oldpars['scale']
            oldpars['scale'] = 1.0
        if 'scale.lower' in oldpars:
            oldpars['volfraction.lower'] = oldpars['scale.lower']
        if 'scale.upper' in oldpars:
            oldpars['volfraction.upper'] = oldpars['scale.upper']
        if 'scale.fittable' in oldpars:
            oldpars['volfraction.fittable'] = oldpars['scale.fittable']
        if 'scale.std' in oldpars:
            oldpars['volfraction.std'] = oldpars['scale.std']
        if 'scale.units' in oldpars:
            oldpars['volfraction.units'] = oldpars['scale.units']
    elif name == 'pearl_necklace':
        pass
        #_remove_pd(oldpars, 'num_pearls', name)
        #_remove_pd(oldpars, 'thick_string', name)
    elif name == 'polymer_micelle':
        if 'ndensity' in oldpars:
            oldpars['ndensity'] /= 1e15
        if 'ndensity.lower' in oldpars:
            oldpars['ndensity.lower'] /= 1e15
        if 'ndensity.upper' in oldpars:
            oldpars['ndensity.upper'] /= 1e15
    elif name == 'rpa':
        # convert scattering lengths from femtometers to centimeters
        for p in "L1", "L2", "L3", "L4":
            if p in oldpars:
                oldpars[p] /= 1e-13
            if p + ".lower" in oldpars:
                oldpars[p + ".lower"] /= 1e-13
            if p + ".upper" in oldpars:
                oldpars[p + ".upper"] /= 1e-13
    elif name == 'spherical_sld':
        j = 0
        while "func_inter" + str(j) in oldpars:
            name = "func_inter" + str(j)
            new_name = "shape" + str(j + 1)
            if oldpars[name] == 'Erf(|nu|*z)':
                oldpars[new_name] = int(0)
            elif oldpars[name] == 'RPower(z^|nu|)':
                oldpars[new_name] = int(1)
            elif oldpars[name] == 'LPower(z^|nu|)':
                oldpars[new_name] = int(2)
            elif oldpars[name] == 'RExp(-|nu|*z)':
                oldpars[new_name] = int(3)
            elif oldpars[name] == 'LExp(-|nu|*z)':
                oldpars[new_name] = int(4)
            else:
                oldpars[new_name] = int(0)
            oldpars.pop(name)
            oldpars['n_shells'] = str(j + 1)
            j += 1
    elif name == 'teubner_strey':
        # basically undoing the entire Teubner-Strey calculations here.
        #    drho = (sld_a - sld_b)
        #    k = 2.0*math.pi*xi/d
        #    a2 = (1.0 + k**2)**2
        #    c1 = 2.0 * xi**2 * (1.0 - k**2)
        #    c2 = xi**4
        #    prefactor = 8.0*math.pi*phi*(1.0-phi)*drho**2*c2/xi
        #    scale = 1e-4*prefactor
        #    oldpars['scale'] = a2/scale
        #    oldpars['c1'] = c1/scale
        #    oldpars['c2'] = c2/scale

        # need xi, d, sld_a, sld_b, phi=volfraction_a
        # assume contrast is 1.0e-6, scale=1, background=0
        sld_a, sld_b = 1.0, 0.
        drho = sld_a - sld_b

        # find xi
        p_scale = oldpars['scale']
        p_c1 = oldpars['c1']
        p_c2 = oldpars['c2']
        i_1 = 0.5*p_c1/p_c2
        i_2 = math.sqrt(math.fabs(p_scale/p_c2))
        i_3 = 2/(i_1 + i_2)
        xi = math.sqrt(math.fabs(i_3))

        # find d from xi
        k = math.sqrt(math.fabs(1 - 0.5*p_c1/p_c2*xi**2))
        d = 2*math.pi*xi/k

        # solve quadratic phi (1-phi) = xi/(1e-4 8 pi drho^2 c2)
        # favour volume fraction in [0, 0.5]
        c = xi / (1e-4 * 8.0 * math.pi * drho**2 * p_c2)
        phi = 0.5 - math.sqrt(0.25 - c)

        # scale sld_a by 1e-6 because the translator will scale it back
        oldpars.update(volfraction_a=phi, xi=xi, d=d, sld_a=sld_a*1e-6,
                       sld_b=sld_b, scale=1.0)
        oldpars.pop('c1')
        oldpars.pop('c2')

    return oldpars

def convert_model(name, pars, use_underscore=False, model_version=(3, 1, 2)):
    """
    Convert model from old style parameter names to new style.
    """
    newpars = pars
    keys = sorted(CONVERSION_TABLE.keys())
    for i, version in enumerate(keys):
        # Don't allow indices outside list
        next_i = i + 1
        if next_i == len(keys):
            next_i = i
        # If the save state is from a later version, skip the check
        if model_version <= keys[next_i]:
            newname = _conversion_target(name, version)
        else:
            newname = None
        # If no conversion is found, move on
        if newname is None:
            newname = name
            continue
        if ':' in newname:   # core_shell_ellipsoid:1
            model_info = load_model_info(newname[:-2])
            # Know the table exists and isn't multiplicity so grab it directly
            # Can't use _get_translation_table since that will return the 'bare'
            # version.
            translation = CONVERSION_TABLE.get(version, {})[newname][1]
        else:
            model_info = load_model_info(newname)
            translation = _get_translation_table(model_info, version)
        newpars = _hand_convert(newname, newpars, version)
        newpars = _convert_pars(newpars, translation)
        # TODO: Still not convinced this is the best check
        if not model_info.structure_factor and version == (3, 1, 2):
            newpars = _rescale_sld(model_info, newpars, 1e6)
        newpars.setdefault('scale', 1.0)
        newpars.setdefault('background', 0.0)
        if use_underscore:
            newpars = _pd_to_underscores(newpars)
        name = newname
    return newname, newpars

# ========= BACKWARD CONVERSION sasmodels => sasview 3.x ===========

def _revert_pars(pars, mapping):
    """
    Rename the parameters and any associated polydispersity attributes.
    """
    newpars = pars.copy()

    for new, old in mapping.items():
        for underscore, dot in PD_DOT:
            if old and old+underscore == new+dot:
                continue
            if new+underscore in newpars:
                if old is not None:
                    newpars[old+dot] = pars[new+underscore]
                del newpars[new+underscore]
    for k in list(newpars.keys()):
        for underscore, dot in PD_DOT[1:]:  # skip "" => ""
            if k.endswith(underscore):
                newpars[k[:-len(underscore)]+dot] = newpars[k]
                del newpars[k]
    return newpars

def revert_name(model_info):
    oldname, _ = CONVERSION_TABLE.get(model_info.id, [None, {}])
    return oldname

def _remove_pd(pars, key, name):
    """
    Remove polydispersity from the parameter list.

    Note: operates in place
    """
    # Bumps style parameter names
    width = pars.pop(key+".width", 0.0)
    n_points = pars.pop(key+".npts", 0)
    if width != 0.0 and n_points != 0:
        warnings.warn("parameter %s not polydisperse in sasview %s"%(key, name))
    pars.pop(key+".nsigmas", None)
    pars.pop(key+".type", None)
    return pars

def _trim_vectors(model_info, pars, oldpars):
    _, translation = CONVERSION_TABLE.get(model_info.id, [None, {}])
    for p in model_info.parameters.kernel_parameters:
        if p.length_control is not None:
            n = int(pars[p.length_control])
            oldname = translation.get(p.id, p.id)
            for k in range(n+1, p.length+1):
                for _, old in PD_DOT:
                    oldpars.pop(oldname+str(k)+old, None)
    return oldpars

def revert_pars(model_info, pars):
    """
    Convert model from new style parameter names to old style.
    """
    if model_info.composition is not None:
        composition_type, parts = model_info.composition
        if composition_type == 'product':
            translation = _get_translation_table(parts[0])
            # structure factor models include scale:scale_factor mapping
            translation.update(_get_translation_table(parts[1]))
        else:
            raise NotImplementedError("cannot convert to sasview sum")
    else:
        translation = _get_translation_table(model_info)
    oldpars = _revert_pars(_rescale_sld(model_info, pars, 1e-6), translation)
    oldpars = _trim_vectors(model_info, pars, oldpars)

    # Make sure the control parameter is an integer
    if "CONTROL" in oldpars:
        oldpars["CONTROL"] = int(oldpars["CONTROL"])

    # Note: update compare.constrain_pars to match
    name = model_info.id
    if name in MODELS_WITHOUT_SCALE or model_info.structure_factor:
        if oldpars.pop('scale', 1.0) != 1.0:
            warnings.warn("parameter scale not used in sasview %s"%name)
    if name in MODELS_WITHOUT_BACKGROUND or model_info.structure_factor:
        if oldpars.pop('background', 0.0) != 0.0:
            warnings.warn("parameter background not used in sasview %s"%name)

    # Remove magnetic parameters from non-magnetic sasview models
    if name not in MAGNETIC_SASVIEW_MODELS:
        oldpars = dict((k, v) for k, v in oldpars.items() if ':' not in k)

    # If it is a product model P*S, then check the individual forms for special
    # cases.  Note: despite the structure factor alone not having scale or
    # background, the product model does, so this is below the test for
    # models without scale or background.
    namelist = name.split('*') if '*' in name else [name]
    for name in namelist:
        if name in MODELS_WITHOUT_VOLFRACTION:
            del oldpars['volfraction']
        elif name == 'core_multi_shell':
            # kill extra shells
            for k in range(5, 11):
                oldpars.pop('sld_shell'+str(k), 0)
                oldpars.pop('thick_shell'+str(k), 0)
                oldpars.pop('mtheta:sld'+str(k), 0)
                oldpars.pop('mphi:sld'+str(k), 0)
                oldpars.pop('M0:sld'+str(k), 0)
                _remove_pd(oldpars, 'sld_shell'+str(k), 'sld')
                _remove_pd(oldpars, 'thick_shell'+str(k), 'thickness')
        elif name == 'core_shell_parallelepiped':
            _remove_pd(oldpars, 'rimA', name)
            _remove_pd(oldpars, 'rimB', name)
            _remove_pd(oldpars, 'rimC', name)
        elif name == 'hollow_cylinder':
            # now uses radius and thickness
            thickness = oldpars['core_radius']
            oldpars['radius'] += thickness
            oldpars['radius.width'] *= thickness/oldpars['radius']
        #elif name in ['mono_gauss_coil', 'poly_gauss_coil']:
        #    del oldpars['i_zero']
        elif name == 'onion':
            oldpars.pop('n_shells', None)
        elif name == 'pearl_necklace':
            _remove_pd(oldpars, 'num_pearls', name)
            _remove_pd(oldpars, 'thick_string', name)
        elif name == 'polymer_micelle':
            if 'ndensity' in oldpars:
                oldpars['ndensity'] *= 1e15
        elif name == 'rpa':
            # convert scattering lengths from femtometers to centimeters
            for p in "L1", "L2", "L3", "L4":
                if p in oldpars: oldpars[p] *= 1e-13
            if pars['case_num'] < 2:
                for k in ("a", "b"):
                    for p in ("L", "N", "Phi", "b", "v"):
                        oldpars.pop(p+k, None)
                for k in "Kab,Kac,Kad,Kbc,Kbd".split(','):
                    oldpars.pop(k, None)
            elif pars['case_num'] < 5:
                for k in ("a",):
                    for p in ("L", "N", "Phi", "b", "v"):
                        oldpars.pop(p+k, None)
                for k in "Kab,Kac,Kad".split(','):
                    oldpars.pop(k, None)
        elif name == 'spherical_sld':
            oldpars["CONTROL"] -= 1
            # remove polydispersity from shells
            for k in range(1, 11):
                _remove_pd(oldpars, 'thick_flat'+str(k), 'thickness')
                _remove_pd(oldpars, 'thick_inter'+str(k), 'interface')
            # remove extra shells
            for k in range(int(pars['n_shells']), 11):
                oldpars.pop('sld_flat'+str(k), 0)
                oldpars.pop('thick_flat'+str(k), 0)
                oldpars.pop('thick_inter'+str(k), 0)
                oldpars.pop('func_inter'+str(k), 0)
                oldpars.pop('nu_inter'+str(k), 0)
        elif name == 'stacked_disks':
            _remove_pd(oldpars, 'n_stacking', name)
        elif name == 'teubner_strey':
            # basically redoing the entire Teubner-Strey calculations here.
            volfraction = oldpars.pop('volfraction_a')
            xi = oldpars.pop('xi')
            d = oldpars.pop('d')
            sld_a = oldpars.pop('sld_a')
            sld_b = oldpars.pop('sld_b')
            drho = 1e6*(sld_a - sld_b)  # conversion autoscaled these
            k = 2.0*math.pi*xi/d
            a2 = (1.0 + k**2)**2
            c1 = 2.0 * xi**2 * (1.0 - k**2)
            c2 = xi**4
            prefactor = 8.0*math.pi*volfraction*(1.0-volfraction)*drho**2*c2/xi
            scale = 1e-4*prefactor
            oldpars['scale'] = a2/scale
            oldpars['c1'] = c1/scale
            oldpars['c2'] = c2/scale

    #print("convert from",list(sorted(pars)))
    #print("convert to",list(sorted(oldpars.items())))
    return oldpars

def constrain_new_to_old(model_info, pars):
    """
    Restrict parameter values to those that will match sasview.
    """
    name = model_info.id
    # Note: update convert.revert_model to match
    if name in MODELS_WITHOUT_SCALE or model_info.structure_factor:
        pars['scale'] = 1
    if name in MODELS_WITHOUT_BACKGROUND or model_info.structure_factor:
        pars['background'] = 0
    # sasview multiplies background by structure factor
    if '*' in name:
        pars['background'] = 0

    # Shut off magnetism when comparing non-magnetic sasview models
    if name not in MAGNETIC_SASVIEW_MODELS:
        suppress_magnetism = False
        for key in pars.keys():
            if key.startswith("M0:"):
                suppress_magnetism = suppress_magnetism or (pars[key] != 0)
                pars[key] = 0
        if suppress_magnetism:
            warnings.warn("suppressing magnetism for comparison with sasview")

    # Shut off theta polydispersity since algorithm has changed
    if 'theta_pd_n' in pars:
        if pars['theta_pd_n'] != 0:
            warnings.warn("suppressing theta polydispersity for comparison with sasview")
        pars['theta_pd_n'] = 0

    # If it is a product model P*S, then check the individual forms for special
    # cases.  Note: despite the structure factor alone not having scale or
    # background, the product model does, so this is below the test for
    # models without scale or background.
    namelist = name.split('*') if '*' in name else [name]
    for name in namelist:
        if name in MODELS_WITHOUT_VOLFRACTION:
            pars['volfraction'] = 1
        if name == 'core_multi_shell':
            pars['n'] = min(math.ceil(pars['n']), 4)
        elif name == 'gel_fit':
            pars['scale'] = 1
        elif name == 'line':
            pars['scale'] = 1
            pars['background'] = 0
        elif name == 'mono_gauss_coil':
            pars['scale'] = 1
        elif name == 'onion':
            pars['n_shells'] = math.ceil(pars['n_shells'])
        elif name == 'pearl_necklace':
            pars['string_thickness_pd_n'] = 0
            pars['number_of_pearls_pd_n'] = 0
        elif name == 'poly_gauss_coil':
            pars['scale'] = 1
        elif name == 'rpa':
            pars['case_num'] = int(pars['case_num'])
        elif name == 'spherical_sld':
            pars['n_shells'] = math.ceil(pars['n_shells'])
            pars['n_steps'] = math.ceil(pars['n_steps'])
            for k in range(1, 11):
                pars['shape%d'%k] = math.trunc(pars['shape%d'%k]+0.5)
            for k in range(2, 11):
                pars['thickness%d_pd_n'%k] = 0
                pars['interface%d_pd_n'%k] = 0
        elif name == 'teubner_strey':
            pars['scale'] = 1
            if pars['volfraction_a'] > 0.5:
                pars['volfraction_a'] = 1.0 - pars['volfraction_a']
        elif name == 'unified_power_Rg':
            pars['level'] = int(pars['level'])

def _check_one(name, seed=None):
    """
    Generate a random set of parameters for *name*, and check that they can
    be converted back to SasView 3.x and forward again to sasmodels.  Raises
    an error if the parameters are changed.
    """
    from . import compare

    model_info = load_model_info(name)

    old_name = revert_name(model_info)
    if old_name is None:
        return

    pars = compare.get_pars(model_info, use_demo=False)
    if seed is not None:
        np.random.seed(seed)
    pars = compare.randomize_pars(model_info, pars)
    if name == "teubner_strey":
        # T-S model is underconstrained, so fix the assumptions.
        pars['sld_a'], pars['sld_b'] = 1.0, 0.0
    compare.constrain_pars(model_info, pars)
    constrain_new_to_old(model_info, pars)
    old_pars = revert_pars(model_info, pars)
    new_name, new_pars = convert_model(old_name, old_pars, use_underscore=True)
    if 1:
        print("==== %s in ====="%name)
        print(str(compare.parlist(model_info, pars, True)))
        print("==== %s ====="%old_name)
        for k, v in sorted(old_pars.items()):
            print(k, v)
        print("==== %s out ====="%new_name)
        print(str(compare.parlist(model_info, new_pars, True)))
    assert name == new_name, "%r != %r"%(name, new_name)
    for k, v in new_pars.items():
        assert k in pars, "%s: %r appeared from conversion"%(name, k)
        if isinstance(v, float):
            assert abs(v-pars[k]) <= abs(1e-12*v), \
                "%s: %r  %s != %s"%(name, k, v, pars[k])
        else:
            assert v == pars[k], "%s: %r  %s != %s"%(name, k, v, pars[k])
    for k, v in pars.items():
        assert k in pars, "%s: %r not converted"%(name, k)

def test_backward_forward():
    from .core import list_models
    L = lambda name: _check_one(name, seed=1)
    for name in list_models('all'):
        yield L, name
