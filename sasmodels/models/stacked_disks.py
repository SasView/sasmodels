r"""
This model provides the form factor, $P(q)$, for stacked discs (tactoids)
with a core/layer structure where the form factor is normalized by the volume
of the cylinder. Assuming the next neighbor distance (d-spacing) in a stack
of parallel discs obeys a Gaussian distribution, a structure factor $S(q)$
proposed by Kratky and Porod in 1949 is used in this function.

Note that the resolution smearing calculation uses 76 Gauss quadrature points
to properly smear the model since the function is HIGHLY oscillatory,
especially around the q-values that correspond to the repeat distance of
the layers.

The 2D scattering intensity is the same as 1D, regardless of the orientation
of the q vector which is defined as

.. math:: q = \sqrt{q_x^2 + q_y^2}

Definition
----------

.. figure:: img/stacked_disks_geometry.png

The scattered intensity $I(q)$ is calculated as

.. math::

    I(q) = N\int_{0}^{\pi /2}\left[ \Delta \rho_t
    \left( V_t f_t(q) - V_c f_c(q)\right) + \Delta \rho_c V_c f_c(q)
    \right]^2 S(q)\sin{\alpha}\ d\alpha + \text{background}

where the contrast

.. math::

    \Delta \rho_i = \rho_i - \rho_\text{solvent}

and $N$ is the number of discs per unit volume,
$\alpha$ is the angle between the axis of the disc and $q$,
and $V_t$ and $V_c$ are the total volume and the core volume of
a single disc, respectively.

.. math::

    \left\langle f_{t}^2(q)\right\rangle_{\alpha} =
    \int_{0}^{\pi/2}\left[
    \left(\frac{\sin(q(d+h)\cos{\alpha})}{q(d+h)\cos{\alpha}}\right)
    \left(\frac{2J_1(qR\sin{\alpha})}{qR\sin{\alpha}} \right)
    \right]^2 \sin{\alpha}\ d\alpha

    \left\langle f_{c}^2(q)\right\rangle_{\alpha} =
    \int_{0}^{\pi/2}\left[
    \left(\frac{\sin(qh)\cos{\alpha})}{qh\cos{\alpha}}\right)
    \left(\frac{2J_1(qR\sin{\alpha})}{qR\sin{\alpha}}\right)
    \right]^2 \sin{\alpha}\ d\alpha

where $d$ = thickness of the layer (*thick_layer*),
$2h$ = core thickness (*thick_core*), and $R$ = radius of the disc (*radius*).

.. math::

    S(q) = 1 + \frac{1}{2}\sum_{k=1}^n(n-k)\cos{(kDq\cos{\alpha})}
    \exp\left[ -k(q\cos{\alpha})^2\sigma_d/2\right]

where $n$ is the total number of the disc stacked (*n_stacking*),
$D = 2(d+h)$ is the next neighbor center-to-center distance (d-spacing),
and $\sigma_d$ = the Gaussian standard deviation of the d-spacing (*sigma_d*).

.. note::
    Each assembly in the stack is layer/core/layer, so the spacing of the
    cores is core plus two layers. The 2nd virial coefficient of the cylinder
    is calculated based on the *radius* and *length*
    = *n_stacking* * (*thick_core* + 2 * *thick_layer*)
    values, and used as the effective radius for $S(Q)$ when $P(Q) * S(Q)$
    is applied.

To provide easy access to the orientation of the stacked disks, we define
the axis of the cylinder using two angles $\theta$ and $\varphi$.

.. figure:: img/cylinder_angle_definition.jpg

    Examples of the angles against
    the detector plane.


Our model uses the form factor calculations implemented in a c-library provided
by the NIST Center for Neutron Research (Kline, 2006)

References
----------

A Guinier and G Fournet, *Small-Angle Scattering of X-Rays*,
John Wiley and Sons, New York, 1955

O Kratky and G Porod, *J. Colloid Science*, 4, (1949) 35

J S Higgins and H C Benoit, *Polymers and Neutron Scattering*,
Clarendon, Oxford, 1994

**Author:** NIST IGOR/DANSE **on:** pre 2010

**Last Modified by:** Piotr Rozyczko **on:** February 18, 2016

**Last Reviewed by:** Richard Heenan **on:** March 22, 2016
"""

from numpy import inf

name = "stacked_disks"
title = ""
description = """\
    One layer of disk consists of a core, a top layer, and a bottom layer.
    radius =  the radius of the disk
    thick_core = thickness of the core
    thick_layer = thickness of a layer
    sld_core = the SLD of the core
    sld_layer = the SLD of the layers
    n_stacking = the number of the disks
    sigma_d =  Gaussian STD of d-spacing
    sld_solvent = the SLD of the solvent
    """
category = "shape:cylinder"

# pylint: disable=bad-whitespace, line-too-long
#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["thick_core",  "Ang",        10.0, [0, inf],    "volume",      "Thickness of the core disk"],
    ["thick_layer", "Ang",        10.0, [0, inf],    "volume",      "Thickness of layer each side of core"],
    ["radius",      "Ang",        15.0, [0, inf],    "volume",      "Radius of the stacked disk"],
    ["n_stacking",  "",            1.0, [0, inf],    "volume",      "Number of stacked layer/core/layer disks"],
    ["sigma_d",     "Ang",         0,   [0, inf],    "",            "Sigma of nearest neighbor spacing"],
    ["sld_core",    "1e-6/Ang^2",  4,   [-inf, inf], "sld",         "Core scattering length density"],
    ["sld_layer",   "1e-6/Ang^2",  0.0, [-inf, inf], "sld",         "Layer scattering length density"],
    ["sld_solvent", "1e-6/Ang^2",  5.0, [-inf, inf], "sld",         "Solvent scattering length density"],
    ["theta",       "degrees",     0,   [-inf, inf], "orientation", "Orientation of the stacked disk axis w/respect incoming beam"],
    ["phi",         "degrees",     0,   [-inf, inf], "orientation", "Orientation of the stacked disk in the plane of the detector"],
    ]
# pylint: enable=bad-whitespace, line-too-long

source = ["lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c", "stacked_disks.c"]

demo = dict(background=0.001,
            scale=0.01,
            thick_core=10.0,
            thick_layer=10.0,
            radius=15.0,
            n_stacking=1,
            sigma_d=0,
            sld_core=4,
            sld_layer=0.0,
            sld_solvent=5.0,
            theta=0,
            phi=0)

tests = [
    # Accuracy tests based on content in test/utest_extra_models.py.
    # Added 2 tests with n_stacked = 5 using SasView 3.1.2 - PDB
    [{'thick_core': 10.0,
      'thick_layer': 15.0,
      'radius': 3000.0,
      'n_stacking': 1.0,
      'sigma_d': 0.0,
      'sld_core': 4.0,
      'sld_layer': -0.4,
      'solvent_sd': 5.0,
      'theta': 0.0,
      'phi': 0.0,
      'scale': 0.01,
      'background': 0.001,
     }, 0.001, 5075.12],

    [{'thick_core': 10.0,
      'thick_layer': 15.0,
      'radius': 3000.0,
      'n_stacking': 5.0,
      'sigma_d': 0.0,
      'sld_core': 4.0,
      'sld_layer': -0.4,
      'solvent_sd': 5.0,
      'theta': 0.0,
      'phi': 0.0,
      'scale': 0.01,
      'background': 0.001,
     }, 0.001, 5065.12793824],

    [{'thick_core': 10.0,
      'thick_layer': 15.0,
      'radius': 3000.0,
      'n_stacking': 5.0,
      'sigma_d': 0.0,
      'sld_core': 4.0,
      'sld_layer': -0.4,
      'solvent_sd': 5.0,
      'theta': 0.0,
      'phi': 0.0,
      'scale': 0.01,
      'background': 0.001,
     }, 0.164, 0.0127673597265],

    [{'thick_core': 10.0,
      'thick_layer': 15.0,
      'radius': 3000.0,
      'n_stacking': 1.0,
      'sigma_d': 0.0,
      'sld_core': 4.0,
      'sld_layer': -0.4,
      'solvent_sd': 5.0,
      'theta': 0.0,
      'phi': 0.0,
      'scale': 0.01,
      'background': 0.001,
     }, [0.001, 90.0], [5075.12, 0.001]],

    [{'thick_core': 10.0,
      'thick_layer': 15.0,
      'radius': 3000.0,
      'n_stacking': 1.0,
      'sigma_d': 0.0,
      'sld_core': 4.0,
      'sld_layer': -0.4,
      'solvent_sd': 5.0,
      'theta': 0.0,
      'phi': 0.0,
      'scale': 0.01,
      'background': 0.001,
     }, ([0.4, 0.5]), [0.00105074, 0.00121761]],

    [{'thick_core': 10.0,
      'thick_layer': 15.0,
      'radius': 3000.0,
      'n_stacking': 1.0,
      'sigma_d': 0.0,
      'sld_core': 4.0,
      'sld_layer': -0.4,
      'solvent_sd': 5.0,
      'theta': 0.0,
      'phi': 0.0,
      'scale': 0.01,
      'background': 0.001,
     }, ([1.3, 1.57]), [0.0010039, 0.0010038]],
    ]
# 21Mar2016   RKH notes that unit tests all have n_stacking=1, ought to test other values

