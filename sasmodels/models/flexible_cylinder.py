r"""
This model provides the form factor, $P(q)$, for a flexible cylinder
where the form factor is normalized by the volume of the cylinder.
**Inter-cylinder interactions are NOT provided for.**

.. math::

    P(q) = \text{scale} \left<F^2\right>/V + \text{background}

where the averaging $\left<\ldots\right>$ is applied only for the 1D
calculation

The 2D scattering intensity is the same as 1D, regardless of the orientation of
the q vector which is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}

Definitions
-----------

.. figure:: img/flexible_cylinder_geometry.jpg


The chain of contour length, $L$, (the total length) can be described as a
chain of some number of locally stiff segments of length $l_p$, the persistence
length (the length along the cylinder over which the flexible cylinder can be
considered a rigid rod).
The Kuhn length $(b = 2*l_p)$ is also used to describe the stiffness of a chain.

The returned value is in units of $cm^{-1}$, on absolute scale.

In the parameters, the sldCyl and sldSolv represent the SLD of the chain/cylinder
and solvent respectively.

Our model uses the form factor calculations implemented in a c-library provided
by the NIST Center for Neutron Research (Kline, 2006).


From the reference:

    'Method 3 With Excluded Volume' is used.
    The model is a parametrization of simulations of a discrete representation
    of the worm-like chain model of Kratky and Porod applied in the
    pseudocontinuous limit.
    See equations (13,26-27) in the original reference for the details.

References
----------

J S Pedersen and P Schurtenberger. *Scattering functions of semiflexible
polymers with and without excluded volume effects.* Macromolecules,
29 (1996) 7602-7612

Correction of the formula can be found in

W R Chen, P D Butler and L J Magid, *Incorporating Intermicellar Interactions
in the Fitting of SANS Data from Cationic Wormlike Micelles.* Langmuir,
22(15) 2006 6539-6548
"""
from numpy import inf

name = "flexible_cylinder"
title = "Flexible cylinder where the form factor is normalized by the volume" \
        "of the cylinder."
description = """Note : scale and contrast = (sld - sld_solvent) are both
                multiplicative factors in the model and are perfectly
                correlated. One or both of these parameters must be held fixed
                during model fitting.
              """

category = "shape:cylinder"
single = False  # double precision only!

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["length",      "Ang",       1000.0, [0, inf],    "volume", "Length of the flexible cylinder"],
    ["kuhn_length", "Ang",        100.0, [0, inf],    "volume", "Kuhn length of the flexible cylinder"],
    ["radius",      "Ang",         20.0, [0, inf],    "volume", "Radius of the flexible cylinder"],
    ["sld",         "1e-6/Ang^2",   1.0, [-inf, inf], "",       "Cylinder scattering length density"],
    ["sld_solvent", "1e-6/Ang^2",   6.3, [-inf, inf], "",       "Solvent scattering length density"],
    ]
# pylint: enable=bad-whitespace, line-too-long
source = ["lib/J1.c", "lib/J1c.c", "lib/wrc_cyl.c", "flexible_cylinder.c"]

demo = dict(scale=1.0, background=0.0001,
            length=1000.0,
            kuhn_length=100.0,
            radius=20.0,
            sld=1.0,
            sld_solvent=6.3)

oldname = 'FlexibleCylinderModel'
oldpars = dict(sld='sldCyl', sld_solvent='sldSolv')


tests = [
    # Accuracy tests based on content in test/utest_other_models.py
    # Currently fails in OCL
    # [{'length':     1000.0,
    #  'kuhn_length': 100.0,
    #  'radius':       20.0,
    #  'sld':           1.0,
    #  'sld_solvent':   6.3,
    #  'background':    0.0001,
    #  }, 0.001, 3509.2187],

    # Additional tests with larger range of parameters
    [{'length':    1000.0,
      'kuhn_length': 100.0,
      'radius':       20.0,
      'sld':           1.0,
      'sld_solvent':   6.3,
      'background':    0.0001,
     }, 1.0, 0.000595345],
    [{'length':        10.0,
      'kuhn_length': 800.0,
      'radius':        2.0,
      'sld':           6.0,
      'sld_solvent':  12.3,
      'background':    0.001,
     }, 0.1, 1.55228],
    [{'length':        100.0,
      'kuhn_length': 800.0,
      'radius':       50.0,
      'sld':           0.1,
      'sld_solvent':   5.1,
      'background':    0.0,
     }, 1.0, 0.000938456]
    ]

