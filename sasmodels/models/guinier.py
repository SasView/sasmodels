r"""
Definition
----------

This model fits the Guinier function

.. math::

    I(q) = \text{scale} \cdot \exp{\left[ \frac{-Q^2R_g^2}{3} \right]}
            + \text{background}

to the data directly without any need for linearisation (*cf*. the usual
plot of $\ln I(q)$ vs $q^2$\ ). Note that you may have to restrict the data
range to include small q only, where the Guinier approximation actually
applies. See also the guinier_porod model.

For 2D data the scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math:: q = \sqrt{q_x^2 + q_y^2}

References
----------

A Guinier and G Fournet, *Small-Angle Scattering of X-Rays*,
John Wiley & Sons, New York (1955)
"""

from numpy import inf

name = "guinier"
title = ""
description = """
 I(q) = scale.exp ( - rg^2 q^2 / 3.0 )

    List of default parameters:
    scale = scale
    rg = Radius of gyration
"""
category = "shape-independent"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["rg", "Ang", 60.0, [0, inf], "", "Radius of Gyration"]]

Iq = """
    double exponent = rg*rg*q*q/3.0;
    double value = exp(-exponent);
    return value;
"""

def random():
    import numpy as np
    scale = 10**np.random.uniform(1, 5)
    # Note: compare.py has rg cutoff for guinier, so use that
    q_max = 1.0
    rg_max = np.sqrt(90*np.log(10) + 3*np.log(scale))/q_max
    rg = 10**np.random.uniform(0, np.log10(rg_max))
    pars = dict(
        #background=0,
        scale=scale,
        rg=rg,
    )
    return pars

# parameters for demo
demo = dict(scale=1.0, rg=60.0)

# parameters for unit tests
tests = [[{'rg' : 31.5}, 0.005, 0.992756]]
