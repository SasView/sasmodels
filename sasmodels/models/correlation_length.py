#correlation length model
# Note: model title and parameter table are inserted automatically
r"""
Definition
----------

The scattering intensity I(q) is calculated as

.. math::
    I(Q) = \frac{A}{Q^n} + \frac{C}{1 + (Q\xi)^m} + B

The first term describes Porod scattering from clusters (exponent = n) and the
second term is a Lorentzian function describing scattering from polymer chains
(exponent = m). This second term characterizes the polymer/solvent interactions
and therefore the thermodynamics. The two multiplicative factors A and C, the
incoherent background B and the two exponents n and m are used as fitting
parameters. (Respectively $porod\_scale$, $lorentz\_scale$, $background$, $exponent\_p$ and 
$exponent\_l$ in the parameter list.) The remaining parameter \ |xi|\  is a correlation 
length for the polymer chains. Note that when m=2 this functional form becomes the 
familiar Lorentzian function. Some interpretation of the values of A and C may be 
possible depending on the values of m and n.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D,
where the q vector is defined as

.. math::
    q = \sqrt{q_x^2 + q_y^2}

References
----------

B Hammouda, D L Ho and S R Kline, Insight into Clustering in
Poly(ethylene oxide) Solutions, Macromolecules, 37 (2004) 6932-6937
"""

from numpy import inf, sqrt

name = "correlation_length"
title = """Calculates an empirical functional form for SAS data characterized
by a low-Q signal and a high-Q signal."""
description = """
"""
category = "shape-independent"
# pylint: disable=bad-continuation, line-too-long
#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [
              ["lorentz_scale", "", 10.0, [0, inf], "", "Lorentzian Scaling Factor"],
              ["porod_scale", "", 1e-06, [0, inf], "", "Porod Scaling Factor"],
              ["cor_length", "Ang", 50.0, [0, inf], "", "Correlation length, xi, in Lorentzian"],
              ["exponent_p", "", 3.0, [0, inf], "", "Porod Exponent, n, in q^-n"],
              ["exponent_l", "1/Ang^2", 2.0, [0, inf], "", "Lorentzian Exponent, m, in 1/( 1 + (q.xi)^m)"],
             ]
# pylint: enable=bad-continuation, line-too-long

def Iq(q, lorentz_scale, porod_scale, cor_length, exponent_p, exponent_l):
    """
    1D calculation of the Correlation length model
    """
    porod = porod_scale / pow(q, exponent_p)
    lorentz = lorentz_scale / (1.0 + pow(q * cor_length, exponent_l))
    inten = porod + lorentz
    return inten

def Iqxy(qx, qy, lorentz_scale, porod_scale, cor_length, exponent_p, exponent_l):
    """
    2D calculation of the Correlation length model
    There is no orientation contribution.
    """
    q = sqrt(qx ** 2 + qy ** 2)
    return Iq(q, lorentz_scale, porod_scale, cor_length, exponent_p, exponent_l)

# parameters for demo
demo = dict(lorentz_scale=10.0, porod_scale=1.0e-06, cor_length=50.0,
            exponent_p=3.0, exponent_l=2.0, background=0.1,
           )

tests = [[{}, 0.001, 1009.98],
         [{}, 0.150141, 0.175645],
         [{}, 0.442528, 0.0213957]]
