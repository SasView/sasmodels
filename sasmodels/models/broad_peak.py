r"""
This model calculates an empirical functional form for SAS data characterized
by a broad scattering peak. Many SAS spectra are characterized by a broad peak
even though they are from amorphous soft materials. For example, soft systems
that show a SAS peak include copolymers, polyelectrolytes, multiphase systems,
layered structures, etc.

The d-spacing corresponding to the broad peak is a characteristic distance
between the scattering inhomogeneities (such as in lamellar, cylindrical, or
spherical morphologies, or for bicontinuous structures).

The returned value is scaled to units of |cm^-1|, absolute scale.

Definition
----------

The scattering intensity *I(q)* is calculated as

.. math:

    I(q) = \frac{A}{Q^n} + \frac{C}{1 + (Q\xi}^m} + B

Here the peak position is related to the d-spacing as *Q0* = 2|pi| / *d0*.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D,
where the *q* vector is defined as

.. math:

    q = \sqrt{q_x^2 + q_y^2}


.. image:: img/image175.jpg

*Figure. 1D plot using the default values (w/200 data point).*

REFERENCE
---------

None.

*2013/09/09 - Description reviewed by King, S and Parker, P.*

"""

from numpy import inf, sqrt

name = "broad_peak"
title = "Broad Lorentzian type peak on top of a power law decay"
description = """\
      I(q) = scale_p/pow(q,exponent)+scale_l/
      (1.0 + pow((fabs(q-q_peak)*length_l),exponent_l) )+ background

      List of default parameters:
      porod_scale = Porod term scaling
      porod_exp = Porod exponent
      lorentz_scale = Lorentzian term scaling
      lorentz_length = Lorentzian screening length [A]
      peak_pos = peak location [1/A]
      lorentz_exp = Lorentzian exponent
      background = Incoherent background"""
category = "shape-independent"

#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["porod_scale", "", 1.0e-05, [-inf, inf], "", "Power law scale factor"],
              ["porod_exp", "", 3.0, [-inf, inf], "", "Exponent of power law"],
              ["lorentz_scale", "", 10.0, [-inf, inf], "", "Scale factor for broad Lorentzian peak"],
              ["lorentz_length", "Ang", 50.0, [-inf, inf], "", "Lorentzian screening length"],
              ["peak_pos", "1/Ang", 0.1, [-inf, inf], "", "Peak postion in q"],
              ["lorentz_exp", "", 2.0, [-inf, inf], "", "exponent of Lorentz function"],
             ]

def Iq(q, porod_scale, porod_exp, lorentz_scale, lorentz_length, peak_pos, lorentz_exp):
    inten = (porod_scale / q ** porod_exp + lorentz_scale
             / (1.0 + (abs(q - peak_pos) * lorentz_length) ** lorentz_exp))
    return inten
Iq.vectorized = True  # Iq accepts an array of Q values

def Iqxy(qx, qy, *args):
    return Iq(sqrt(qx ** 2 + qy ** 2), *args)
Iqxy.vectorized = True # Iqxy accepts an array of Qx, Qy values


demo = dict(scale=1, background=0,
            porod_scale=1.0e-05, porod_exp=3,
            lorentz_scale=10, lorentz_length=50, peak_pos=0.1, lorentz_exp=2)

oldname = "BroadPeakModel"
oldpars = dict(porod_scale='scale_p', porod_exp='exponent_p',
               lorentz_scale='scale_l', lorentz_length='length_l', peak_pos='q_peak',
               lorentz_exp='exponent_l', scale=None)
