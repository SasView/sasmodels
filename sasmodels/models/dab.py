r"""

Scattering model class for the DAB (Debye Anderson Brumberger) Model

Definition
----------

Calculates the scattering from a randomly distributed, two-phase system based on
the Debye-Anderson-Brumberger (DAB) model for such systems. The two-phase system
is characterized by a single length scale, the correlation length, which is a
measure of the average spacing between regions of phase 1 and phase 2. **The
model also assumes smooth interfaces between the phases** and hence exhibits
Porod behavior (I ~ *q*\ :sup:`-4`) at large *q* (*QL* >> 1).

The DAB model is ostensibly a development of the earlier Debye-Bueche model.

*Definition*

.. math:: I(q) = \text{scale}\cdot\frac{L^3}{(1 + (q \cdot L)^2)^2} + \text{background}

where scale is

.. math:: \text{scale} = 8 \pi \phi (1-\phi)(\Delta \rho)^2

and the parameter *L* is the correlation length.

For 2D data: The 2D scattering intensity is calculated in the same way as 1D,
where the *q* vector is defined as

.. math:: q = \sqrt{q_x^2 + q_y^2}

==============  ========  =============
Parameter name  Units     Default value
==============  ========  =============
scale           None      1.0
corr length L   |Ang|     50.0
background      |cm^-1|   0.0
==============  ========  =============


.. figure:: img/dab_1d.jpg

   1D plot using the default values (w/200 data point).


Reference
---------

P Debye, H R Anderson, H Brumberger, *Scattering by an Inhomogeneous Solid. II.
The Correlation Function and its Application*, *J. Appl. Phys.*, 28(6) (1957) 679

P Debye, A M Bueche, *Scattering by an Inhomogeneous Solid*, *J. Appl. Phys.*,
20 (1949) 518

*2013/09/09 - Description reviewed by King, S and Parker, P.*

"""

from numpy import inf

name = "dab"
title = "DAB (Debye Anderson Brumberger) Model"
description = """\

F(q)= scale * L^3/(1 + (q*L)^2)^2

L: the correlation length

"""
category = "shape-independent"

#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["length", "Ang", 50.0, [0, inf], "", "correlation length"],
             ]

Iq = """
    double numerator   = pow(length, 3);
    double denominator = pow(1 + pow(q*length,2), 2);
    
    return numerator / denominator ;
    """

Iqxy = """
    // never called since no orientation or magnetic parameters.
    //return -1.0;
    return Iq(sqrt(qx*qx + qy*qy), length);
    """

# ER defaults to 1.0

# VR defaults to 1.0

demo = dict(scale=1, background=0, length=50)
oldname = "DABModel"
oldpars = dict(length='length')
