#poly_gauss_coil model
#conversion of Poly_GaussCoil.py
#converted by Steve King, Mar 2016


 
r"""
This empirical model describes the scattering from *polydisperse* polymer chains in theta solvents or polymer melts, assuming a Schulz-Zimm type molecular weight distribution.

To describe the scattering from *monodisperse* polymer chains, see the :ref:`mono_gauss_coil <mono-gauss-coil>` model.

Definition
----------

     *I(q)* = *scale* |cdot| *I* \ :sub:`0` |cdot| *P(q)* + *background*
	 
where

     *I*\ :sub:`0` = |phi|\ :sub:`poly` |cdot| *V* |cdot| (|rho|\ :sub:`poly` - |rho|\ :sub:`solv`)\ :sup:`2`

     *P(q)* = 2 [(1 + UZ)\ :sup:`-1/U` + Z - 1] / [(1 + U) Z\ :sup:`2`]
	 
	 *Z* = [(*q R*\ :sub:`g`)\ :sup:`2`] / (1 + 2U)
	 
	 *U* = (Mw / Mn) - 1 = (*polydispersity ratio*) - 1

and

	 *V* = *M* / (*N*\ :sub:`A` |delta|)
	 
Here, |phi|\ :sub:`poly`, is the volume fraction of polymer, *V* is the volume of a polymer coil, *M* is the molecular weight of the polymer, *N*\ :sub:`A` is Avogadro's Number, |delta| is the bulk density of the polymer, |rho|\ :sub:`poly` is the sld of the polymer, |rho|\ :sub:`solv` is the sld of the solvent, and *R*\ :sub:`g` is the radius of gyration of the polymer coil.

The 2D scattering intensity is calculated in the same way as the 1D, but where the *q* vector is redefined as

.. image:: img/2d_q_vector.gif

References
----------

O Glatter and O Kratky (editors), *Small Angle X-ray Scattering*, Academic Press, (1982)
Page 404.

J S Higgins, H C Benoit, *Polymers and Neutron Scattering*, Oxford Science Publications, (1996).

S M King, *Small Angle Neutron Scattering* in *Modern Techniques for Polymer Characterisation*, Wiley, (1999).

http://www.ncnr.nist.gov/staff/hammouda/distance_learning/chapter_28.pdf
"""

from numpy import inf, sqrt, power

name =  "poly_gauss_coil"
title =  "Scattering from polydisperse polymer coils"

description =  """
    Evaluates the scattering from 
	polydisperse polymer chains.
    """
category =  "shape-independent"

#             ["name", "units", default, [lower, upper], "type", "description"],
parameters =  [["i_zero", "1/cm", 1.0, [-inf, inf], "", "Intensity at q=0"],
               ["radius_gyration", "Ang", 50.0, [0.0, inf], "", "Radius of gyration"],
               ["polydispersity", "None", 2.0, [1.0, inf], "", "Polymer Mw/Mn"]]

# NB: Scale and Background are implicit parameters on every model
def Iq(q, radius_gyration, polydispersity):
    # pylint: disable = missing-docstring
	u = polydispersity - 1.0
    # TO DO
	# should trap the case of polydispersity = 1 by switching to a taylor expansion
	minusoneonu = -1.0 / u
	z = ((x * radius_gyration) * (x * radius_gyration)) / (1.0 + 2.0 * u)
	if x == 0:
	   inten = i_zero * 1.0
	else:
	   inten = i_zero * 2.0 * (power((1.0 + u * z),minusoneonu) + z - 1.0 ) / ((1.0 + u) * (z * z))
	return inten
Iq.vectorized =  True  # Iq accepts an array of q values

def Iqxy(qx, qy, *args):
    # pylint: disable = missing-docstring
    return Iq(sqrt(qx ** 2 + qy ** 2), *args)
Iqxy.vectorized =  True # Iqxy accepts an array of qx, qy values

demo =  dict(scale = 1.0,
            i_zero = 1.0,
            radius_gyration = 50.0,
            polydispersity = 2.0,
            background = 0.0)

oldname =  "Poly_GaussCoil"
oldpars =  dict(scale = 'scale',
               radius_gyration = 'rg',
               polydispersity = 'poly_m',
               background = 'background')

tests =  [
    [{'scale': 1.0, 'radius_gyration': 50.0, 'polydispersity': 2.0, 'background': 0.0},
     [0.0106939, 0.469418], [0.912993, 0.0054163]],
    ]
