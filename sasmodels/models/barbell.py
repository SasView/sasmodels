#barbell model
# cylinder model
# Note: model title and parameter table are inserted automatically
r"""

Calculates the scattering from a barbell-shaped cylinder (This model simply
becomes the DumBellModel when the length of the cylinder, *L*, is set to zero).
That is, a sphereocylinder with spherical end caps that have a radius larger
than that of the cylinder and the center of the end cap radius lies outside
of the cylinder. All dimensions of the BarBell are considered to be
monodisperse. See the diagram for the details of the geometry and restrictions
on parameter values.

Definition
----------

The returned value is scaled to units of |cm^-1|\ |sr^-1|, absolute scale.

The barbell geometry is defined as

.. image:: img/barbell_geometry.jpg

where *r* is the radius of the cylinder. All other parameters are as defined
in the diagram.

Since the end cap radius
*R* >= *r* and by definition for this geometry *h* < 0, *h* is then
defined by *r* and *R* as

*h* = -1 \* sqrt(*R*\ :sup:`2` - *r*\ :sup:`2`)

The scattered intensity *I(q)* is calculated as

.. math::

    I(Q) = \frac{(\Delta \rho)^2}{V} \left< A^2(Q)\right>

where the amplitude *A(q)* is given as

.. math::

    A(Q) =&\ \pi r^2L
        {\sin\left(\tfrac12 QL\cos\theta\right)
            \over \tfrac12 QL\cos\theta}
        {2 J_1(Qr\sin\theta) \over Qr\sin\theta} \\
        &\ + 4 \pi R^3 \int_{-h/R}^1 dt
        \cos\left[ Q\cos\theta
            \left(Rt + h + {\tfrac12} L\right)\right]
        \times (1-t^2)
        {J_1\left[QR\sin\theta \left(1-t^2\right)^{1/2}\right]
             \over QR\sin\theta \left(1-t^2\right)^{1/2}}

The < > brackets denote an average of the structure over all orientations.
<*A* :sup:`2`\ *(q)*> is then the form factor, *P(q)*. The scale factor is
equivalent to the volume fraction of cylinders, each of volume, *V*. Contrast
is the difference of scattering length densities of the cylinder and the
surrounding solvent.

The volume of the barbell is

.. math::

    V = \pi r_c^2 L + 2\pi\left(\tfrac23R^3 + R^2h-\tfrac13h^3\right)


and its radius-of-gyration is

.. math::

    R_g^2 =&\ \left[ \tfrac{12}{5}R^5
        + R^4\left(6h+\tfrac32 L\right)
        + R^2\left(4h^2 + L^2 + 4Lh\right)
        + R^2\left(3Lh^2 + \tfrac32 L^2h\right) \right. \\
        &\ \left. + \tfrac25 h^5 - \tfrac12 Lh^4 - \tfrac12 L^2h^3
        + \tfrac14 L^3r^2 + \tfrac32 Lr^4 \right]
        \left( 4R^3 6R^2h - 2h^3 + 3r^2L \right)^{-1}

**The requirement that** *R* >= *r* **is not enforced in the model!** It is
up to you to restrict this during analysis.

This example dataset is produced by running the Macro PlotBarbell(),
using 200 data points, *qmin* = 0.001 |Ang^-1|, *qmax* = 0.7 |Ang^-1|,
*sld* = 4e-6 |Ang^-2| and the default model values.

.. image:: img/barbell_1d.jpg

*Figure. 1D plot using the default values (w/256 data point).*

For 2D data: The 2D scattering intensity is calculated similar to the 2D
cylinder model. For example, for |theta| = 45 deg and |phi| = 0 deg with
default values for other parameters

.. image:: img/barbell_2d.jpg

*Figure. 2D plot (w/(256X265) data points).*

.. image:: img/orientation.jpg

Figure. Definition of the angles for oriented 2D barbells.

.. image:: img/orientation2.jpg

*Figure. Examples of the angles for oriented pp against the detector plane.*

REFERENCE
---------

H Kaya, *J. Appl. Cryst.*, 37 (2004) 37 223-230

H Kaya and N R deSouza, *J. Appl. Cryst.*, 37 (2004) 508-509 (addenda and errata)

"""
from numpy import inf

name = "barbell"
title = "Cylinder with spherical end caps"
description = """
        Calculates the scattering from a barbell-shaped cylinder. That is a sphereocylinder with spherical end caps that have a radius larger than that of the cylinder and the center of the end cap
        radius lies outside of the cylinder.
        Note: As the length of cylinder(bar) -->0,it becomes a dumbbell. And when rad_bar = rad_bell, it is a spherocylinder.
        It must be that rad_bar <(=) rad_bell.
"""
category = "shape:cylinder"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["sld", "1e-6/Ang^2", 4, [-inf, inf], "", "Barbell scattering length density"],
              ["solvent_sld", "1e-6/Ang^2", 1, [-inf, inf], "", "Solvent scattering length density"],
              ["bell_radius", "Ang", 40, [0, inf], "volume", "Spherical bell radius"],
              ["radius", "Ang", 20, [0, inf], "volume", "Cylindrical bar radius"],
              ["length", "Ang", 400, [0, inf], "volume", "Cylinder bar length"],
              ["theta", "degrees", 60, [-inf, inf], "orientation", "In plane angle"],
              ["phi", "degrees", 60, [-inf, inf], "orientation", "Out of plane angle"],
             ]

source = ["lib/J1.c", "lib/gauss76.c", "barbell.c"]

# parameters for demo
demo = dict(scale=1, background=0,
            sld=6, solvent_sld=1,
            bell_radius=40, radius=20, length=400,
            theta=60, phi=60,
            radius_pd=.2, radius_pd_n=5,
            length_pd=.2, length_pd_n=5,
            theta_pd=15, theta_pd_n=0,
            phi_pd=15, phi_pd_n=0,
           )

# For testing against the old sasview models, include the converted parameter
# names and the target sasview model name.
oldname = 'BarBellModel'
oldpars = dict(sld='sld_barbell',
               solvent_sld='sld_solv', bell_radius='rad_bell',
               radius='rad_bar', length='len_bar')
