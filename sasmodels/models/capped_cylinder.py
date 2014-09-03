r"""
Calculates the scattering from a cylinder with spherical section end-caps.
This model simply becomes the a convex lens when the length of the cylinder
$L=0$, that is, a sphereocylinder with end caps that have a radius larger
than that of the cylinder and the center of the end cap radius lies within
the cylinder. See the diagram for the details of the geometry and
restrictions on parameter values.

Definitions
-----------

The returned value is scaled to units of |cm^-1|\ |sr^-1|, absolute scale.

The capped cylinder geometry is defined as

.. image:: img/capped_cylinder_geometry.jpg

where $r$ is the radius of the cylinder. All other parameters are as defined
in the diagram. Since the end cap radius $R \ge r$ and by definition for this
geometry $h < 0$, $h$ is then defined by $r$ and $R$ as

.. math::

    h = - \sqrt{R^2 - r^2}

The scattered intensity $I(Q)$ is calculated as

.. math::

    I(Q) = \frac{(\Delta \rho)^2}{V} \left< A^2(Q)\right>

where the amplitude $A(Q)$ is given as

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

The $\left< \ldots \right>$ brackets denote an average of the structure over
all orientations. $\left< A^2(Q)\right>$ is then the form factor, $P(Q)$.
The scale factor is equivalent to the volume fraction of cylinders, each of
volume, $V$. Contrast is the difference of scattering length densities of
the cylinder and the surrounding solvent.

The volume of the capped cylinder is (with $h$ as a positive value here)

.. math::

    V = \pi r_c^2 L + \tfrac{2\pi}{3}(R-h)^2(2R + h)


and its radius-of-gyration is

.. math::

    R_g^2 =&\ \left[ \tfrac{12}{5}R^5
        + R^4\left(6h+\tfrac32 L\right)
        + R^2\left(4h^2 + L^2 + 4Lh\right)
        + R^2\left(3Lh^2 + \tfrac32 L^2h\right) \right. \\
        &\ \left. + \tfrac25 h^5 - \tfrac12 Lh^4 - \tfrac12 L^2h^3
        + \tfrac14 L^3r^2 + \tfrac32 Lr^4 \right]
        \left( 4R^3 6R^2h - 2h^3 + 3r^2L \right)^{-1}


.. note::

    The requirement that $R \ge r$ is not enforced in the model!
    It is up to you to restrict this during analysis.

:num:`Figure #capped-cylinder-1d` shows the output produced by
a running the 1D capped cylinder model, using *qmin* = 0.001 |Ang^-1|,
*qmax* = 0.7 |Ang^-1| and  the default values of the parameters.

.. _capped-cylinder-1d:

.. figure:: img/capped_cylinder_1d.jpg

    1D plot using the default values (w/256 data point).

The 2D scattering intensity is calculated similar to the 2D cylinder model.
:num:`Figure #capped-cylinder-2d` shows the output for $\theta=45^\circ$
and $\phi=0^\circ$ with default values for the other parameters.

.. _capped-cylinder-2d:

.. figure:: img/capped_cylinder_2d.jpg

    2D plot (w/(256X265) data points).

.. figure:: img/orientation.jpg

    Definition of the angles for oriented 2D cylinders.

.. figure:: img/orientation2.jpg

    Examples of the angles for oriented pp against the detector plane.

REFERENCE

H Kaya, *J. Appl. Cryst.*, 37 (2004) 223-230

H Kaya and N-R deSouza, *J. Appl. Cryst.*, 37 (2004) 508-509 (addenda and errata)
"""
from numpy import inf

name = "capped_cylinder"
title = "Right circular cylinder with spherical end caps and uniform SLD"
description = """That is, a sphereocylinder
	with end caps that have a radius larger than
	that of the cylinder and the center of the
	end cap radius lies within the cylinder.
	Note: As the length of cylinder -->0,
	it becomes a ConvexLens.
	It must be that radius <(=) cap_radius.
	[Parameters];
	scale: volume fraction of spheres,
	background:incoherent background,
	radius: radius of the cylinder,
	length: length of the cylinder,
	cap_radius: radius of the semi-spherical cap,
	sld: SLD of the capped cylinder,
	solvent_sld: SLD of the solvent.
"""

parameters = [
    #   [ "name", "units", default, [lower, upper], "type",
    #     "description" ],
    [ "sld", "1e-6/Ang^2", 4, [-inf,inf], "",
      "Cylinder scattering length density" ],
    [ "solvent_sld", "1e-6/Ang^2", 1, [-inf,inf], "",
      "Solvent scattering length density" ],
    [ "radius", "Ang",  20, [0, inf], "volume",
      "Cylinder radius" ],
    # TODO: use an expression for cap radius with fixed bounds.
    # The current form requires cap radius R bigger than cylinder radius r.
    # Could instead use R/r in [1,inf], r/R in [0,1], or the angle between
    # cylinder and cap in [0,90].  The problem is similar for the barbell
    # model.  Propose r/R in [0,1] in both cases, with the model specifying
    # cylinder radius in the capped cylinder model and sphere radius in the
    # barbell model.  This leads to the natural value of zero for no cap
    # in the capped cylinder, and zero for no bar in the barbell model.  In
    # both models, one would be a pill.
    [ "cap_radius", "Ang", 20, [0, inf], "volume",
      "Cap radius" ],
    [ "length", "Ang",  400, [0, inf], "volume",
      "Cylinder length" ],
    [ "theta", "degrees", 60, [-inf, inf], "orientation",
      "In plane angle" ],
    [ "phi", "degrees", 60, [-inf, inf], "orientation",
      "Out of plane angle" ],
    ]

source = [ "lib/J1.c", "lib/gauss76.c", "capped_cylinder.c" ]



