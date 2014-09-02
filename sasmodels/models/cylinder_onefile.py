# This is how a model would look if we were to stick it all into one file.
r"""
CylinderModel
=============

This model provides the form factor for a right circular cylinder with uniform
scattering length density. The form factor is normalized by the particle volume.

For information about polarised and magnetic scattering, click here_.

Definition
----------

The output of the 2D scattering intensity function for oriented cylinders is
given by (Guinier, 1955)

.. math::

    P(q,\alpha) = \frac{\text{scale}}{V}f^2(q) + \text{bkg}

where

.. math::

    f(q) = 2 (\Delta \rho) V
           \frac{\sin (q L/2 \cos \alpha)}{q L/2 \cos \alpha}
           \frac{J_1 (q r \sin \alpha)}{q r \sin \alpha}

and $\alpha$ is the angle between the axis of the cylinder and $\vec q$, $V$
is the volume of the cylinder, $L$ is the length of the cylinder, $r$ is the
radius of the cylinder, and $d\rho$ (contrast) is the scattering length density
difference between the scatterer and the solvent. $J_1$ is the first order
Bessel function.

To provide easy access to the orientation of the cylinder, we define the
axis of the cylinder using two angles $\theta$ and $\phi$. Those angles
are defined in Figure :num:`figure #CylinderModel-orientation`.

.. _CylinderModel-orientation:

.. figure:: img/image061.JPG

    Definition of the angles for oriented cylinders.

.. figure:: img/image062.JPG

    Examples of the angles for oriented pp against the detector plane.

NB: The 2nd virial coefficient of the cylinder is calculated based on the
radius and length values, and used as the effective radius for $S(Q)$
when $P(Q) \cdot S(Q)$ is applied.

The returned value is scaled to units of |cm^-1| and the parameters of
the CylinderModel are the following:

%(parameters)s

The output of the 1D scattering intensity function for randomly oriented
cylinders is then given by

.. math::

    P(q) = \frac{\text{scale}}{V}
        \int_0^{\pi/2} f^2(q,\alpha) \sin \alpha d\alpha + \text{background}

The *theta* and *phi* parameters are not used for the 1D output. Our
implementation of the scattering kernel and the 1D scattering intensity
use the c-library from NIST.

Validation of the CylinderModel
-------------------------------

Validation of our code was done by comparing the output of the 1D model
to the output of the software provided by the NIST (Kline, 2006).
Figure :num:`figure #CylinderModel-compare` shows a comparison of
the 1D output of our model and the output of the NIST software.

.. _CylinderModel-compare:

.. figure:: img/image065.JPG

    Comparison of the SasView scattering intensity for a cylinder with the
    output of the NIST SANS analysis software.
    The parameters were set to: *Scale* = 1.0, *Radius* = 20 |Ang|,
    *Length* = 400 |Ang|, *Contrast* = 3e-6 |Ang^-2|, and
    *Background* = 0.01 |cm^-1|.

In general, averaging over a distribution of orientations is done by
evaluating the following

.. math::

    P(q) = \int_0^{\pi/2} d\phi
        \int_0^\pi p(\theta, \phi) P_0(q,\alpha) \sin \theta d\theta


where $p(\theta,\phi)$ is the probability distribution for the orientation
and $P_0(q,\alpha)$ is the scattering intensity for the fully oriented
system. Since we have no other software to compare the implementation of
the intensity for fully oriented cylinders, we can compare the result of
averaging our 2D output using a uniform distribution $p(\theta, \phi) = 1.0$.
Figure :num:`figure #CylinderModel-crosscheck` shows the result of
such a cross-check.

.. _CylinderModel-crosscheck:

.. figure:: img/image066.JPG

    Comparison of the intensity for uniformly distributed cylinders
    calculated from our 2D model and the intensity from the NIST SANS
    analysis software.
    The parameters used were: *Scale* = 1.0, *Radius* = 20 |Ang|,
    *Length* = 400 |Ang|, *Contrast* = 3e-6 |Ang^-2|, and
    *Background* = 0.0 |cm^-1|.
"""

from numpy import pi, inf

name = "cylinder"
title = "Cylinder with uniform scattering length density"
description = """
     f(q)= 2*(sldCyl - sldSolv)*V*sin(qLcos(alpha/2))
            /[qLcos(alpha/2)]*J1(qRsin(alpha/2))/[qRsin(alpha)]

            P(q,alpha)= scale/V*f(q)^(2)+background
            V: Volume of the cylinder
            R: Radius of the cylinder
            L: Length of the cylinder
            J1: The bessel function
            alpha: angle betweenthe axis of the
            cylinder and the q-vector for 1D
            :the ouput is P(q)=scale/V*integral
            from pi/2 to zero of...
            f(q)^(2)*sin(alpha)*dalpha+ bkg
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
    [ "length", "Ang",  400, [0, inf], "volume",
      "Cylinder length" ],
    [ "theta", "degrees", 60, [-inf, inf], "orientation",
      "In plane angle" ],
    [ "phi", "degrees", 60, [-inf, inf], "orientation",
      "Out of plane angle" ],
    ]

source = [ "lib/J1.c", "lib/gauss76.c", "lib/cylkernel.c" ]

form_volume = """
    return M_PI*radius*radius*length;
    """

Iq = """
    const real h = REAL(0.5)*length;
    real summ = REAL(0.0);
    for (int i=0; i<76 ;i++) {
        //const real zi = ( Gauss76Z[i]*(uplim-lolim) + uplim + lolim )/2.0;
        const real zi = REAL(0.5)*(Gauss76Z[i]*M_PI_2 + M_PI_2);
        summ += Gauss76Wt[i] * CylKernel(q, radius, h, zi);
    }
    //const real form = (uplim-lolim)/2.0*summ;
    const real form = summ * M_PI_4;

    // Multiply by contrast^2, normalize by cylinder volume and convert to cm-1
    // NOTE that for this (Fournet) definition of the integral, one must MULTIPLY by Vcyl
    // The additional volume factor is for polydisperse volume normalization.
    const real s = (sld - solvent_sld) * form_volume(radius, length);
    return REAL(1.0e-4) * form * s * s;
    """

Iqxy = """
    real sn, cn; // slots to hold sincos function output

    // Compute angle alpha between q and the cylinder axis
    SINCOS(theta*M_PI_180, sn, cn);
    // # The following correction factor exists in sasview, but it can't be
    // # right, so we are leaving it out for now.
    // const real correction = fabs(cn)*M_PI_2;
    const real q = sqrt(qx*qx+qy*qy);
    const real cos_val = cn*cos(phi*M_PI_180)*(qx/q) + sn*(qy/q);
    const real alpha = acos(cos_val);

    // The following is CylKernel() / sin(alpha), but we are doing it in place
    // to avoid sin(alpha)/sin(alpha) for alpha = 0.  It is also a teensy bit
    // faster since we don't mulitply and divide sin(alpha).
    SINCOS(alpha, sn, cn);
    const real besarg = q*radius*sn;
    const real siarg = REAL(0.5)*q*length*cn;
    // lim_{x->0} J1(x)/x = 1/2,   lim_{x->0} sin(x)/x = 1
    const real bj = (besarg == REAL(0.0) ? REAL(0.5) : J1(besarg)/besarg);
    const real si = (siarg == REAL(0.0) ? REAL(1.0) : sin(siarg)/siarg);
    const real form = REAL(4.0)*bj*bj*si*si;

    // Multiply by contrast^2, normalize by cylinder volume and convert to cm-1
    // NOTE that for this (Fournet) definition of the integral, one must MULTIPLY by Vcyl
    // The additional volume factor is for polydisperse volume normalization.
    const real s = (sld - solvent_sld) * form_volume(radius, length);
    return REAL(1.0e-4) * form * s * s; // * correction;
    """

def ER(radius, length):
    ddd = 0.75*radius*(2*radius*length + (length+radius)*(length+pi*radius))
    return 0.5 * (ddd)**(1./3.)

