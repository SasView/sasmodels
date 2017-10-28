.. _orientation:

Oriented particles
==================

With two dimensional small angle diffraction data SasView will calculate scattering from
oriented particles, applicable for example to shear flow or orientation in a magnetic field.

In general we first need to define the mean, or a reference orientation of the particles with respect 
to the incoming neutron or X-ray beam. This is done using three angles: $\theta$ and $\phi$ define the 
orientation of the axis of the particle, angle $\Psi$ is defined as the orientation of the major
axis of the particle cross section with respect to its starting position along the beam direction.
The figures below are for an elliptical cross section cylinder,
but may be applied analogously to other shapes of particle.

.. note::
    It is very important to note that these angles, in particular $\theta$ and $\phi$, are NOT in general
    the same as the $\theta$ and $\phi$ appearing in equations for the scattering form factor which gives 
    the scattered intensity or indeed in the equation for scattering vector $Q$.
    The $\theta$ rotation must be applied before the $\phi$ rotation, else there is an ambiguity.

.. figure::
    orient_img/elliptical_cylinder_angle_definition.png

    Definition of angles for oriented elliptical cylinder, where axis_ratio b/a is shown >1,
    Note that rotation $\theta$, initially in the $xz$ plane, is carried out first, then
    rotation $\phi$ about the $z$ axis, finally rotation $\Psi$ is around the axis of the cylinder.
    The neutron or X-ray beam is along the $z$ axis.

.. figure::
    orient_img/elliptical_cylinder_angle_projection.png

    Some examples of the orientation angles for an elliptical cylinder, with $\Psi$ = 0.

Having established the mean direction of the particle we can then apply angular orientation distributions.
This is done by a numerical integration over a range of angles in a similar way to polydispersity for particle size.
In the current version of sasview the orientational dispersity is defined with respect to the axes of the particle.

The $\theta$ and $\phi$ parameters to orient the cylinder only appear in the model when fitting 2d data.
On introducing "Orientational Distribution" in the angles, "distribution of theta" and "distribution of phi" parameters will
appear. These are actually rotations about the axes $\delta_1$ and $\delta_2$ of the cylinder, the $b$ and $a$ axes of the
cylinder cross section. (When $\theta = \phi = 0$ these are parallel to the $Y$ and $X$ axes of the instrument.)
The third orientation distribution, in $\psi$, is about the $c$ axis of the particle. Some experimentation may be required to
understand the 2d patterns fully. A number of different shapes of distribution are available, as described for 
polydispersity, see :ref:`polydispersityhelp` .

Earlier versions of SasView had numerical integration issues in some circumstances when 
distributions passed through 90 degrees. The distributions in particle coordinates are more robust, but should still be approached 
with care for large ranges of angle.

Note that the form factors for asymmetric particles are also performing numerical integrations over one or more variables, so 
care should be taken, especially with very large particles or more extreme aspect ratios. Users can experiment with the 
values of Npts and Nsigs, the number of steps used in the integration and the range spanned in number of standard deviations.
The standard deviation is entered in units of degrees. For a rectangular (uniform) distribution the full width 
should be $\pm\sqrt(3)$ ~ 1.73 standard deviations (this may be changed soon).

Where appropriate, for best numerical results, keep $a < b < c$ and the $\theta$ distribution narrower than the $\phi$ distribution.

Some more detailed technical notes are provided in the developer section of this manual :ref:`orientation_developer` .
    
*Document History*

| 2017-10-27 Richard Heenan 