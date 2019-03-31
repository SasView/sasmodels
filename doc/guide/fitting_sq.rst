.. fitting_sq.rst

.. Much of the following text was scraped from product.py

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

.. _Product_Models:

Fitting Models with Structure Factors
-------------------------------------

.. note::

   This help document is under development

.. figure:: p_and_s_buttons.png

**Product models**, or $P@S$ models for short, multiply the form factor
$P(Q)$ by the structure factor $S(Q)$, modulated by the **effective radius**
of the form factor.


Scattering at vector $\mathbf Q$ for an individual particle with
shape parameters $\mathbf\xi$ and contrast $\rho_c(\mathbf r, \mathbf\xi)$
is computed from the square of the amplitude, $F(\mathbf Q, \mathbf\xi)$, as

.. math::
    I(\mathbf Q) = F(\mathbf Q, \mathbf\xi) F^*(\mathbf Q, \mathbf\xi)
        \big/ V(\mathbf\xi)

with particle volume $V(\mathbf \xi)$ and

.. math::
    F(\mathbf Q, \mathbf\xi) = \int_{\mathbb R^3} \rho_c(\mathbf r, \mathbf\xi)
        e^{i \mathbf Q \cdot \mathbf r} \,\mathrm d \mathbf r

The 1-D scattering pattern for monodisperse particles uses the orientation
average in spherical coordinates,

.. math::
    I(Q) = n \langle F F^*\rangle = \frac{n}{4\pi}
    \int_{\theta=0}^{\pi} \int_{\phi=0}^{2\pi}
    F F^* \sin(\theta) \,\mathrm d\phi \mathrm d\theta

where $F(\mathbf Q,\mathbf\xi)$ uses
$\mathbf Q = [Q \sin\theta\cos\phi, Q \sin\theta\sin\phi, Q \cos\theta]^T$.
A $u$-substitution may be used, with $\alpha = \cos \theta$,
$\surd(1 - \alpha^2) = \sin \theta$, and
$\mathrm d\alpha = -\sin\theta\,\mathrm d\theta$.
Here,

.. math:: n = V_f/V(\mathbf\xi)

is the number density of scatterers estimated from the volume fraction
of particles in solution. In this formalism, each incoming
wave interacts with exactly one particle before being scattered into the
detector. All interference effects are within the particle itself.
The detector accumulates counts in proportion to the relative probability
at each pixel. The extension to heterogeneous systems is simply a matter of
adding the scattering patterns in proportion to the number density of each
particle. That is, given shape parameters $\mathbf\xi$ with probability
$P_\mathbf{\xi}$,

.. math::

    I(Q) = \int_\Xi n(\mathbf\xi) \langle F F^* \rangle \,\mathrm d\xi
         = V_f\frac{\int_\Xi P_\mathbf{\xi} \langle F F^* \rangle
         \,\mathrm d\mathbf\xi}{\int_\Xi P_\mathbf\xi V(\mathbf\xi)\,\mathrm d\mathbf\xi}

This approximation is valid in the dilute limit, where particles are
sufficiently far apart that the interaction between them can be ignored.

As concentration increases, a structure factor term $S(Q)$ can be included,
giving the monodisperse approximation for the interaction between particles,
with

.. math:: I(Q) = n \langle F F^* \rangle S(Q)

For particles without spherical symmetry, the decoupling approximation (DA)
is more accurate, with

.. math::

    I(Q) = n [\langle F F^* \rangle
        + \langle F \rangle \langle F \rangle^* (S(Q) - 1)]

Or equivalently,

.. math:: I(Q) = P(Q)[1 + \beta\,(S(Q) - 1)]

with form factor $P(Q) = n \langle F F^* \rangle$ and
$\beta = \langle F \rangle \langle F \rangle^* \big/ \langle F F^* \rangle$.
These approximations can be extended to heterogeneous systems using averages
over size, $\langle \cdot \rangle_\mathbf\xi = \int_\Xi P_\mathbf\xi \langle\cdot\rangle\,\mathrm d\mathbf\xi \big/ \int_\Xi P_\mathbf\xi \,\mathrm d\mathbf\xi$ and setting
$n = V_f\big/\langle V \rangle_\mathbf\xi$.
Further improvements can be made using the local monodisperse
approximation (LMA) or using partial structure factors, as described
in \cite{bresler_sasfit:_2015}.

Many parameters are common amongst $P@S$ models, and take on specific meanings:

* *scale*:

    Overall model scale factor.

    To compute number density $n$ the volume fraction $V_f$ is needed.  In
    most $P(Q)$ models $V_f$ is not defined and **scale** is used instead.
    Some $P(Q)$ models, such as *vesicle*, do define **volfraction** and so
    can leave **scale** at 1.0.

    The structure factor model $S(Q)$ has **volfraction**.  This is also used
    as the volume fraction for the form factor model $P(Q)$, replacing the
    **volfraction** parameter if it exists in $P$. This means that
    $P@S$ models can leave **scale** at 1.0.

    If the volume fraction required for $S(Q)$ is *not* the volume fraction
    needed to compute the number density for $P(Q)$, then leave
    **volfraction** as the volume fraction for $S(Q)$ and use
    **scale** to define the volume fraction for $P(Q)$ as
    $V_f$ = **scale**  $\cdot$  **volfraction**.  This situation may
    occur in a mixed phase system where the effective volume
    fraction needed to compute the structure is much higher than the
    true volume fraction.

* *volfraction*:

    The volume fraction of material.

    For hollow shapes, **volfraction** still represents the volume fraction of
    material but the $S(Q)$ calculation needs the volume fraction *enclosed by*
    *the shape.*  Thus the user-specified **volfraction** is scaled by the ratio
    form:shell computed from the average form volume and average shell volume
    returned from the $P(Q)$ calculation when calculating $S(Q)$.  The original
    **volfraction** is divided by the shell volume to compute the number
    density $n$ used in $P@S$ to get the absolute scaling on the final $I(Q)$.

* *radius_effective*:

    The radial distance determining the range of the $S(Q)$ interaction.

    This may be estimated from the "size" parameters $\mathbf \xi$ describing
    the form of the shape.  For example, in a system containing freely-rotating
    cylinders, the volume of space each cylinder requires to tumble will be
    much larger than the volume of the cylinder itself.  Thus the effective
    radius will be larger than either the radius or the half-length of the
    cylinder.  It may be sensible to tie or constrain **radius_effective**
    to one or other of these "size" parameters. **radius_effective** may
    also be specified directly, independent of the estimate from $P(Q)$.

    If it is calculated by $P(Q)$, **radius_effective** will be the
    weighted average of the effective radii computed for the polydisperse
    shape parameters, and that average used to compute $S(Q)$.  When
    specified directly, the value of **radius_effective** may be
    polydisperse, and $S(Q)$ will be averaged over a range of effective
    radii.  Whether this makes any physical sense will depend on the system.

* *radius_effective_mode*:

    Selects the **radius_effective** value to use.

    When **radius_effective_mode = 0** then the **radius_effective**
    parameter in the $P@S$ model is used.  Otherwise
    **radius_effective_mode = k** is the index into the list of
    **radius_effective_modes** defined by the model indicating how the
    effective radius should be computed from the parameters of the shape.
    For example, the *ellipsoid* model defines the following::

        1 => average curvature
        2 => equivalent volume sphere
        3 => min radius
        4 => max radius

    **radius_effective_mode** will only appear in the parameter table if
    the model defines the list of modes, otherwise it will be set permanently
    to 0 for user defined effective radius.

* *structure_factor_mode*:

    The type of structure factor calculation to use.

    If the $P@S$ model supports the $\beta(Q)$ *decoupling correction* [1]
    then **structure_factor_mode** will appear in the parameter table after
    the $S(Q)$ parameters.

    If **structure_factor_mode = 0** then the
    *local monodisperse approximation* will be used, i.e.:

    .. math::
        I(Q) = \text{scale} \frac{V_f}{V} P(Q) S(Q) + \text{background}

    where $P(Q) = \langle F(Q)^2 \rangle$.

    If **structure_factor_mode = 1** then the $\beta(Q)$ correction will be
    used, i.e.:

    .. math::
        I(Q) = \text{scale} \frac{V_f}{V} P(Q) [ 1 + \beta(Q) (S(Q) - 1) ]
        + \text{background}

    The $\beta(Q)$ decoupling approximation has the effect of damping the
    oscillations in the normal (local monodisperse) $S(Q)$. When $\beta(Q) = 1$
    the local monodisperse approximation is recovered.

    More mode options may appear in future as more complicated operations are
    added.

References
^^^^^^^^^^

.. [#] Kotlarchyk, M.; Chen, S.-H. *J. Chem. Phys.*, 1983, 79, 2461

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

*Document History*

| 2019-03-30 Paul Kienzle & Steve King
