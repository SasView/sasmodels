.. fitting_sq.rst

.. Much of the following text was scraped from product.py

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

.. _Interaction_Models:

Fitting Models with Structure Factors
=====================================

**Interaction models** (previously called product models), or $P@S$ models
for short, multiply the form factor $P(Q)$ by the structure factor $S(Q)$,
modulated by the **effective radius** of the form factor. For the theory
behind this, see :ref:`PStheory` .

Parameters
^^^^^^^^^^

**Except for volfraction, when writing your own** $P$ **and** $S$ **models,**
**DO NOT give your model parameters these names!**

Many parameters are common amongst $P@S$ models, but take on specific meanings:

* *scale*:

    Overall model scale factor.

    To compute number density $n$ the volume fraction $V_f$ (parameterised as
    **volfraction**) is needed.  In most $P(Q)$ models $V_f$ is not defined and
    **scale** is used instead. Some $P(Q)$ models, such as the *vesicle*, do
    define **volfraction** and so can leave **scale** at 1.0.

    Structure factor models $S(Q)$ contain **volfraction**. In $P@S$ models
    this is *also* used as the volume fraction for the form factor model
    $P(Q)$, so these models can also leave **scale** at 1.0.  If $P(Q)$ already
    has a **volfraction** parameter, it is tied to the **volfraction** for
    $S(Q)$.

    If the volume fraction required for $S(Q)$ is *not* the volume fraction
    needed to compute the $n$ for $P(Q)$, then leave **volfraction** as the
    $V_f$ for $S(Q)$ and use **scale** to define the $V_f$ for $P(Q)$ as
    $V_f$ = **scale**  $\cdot$  **volfraction**.  This situation may occur in
    a mixed phase system where the effective volume fraction needed to compute
    the structure is much higher than the true volume fraction.

* *volfraction*:

    The volume fraction of material, $V_f$.

    For hollow shapes, **volfraction** still represents the volume fraction of
    material but the $S(Q)$ calculation needs the volume fraction *enclosed by*
    *the shape.*  To remedy this the user-specified **volfraction** is scaled
    by the ratio form:shell computed from the average form volume and average
    shell volume returned from the $P(Q)$ calculation when calculating $S(Q)$.
    The original **volfraction** is divided by the shell volume to compute the
    number density $n$ used in the $P@S$ model to get the absolute scaling on
    the final $I(Q)$.

* *radius_effective*:

    The radial distance determining the range of the $S(Q)$ interaction.

    This may be estimated from the "size" parameters $\mathbf \xi$ describing
    the form of the shape.  For example, in a system containing freely-rotating
    cylinders, the volume of space each cylinder requires to tumble will be
    much larger than the volume of the cylinder itself. Thus the *effective*
    radius of a cylinder will be larger than either its actual radius or half-
    length.

    If **radius_effective_mode = 0** (see below) it may be sensible to tie or
    constrain **radius_effective** to one or other of the "size" parameters
    describing the form of the shape (although the parameter cannot then be
    polydisperse). But **radius_effective** may also be specified directly,
    independent of the estimate from $P(Q)$.

    If **radius_effective** is calculated by $P(Q)$, it will be the
    weighted average of the effective radii computed for the polydisperse
    shape parameters, and that average is used to compute $S(Q)$. When
    specified directly, the value of **radius_effective** may be
    polydisperse, and $S(Q)$ will be averaged over a range of effective
    radii. Whether this makes any physical sense will depend on the system.

.. note::

   The following additional parameters are only available in SasView 4.3 and
   later.

 *radius_effective_mode*:

    Defines how the effective radius (parameter **radius_effective**) should
    be computed from the parameters of the shape.

    When **radius_effective_mode = 0** then the unconstrained
    **radius_effective** parameter in the $S(Q)$ model is used. *This is the*
    *default in SasView versions 4.2 and earlier*. Otherwise, in SasView 4.3
    and later, **radius_effective_mode = k** represents an index in a list of
    alternative **radius_effective** calculations.

    In SasView 4.3 and later **k** must be entered as an integer (and it will
    be necessary to read the source code file to discover what calculations the
    modes represent), but in SasView 5.0 and later the options appear in a
    drop-down box.

    For example, the *ellipsoid* model defines the following
    **radius_effective_modes**::

        1 => average curvature
        2 => equivalent volume sphere
        3 => min radius
        4 => max radius

    Note: **radius_effective_mode** will only appear in the parameter table if
    the model defines the list of modes, otherwise it will be set permanently
    to 0 for the user-defined effective radius.

    **WARNING! If** $P(Q)$ **is multiplied by** $S(Q)$ **in the FitPage,**
    **instead of being generated in the Sum|Multi dialog, the**
    **radius_effective used is constrained (equivalent to**
    **radius_effective_mode = 1)**.

* *structure_factor_mode*:

    The type of structure factor calculation to use.

    If the $P@S$ model supports the $\beta(Q)$ *decoupling correction*
    [#kotlarchyk]_ then **structure_factor_mode** will appear in the
    parameter table after the $S(Q)$ parameters.

    If **structure_factor_mode = 0** then the
    *local monodisperse approximation* will be used, i.e.:

    .. math::
        I(Q) = \text{scale} \frac{V_f}{V} P(Q) S(Q) + \text{background}

    where $P(Q) = \langle F(Q)^2 \rangle$. *This is the default in SasView*
    *versions 4.x and earlier*.

    If **structure_factor_mode = 1** then the $\beta(Q)$ correction will be
    used, i.e.:

    .. math::
        I(Q) = \text{scale} \frac{V_f}{V} P(Q) [ 1 + \beta(Q) (S(Q) - 1) ]
        + \text{background}

    The $\beta(Q)$ decoupling approximation has the effect of damping the
    oscillations in the normal (local monodisperse) $S(Q)$. When $\beta(Q) = 1$
    the local monodisperse approximation is recovered. *This mode is only*
    *available in SasView 4.3 and later*.

    More mode options may appear in future as more complicated operations are
    added.

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

Related sections
^^^^^^^^^^^^^^^^

See also:

:ref:`PStheory`

:ref:`polydispersityhelp`

:ref:`Resolution_Smearing`

:ref:`orientation`

References
^^^^^^^^^^

.. [#kotlarchyk] Kotlarchyk, M.; Chen, S.-H. *J. Chem. Phys.*, 1983, 79, 2461

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

*Document History*

| 2019-03-31 Paul Kienzle, Steve King & Richard Heenan
| 2021-11-03 Steve King
| 2022-10-30 Steve King
