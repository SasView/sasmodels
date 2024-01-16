.. currentmodule:: sasmodels
.. theory.rst

.. Much of the following text was scraped from fitting_sq.py

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

.. _PStheory:

Theory
======

Scattering at vector $\mathbf Q$ for an individual particle with
shape parameters $\mathbf\xi$ and contrast $\rho_c(\mathbf r, \mathbf\xi)$
is computed from the square of the amplitude, $F(\mathbf Q, \mathbf\xi)$, as

.. math::
    I(\mathbf Q) = F(\mathbf Q, \mathbf\xi) F^*(\mathbf Q, \mathbf\xi)
        \big/ V(\mathbf\xi)

with the particle volume $V(\mathbf \xi)$ and

.. math::
    F(\mathbf Q, \mathbf\xi) = \int_{\mathbb R^3} \rho_c(\mathbf r, \mathbf\xi)
        e^{i \mathbf Q \cdot \mathbf r} \,\mathrm d \mathbf r = F

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

is the number density of scatterers estimated from the volume fraction $V_f$
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

For particles without spherical symmetry, the decoupling approximation
is more accurate, with

.. math::

    I(Q) = n [\langle F F^* \rangle
        + \langle F \rangle \langle F \rangle^* (S(Q) - 1)]

Or equivalently,

.. math:: I(Q) = P(Q)[1 + \beta\,(S(Q) - 1)]

with the form factor $P(Q) = n \langle F F^* \rangle$ and
$\beta = \langle F \rangle \langle F \rangle^* \big/ \langle F F^* \rangle$.
These approximations can be extended to heterogeneous systems using averages
over size, $\langle \cdot \rangle_\mathbf\xi = \int_\Xi P_\mathbf\xi \langle\cdot\rangle\,\mathrm d\mathbf\xi \big/ \int_\Xi P_\mathbf\xi \,\mathrm d\mathbf\xi$ and setting
$n = V_f\big/\langle V \rangle_\mathbf\xi$.

Further improvements can be made using the local monodisperse
approximation (LMA) or using partial structure factors as done in [#bressler]_,
but these are not implemented in this code.

For hollow shapes, *volfraction* is computed from the material in the
shell rather than the shell plus solvent inside the shell.  Using
$V_e(\mathbf\xi)$ as the enclosed volume of the shell plus solvent and
$V_c(\mathbf\xi)$ as the core volume of solvent inside the shell, we
can compute the average enclosed and shell volumes as

.. math::
    :nowrap:

    \begin{align*}
    \langle V_e \rangle &= \frac{
        \int_\Xi P_\mathbf\xi V_e(\mathbf\xi)\,\mathrm d\mathbf\xi
    }{ \int_\Xi P_\mathbf\xi\,\mathrm d\mathbf \xi } \\
    \langle V_s \rangle &= \frac{
        \int_\Xi P_\mathbf\xi (V_e(\mathbf\xi) - V_c(\mathbf\xi))\,\mathrm d\mathbf\xi
    }{ \int_\Xi P_\mathbf\xi\,\mathrm d\mathbf \xi }
    \end{align*}

Given $n$ particles and a total solvent volume $V_\text{out}$ outside the
shells, the volume fraction of the shell, $\phi_s$ and the shell plus
enclosed solvent $\phi_e$ are

.. math::
    :nowrap:

    \begin{align*}
    \phi_s &= \frac{n \langle V_s \rangle}{n \langle V_s \rangle + n \langle V_c \rangle + V_\text{out}}
           = \frac{n \langle V_s \rangle}{V_\text{total}} \\
    \phi_e &= \frac{n \langle V_e \rangle}{n \langle V_e \rangle + V_\text{out}}
           = \frac{n \langle V_e \rangle}{V_\text{total}}
    \end{align*}

Dividing gives

.. math::

    \frac{\phi_S}{\phi_P} = \frac{\langle V_e \rangle}{\langle V_s \rangle}

so the enclosed volume fraction can be computed from the shell volume fraction
and the form:shell volume ratio as

.. math::

    \phi_S = \phi_P \langle V_e \rangle \big/ \langle V_s \rangle

.. note::

    Prior to Sasmodels v1.0.5 (Nov 2020), the intermediate $P(Q)$ returned by
    the interaction calculator did not incorporate the volume normalization and
    so $I(Q) \ne P(Q) S(Q)$. This became apparent when $P(Q)$ and $I(Q)$ were
    plotted together. Further details can be found `here <https://github.com/SasView/sasview/issues/1698>`_.

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

Related sections
^^^^^^^^^^^^^^^^

The concepts described above are developed further in the following sections:

:ref:`polydispersityhelp`

:ref:`Resolution_Smearing`

:ref:`Interaction_Models`

:ref:`orientation`

References
^^^^^^^^^^

.. [#bressler] Bressler I., Kohlbrecher J., Thunemann A.F.
   *J. Appl. Crystallogr.* 48 (2015) 1587-1598

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

*Document History*

| 2019-03-31 Paul Kienzle, Steve King & Richard Heenan
| 2021-11-03 Steve King
| 2022-10-29 Steve King
