.. _magnetism:

Polarisation/Magnetic Scattering
=======================================================

In earlier versions of SasView magnetic scattering was implemented in just five 
(2D) models

*  :ref:`sphere`
*  :ref:`core-shell-sphere`
*  :ref:`core-multi-shell`
*  :ref:`cylinder`
*  :ref:`parallelepiped`

From SasView 4.x it is implemented on most models in the 'shape' category.

In general, the scattering length density (SLD = $\beta$) in each region where the
SLD is uniform, is a combination of the nuclear and magnetic SLDs and, for polarised
neutrons, also depends on the spin states of the neutrons.

For magnetic scattering, only the magnetization component $\mathbf{M_\perp}$
perpendicular to the scattering vector $\mathbf{Q}$ contributes to the magnetic
scattering length.


.. figure::
    mag_img/mag_vector.bmp

The magnetic scattering length density is then

.. math::
    \beta_M = \dfrac{\gamma r_0}{2\mu_B}\sigma \cdot
    \mathbf{M_\perp} = D_M\sigma \cdot \mathbf{M_\perp}

where $\gamma = -1.913$ is the gyromagnetic ratio, $\mu_B$ is the
Bohr magneton, $r_0$ is the classical radius of electron, and $\sigma$
is the Pauli spin.

Assuming that incident neutrons are polarized parallel (+) and anti-parallel (-)
to the $x'$ axis, the possible spin states after the sample are then

No spin-flips (+ +) and (- -)

Spin-flips    (+ -) and (- +)

.. figure::
    mag_img/M_angles_pic.bmp

If the angles of the $Q$ vector and the spin-axis $x'$ to the $x$ - axis are
$\phi$ and $\theta_{up}$, respectively, then, depending on the spin state of the
neutrons, the scattering length densities, including the nuclear scattering
length density ($\beta{_N}$) are

.. math::
    \beta_{\pm\pm} =  \beta_N \mp D_M M_{\perp x'}
    \text{ when there are no spin-flips}

and

.. math::
    \beta_{\pm\mp} =  -D_M (M_{\perp y'} \pm iM_{\perp z'})
    \text{ when there are}

where

.. math::
    M_{\perp x'} = M_{0q_x}\cos(\theta_{up})+M_{0q_y}\sin(\theta_{up}) \\
    M_{\perp y'} = M_{0q_y}\cos(\theta_{up})-M_{0q_x}\sin(\theta_{up}) \\
    M_{\perp z'} = M_{0z} \\
    M_{0q_x} = (M_{0x}\cos\phi - M_{0y}\sin\phi)\cos\phi \\
    M_{0q_y} = (M_{0y}\sin\phi - M_{0x}\cos\phi)\sin\phi

Here, $M_{0x}$, $M_{0x}$, $M_{0z}$ are the x, y and z components
of the magnetization vector given in the laboratory xyz frame given by

.. math::
    M_{0x} = M_0\cos\theta_M\cos\phi_M \\
    M_{0y} = M_0\sin\theta_M \\
    M_{0z} = -M_0\cos\theta_M\sin\phi_M

and the magnetization angles $\theta_M$ and $\phi_M$ are defined in
the figure above.

The user input parameters are:

===========   ================================================================
 M0_sld        = $D_M M_0$
 Up_theta      = $\theta_{up}$
 M_theta       = $\theta_M$
 M_phi         = $\phi_M$
 Up_frac_i     = (spin up)/(spin up + spin down) neutrons *before* the sample
 Up_frac_f     = (spin up)/(spin up + spin down) neutrons *after* the sample
===========   ================================================================

.. note::
    The values of the 'Up_frac_i' and 'Up_frac_f' must be in the range 0 to 1.

.. note::
    This help document was last changed by Steve King, 02May2015
