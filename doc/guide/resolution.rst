.. sm_help.rst

.. This is a port of the original SasView html help file to ReSTructured text
.. by S King, ISIS, during SasView CodeCamp-III in Feb 2015.


.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

Resolution Functions
====================

Sometimes the instrumental geometry used to acquire the experimental data has
an impact on the clarity of features in the reduced scattering curve. For
example, peaks or fringes might be slightly broadened. This is known as
*Q resolution smearing*. To compensate for this effect one can either try and
remove the resolution contribution - a process called *desmearing* - or add the
resolution contribution into a model calculation/simulation (which by definition
will be exact) to make it more representative of what has been measured
experimentally - a process called *smearing*. Sasmodels does the latter.

Both smearing and desmearing rely on functions to describe the resolution
effect. Sasmodels provides three smearing algorithms:

*  *Slit Smearing*
*  *Pinhole Smearing*
*  *2D Smearing*

The $Q$ resolution values should be determined by the data reduction software
for the instrument and stored with the data file.  If not, they will need to
be set manually before fitting.


.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

Slit Smearing
-------------

**This type of smearing is normally only encountered with data from X-ray Kratky**
**cameras or X-ray/neutron Bonse-Hart USAXS/USANS instruments.**

The slit-smeared scattering intensity is defined by

.. math::
    I_s = \frac{1}{\text{Norm}}
          \int_{-\infty}^{\infty} dv\, W_v(v)
          \int_{-\infty}^{\infty} du\, W_u(u)\,
          I\left(\sqrt{(q+v)^2 + |u|^2}\right)

where *Norm* is given by

.. math:: \int_{-\infty}^{\infty} dv\, W_v(v) \int_{-\infty}^{\infty} du\, W_u(u)

**[Equation 1]**

The functions $W_v(v)$ and $W_u(u)$ refer to the slit width weighting
function and the slit height weighting determined at the given $q$ point,
respectively. It is assumed that the weighting function is described by a
rectangular function, such that

.. math:: W_v(v) = \delta(|v| \leq \Delta q_v)

**[Equation 2]**

and

.. math:: W_u(u) = \delta(|u| \leq \Delta q_u)

**[Equation 3]**

so that $\Delta q_\alpha = \int_0^\infty d\alpha\, W_\alpha(\alpha)$
for $\alpha$ as $v$ and $u$.

Here $\Delta q_u$ and $\Delta q_v$ stand for the the slit height (FWHM/2)
and the slit width (FWHM/2) in $q$ space.

This simplifies the integral in Equation 1 to

.. math::

    I_s(q) = \frac{2}{\text{Norm}}
             \int_{-\Delta q_v}^{\Delta q_v} dv
             \int_{0}^{\Delta q_u}
             du\, I\left(\sqrt{(q+v)^2 + u^2}\right)

**[Equation 4]**

which may be solved numerically, depending on the nature of
$\Delta q_u$ and $\Delta q_v$.

Solution 1
^^^^^^^^^^

**For** $\Delta q_v = 0$ **and** $\Delta q_u = \text{constant}$

.. math::

    I_s(q) \approx \int_0^{\Delta q_u} du\, I\left(\sqrt{q^2+u^2}\right)
           = \int_0^{\Delta q_u} d\left(\sqrt{q'^2-q^2}\right)\, I(q')

For discrete $q$ values, at the $q$ values of the data points and at the $q$
values extended up to $q_N = q_i + \Delta q_u$ the smeared
intensity can be approximately calculated as

.. math::

    I_s(q_i)
    \approx \sum_{j=i}^{N-1} \left[\sqrt{q_{j+1}^2 - q_i^2} - \sqrt{q_j^2 - q_i^2}\right]\, I(q_j)
            \sum_{j=1}^{N-1} W_{ij}\, I(q_j)

**[Equation 5]**

where $W_{ij} = 0$ for $I_s$ when $j < i$ or $j > N-1$.

Solution 2
^^^^^^^^^^

**For** $\Delta q_v = \text{constant}$ **and** $\Delta q_u = 0$

Similar to Case 1

.. math::

    I_s(q_i)
    \approx \sum_{j=p}^{N-1} [q_{j+1} - q_i]\, I(q_j)
    \approx \sum_{j=p}^{N-1} W_{ij}\, I(q_j)

for $q_p = q_i - \Delta q_v$ and $q_N = q_i + \Delta q_v$

**[Equation 6]**

where $W_{ij} = 0$ for $I_s$ when $j < p$ or $j > N-1$.

Solution 3
^^^^^^^^^^

**For** $\Delta q_v = \text{constant}$ **and** $\Delta q_u = \text{constant}$

In this case, the best way is to perform the integration of Equation 1
numerically for both slit height and slit width. However, the numerical
integration is imperfect unless a large number of iterations, say, at
least 10000 by 10000 for each element of the matrix $W$, is performed.
This is usually too slow for routine use.

An alternative approach is used in sasmodels which assumes
slit width << slit height. This method combines Solution 1 with the
numerical integration for the slit width. Then

.. math::

    I_s(q_i)
    &\approx \sum_{j=p}^{N-1} \sum_{k=-L}^L
            \left[\sqrt{q_{j+1}^2 - (q_i + (k\Delta q_v/L))^2}
                  - \sqrt{q_j^2 - (q_i + (k\Delta q_v/L))^2}\right]
            (\Delta q_v/L)\, I(q_j) \\
    &\approx \sum_{j=p}^{N-1} W_{ij}\,I(q_j)

**[Equation 7]**

for $q_p = q_i - \Delta q_v$ and $q_N = q_i + \Delta q_v$

where $W_{ij} = 0$ for $I_s$ when $j < p$ or $j > N-1$.

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

Pinhole Smearing
----------------

**This is the type of smearing normally encountered with data from synchrotron**
**SAXS cameras and SANS instruments.**

The pinhole smearing computation is performed in a similar fashion to the
slit-smeared case above except that the weight function used is a Gaussian. Thus
Equation 6 becomes

.. math::

    I_s(q_i)
    &\approx \sum_{j=0}^{N-1}[\operatorname{erf}(q_{j+1})
                - \operatorname{erf}(q_j)]\, I(q_j) \\
    &\approx \sum_{j=0}^{N-1} W_{ij}\, I(q_j)

**[Equation 8]**

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

2D Smearing
-----------

The 2D smearing computation is performed in a similar fashion to the 1D pinhole
smearing above except that the weight function used is a 2D elliptical Gaussian.
Thus

.. math::

  I_s(x_0,\, y_0)
  &= A\iint dx'dy'\,
     \exp \left[ -\left(\frac{(x'-x_0')^2}{2\sigma_{x_0'}^2}
                      + \frac{(y'-y_0')^2}{2\sigma_{y_0'}}\right)\right] I(x',\, y') \\
  &= A\sigma_{x_0'}\sigma_{y_0'}\iint dX dY\,
     \exp\left[-\frac{(X^2+Y^2)}{2}\right] I(\sigma_{x_0'}X x_0',\, \sigma_{y_0'} Y + y_0') \\
  &= A\sigma_{x_0'}\sigma_{y_0'}\iint dR d\Theta\,
     R\exp\left(-\frac{R^2}{2}\right) I(\sigma_{x_0'}R\cos\Theta + x_0',\, \sigma_{y_0'}R\sin\Theta+y_0')

**[Equation 9]**

In Equation 9, $x_0 = q \cos(\theta)$, $y_0 = q \sin(\theta)$, and
the primed axes are all in the coordinate rotated by an angle $\theta$ about
the $z$\ -axis (see the figure below) so that
$x'_0 = x_0 \cos(\theta) + y_0 \sin(\theta)$ and
$y'_0 = -x_0 \sin(\theta) + y_0 \cos(\theta)$.
Note that the rotation angle is zero for a $x$-$y$ symmetric
elliptical Gaussian distribution. The $A$ is a normalization factor.

.. figure:: resolution_2d_rotation.png

    Coordinate axis rotation for 2D resolution calculation.

Now we consider a numerical integration where each of the bins in $\theta$
and $R$ are *evenly* (this is to simplify the equation below) distributed
by $\Delta \theta$ and $\Delta R$ respectively, and it is further assumed
that $I(x',y')$ is constant within the bins. Then

.. math::

   I_s(x_0,\, y_0)
    &\approx A \sigma_{x'_0}\sigma_{y'_0}\sum_i^n
        \Delta\Theta\left[\exp\left(\frac{(R_i-\Delta R/2)^2}{2}\right)
                    - \exp\left(\frac{(R_i + \Delta R/2)^2}{2}\right)\right]
                    I(\sigma_{x'_0} R_i\cos\Theta_i+x'_0,\, \sigma_{y'_0}R_i\sin\Theta_i + y'_0) \\
    &\approx \sum_i^n W_i\, I(\sigma_{x'_0} R_i \cos\Theta_i + x'_0,\, \sigma_{y'_0}R_i\sin\Theta_i + y'_0)

**[Equation 10]**

Since the weighting factor on each of the bins is known, it is convenient to
transform $x'$-$y'$ back to $x$-$y$ coordinates (by rotating it
by $-\theta$ around the $z$\ -axis).

Then, for a polar symmetric smear

.. math::

    I_s(x_0,\, y_0) \approx \sum_i^n W_i\,
        I(x'\cos\theta - y'\sin\theta,\, x'sin\theta + y'\cos\theta)

**[Equation 11]**

where

.. math::

    x' &= \sigma_{x'_0} R_i \cos\Theta_i + x'_0 \\
    y' &= \sigma_{y'_0} R_i \sin\Theta_i + y'_0 \\
    x'_0 &= q = \sqrt{x_0^2 + y_0^2} \\
    y'_0 &= 0

while for a $x$-$y$ symmetric smear

.. math::

    I_s(x_0,\, y_0) \approx \sum_i^n W_i\, I(x',\, y')

**[Equation 12]**

where

.. math::

    x' &= \sigma_{x'_0} R_i \cos\Theta_i + x'_0 \\
    y' &= \sigma_{y'_0} R_i \sin\Theta_i + y'_0 \\
    x'_0 &= x_0 = q_x \\
    y'_0 &= y_0 = q_y


The current version of sasmodels uses Equation 11 for 2D smearing, assuming
that all the Gaussian weighting functions are aligned in the polar coordinate.

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

Weighting & Normalization
-------------------------

In all the cases above, the weighting matrix $W$ is calculated on the first
call to a smearing function, and includes ~60 $q$ values (finely and evenly
binned) below (>0) and above the $q$ range of data in order to smear all
data points for a given model and slit/pinhole size. The *Norm* factor is
found numerically with the weighting matrix and applied on the computation
of $I_s$.

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

*Document History*

| 2015-05-01 Steve King
| 2017-05-08 Paul Kienzle
