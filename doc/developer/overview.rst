.. py:currentmodule:: sasmodels

***************************
Code Overview
***************************

Computational kernels
---------------------

* :mod:`core`
* :mod:`modelinfo`
* :mod:`kernel`
* :mod:`product`
* :mod:`mixture`

At the heart of *sasmodels* are the individual computational kernels.  These
functions take a particular $q$ value and a set of parameter values and
return the expected scattering for that $q$. The instructions for writing
a kernel are documented in :ref:`Writing_a_Plugin`.  The source code for
the kernels is stored in :mod:`models`.

The primary interface to the models is through :mod:`core`, which
provides functions for listing available models, loading the model definition
and compiling the model.  Use :func:`core.load_model` to load in
a model definition and compile it.  This makes use of
:func:`core.load_model_info` to load the model definition and
:func:`core.build_model` to turn it into a computational kernel model
:class:`kernel.KernelModel`.

The :class:`modelinfo.ModelInfo` class defines the properties
of the model including the available model parameters
:class:`modelinfo.ParameterTable` with individual parameter attributes
such as units and hard limits defined in :class:`modelinfo.Parameter`.

The :class:`product.ProductModel` and :class:`mixture.MixtureModel` classes
are derived models, created automatically for models with names like
"hardsphere*sphere" and "cylinder+sphere".

Data loaders
------------

* :mod:`data`

In order to test models a minimal set of data management routines is
provided in :mod:`data`.  In particular, it provides mock :class:`data.Data1D`
and :class:`data.Data2D` classes which mimic those classes in *SasView*.
The functions :func:`data.empty_data1D` and :func:`data.empty_data2D`
are handy for creating containers with a particular set of $q$, $\Delta q$
points which can later be evaluated, and :func:`data.plot_theory` to show
the result.  If *SasView* is available on the path then :func:`data.load_data`
can be used to load any data type defined in *SasView*.  The function
:func:`data.plot_data` can plot that data alone without the theory value.

Kernel execution
----------------

* :mod:`resolution`
* :mod:`resolution2d`
* :mod:`sesans`
* :mod:`weights`
* :mod:`details`
* :mod:`direct_model`
* :mod:`bumps_model`
* :mod:`sasview_model`

To execute a computational kernel at a particular set of $q$ values, the
use :meth:`kernel.KernelModel.make_kernel`, which returns a callable
:class:`kernel.Kernel` for that $q$ vector (or a pair of $q_x$, $q_y$
for 2-D datasets).

The calculated $q$ values should include the measured
data points as well as additional $q$ values required to properly compute the
$q$ resolution function.  The *Resolution* subclasses in :mod:`resolution`
define the *q_calc* attribute for this purpose.  These are
:class:`resolution.Perfect1D` for perfect resolution,
:class:`resolution.Pinhole1D` for the usual SANS pinhole aperture,
:class:`resolution.Slit1D` for the usual USANS slit aperture and
:class:`resolution2d.Pinhole2D` for 2-D pinhole data.
In addition, :class:`resolution2d.Slit2D` defines 1-D slit smeared data
for oriented samples, which require calculation at particular $q_x$ and
$q_y$ values instead of $|q|$ as is the case for orientationally averaged
USANS.  The :class:`sesans.SesansTransform` class acts like a 1-D resolution,
having a *q_calc* attribute that defines the calculated $q$ values for
the SANS models that get converted to spin-echo values by the
:meth:`sesnas.SesansTransform.apply` method.

Polydispersity is defined by :class:`weights.Dispersion` classes,
:class:`weights.RectangleDispersion`, :class:`weights.ArrayDispersion`,
:class:`weights.LogNormalDispersion`, :class:`weights.GaussianDispersion`,
:class:`weights.SchulzDispersion`.  The :func:`weights.get_weights`
function creates a dispersion object of the class matching
:attr:`weights.Dispersion.type`, and calls it with the current value
of the parameter.  This returns a vector of values and weights for a
weighted average polydispersity.

In order to call the :class:`kernel.Kernel`, the values and weights for
all parameters must be composed into a :class:`details.CallDetails` object.
This is a compact vector representation of the entire polydispersity
loop that can be passed easily to the kernel.  Additionally, the magnetic
parameters must be converted from polar to cartesian coordinates.  This
work is done by the :func:`details.make_kernel_args` function, which returns
values that can be sent directly to the kernel.  It uses
:func:`details.make_details` to set the details object and
:func:`details.convert_magnetism` for the coordinate transform.

In the end, making a simple theory function evaluation requires a lot of
setup. To make calling them a little easier, the *DirectModel* and
*BumpsModel* interfaces are provided.  See :ref:`Scripting_Interface`
for an example.

The :class:`direct_model.DirectModel` interface accepts a data object
and a kernel model.  Within the class,
the :meth:`direct_model.DataMixin._interpret_data` method is called to
query the data and set the resolution.
The :meth:`direct_model.DataMixin._calc_theory` takes a set of parameter
values, builds the kernel arguments, calls the kernel, and applies the
resolution function, returning the predicted value for the data $q$ values.
The :class:`bumps_model.Experiment` class is like the DirectModel class,
but it defines a Fitness class that can be handed directly to the
bumps optimization and uncertainty analysis program.

The :class:`sasview_model.SasviewModel` class defines a SasView 4.x
compatible interface to the sasmodels definitions, allowing sasmodels
to be used directly from SasView.  Over time the SasView shim should
disappear as SasView access the :class:`modelinfo.ModelInfo` and
computational kernels directly.

Kernel execution
----------------

* :mod:`kernelcl`
* :mod:`kerneldll`
* :mod:`kernelpy`
* :mod:`generate`


The kernel functions for the most part do not define polydispersity,
resolution or magnetism directly.  Instead sasmodels automatically
applies these, calling the computation kernel as needed.

The outermost loop is the resolution calculation.  For the 1-D case
this computes a single vector of $I(q)$ values and applies the convolution
to the resulting set.  Since the same $I(q)$ vector is used to compute the
convolution at each point, it can be precomputed before the convolution,
and so the convolution is reasonably efficient.  The 2-D case is not
that efficient, and instead recomputes the entire shifted/scaled set
of $q_x$, $q_y$ values many times, or very many times depending on the
accuracy requested.

Polydispersity is handled as a mesh over the polydisperse parameters.
This is the next level of the loop.  For C kernels run in a DLL or
using OpenCL, the polydisperisty loop is generated separately for each
model as C code.  Inside the polydispersity loop there is a loop over
the magnetic cross sections for magnetic models, updating the SLD
parameters with the effective magnetic SLD for that particular $q$
value. For OpenCL, each $q$ value loops over the
polydispersity mesh on a separate processor. For DLL, the outer loop
cycles through polydispersity, and the inner loop distributes q values
amongst the processors.  Like the DLL, the Python kernel execution
cycles over the polydisperse parameters and the magnetic cross sections,
calling the computation kernel with a vector of $q$ values.  Assuming
the kernel code accepts vectors, this can be fast enough (though it is
painfully slow if not vectorized).

Further details are provided in the next section,
:ref:`Calculator_Interface`

.. _orientation_developer:

Orientation and Numerical Integration
-------------------------------------

For 2d data from oriented anisotropic particles, the mean particle
orientation is defined by angles $\theta$, $\phi$ and $\Psi$, which are not
in general the same as similarly named angles in many form factors. The
wikipedia page on Euler angles (https://en.wikipedia.org/wiki/Euler_angles)
lists the different conventions available. To quote: "Different authors may
use different sets of rotation axes to define Euler angles, or different
names for the same angles. Therefore, any discussion employing Euler angles
should always be preceded by their definition."

We are using the $z$-$y$-$z$ convention with extrinsic rotations
$\Psi$-$\theta$-$\phi$ for the particle orientation and $x$-$y$-$z$
convention with extrinsic rotations $\Psi$-$\theta$-$\phi$ for jitter, with
jitter applied before particle orientation.

When computing the orientation dispersity integral, the weights for
the individual points depends on the map projection used to translate jitter
angles into latitude/longitude.  The choice of projection is set by
*sasmodels.generate.PROJECTION*, with the default *PROJECTION=1* for
equirectangular and *PROJECTION=2* for sinusoidal.  The more complicated
Guyou and Postel projections are not implemented. See jitter.draw_mesh
for details.

For numerical integration within form factors etc. sasmodels is mostly using
Gaussian quadrature with 20, 76 or 150 points depending on the model. It also
makes use of symmetries such as calculating only over one quadrant rather
than the whole sphere. There is often a U-substitution replacing $\theta$
with $cos(\theta)$ which changes the limits of integration from 0 to $\pi/2$
to 0 to 1 and also conveniently absorbs the $sin(\theta)$ scale factor in the
integration. This can cause confusion if checking equations to include in a
paper or thesis! Most models use the same core kernel code expressed in terms
of the rotated view ($q_a$, $q_b$, $q_c$) for both the 1D and the 2D models,
but there are also historical quirks such as the parallelepiped model, which
has a useless transformation representing $j_0(a q_a)$ as $j_0(b q_a a/b)$.

Useful testing routines include:

The *sascomp* utility is used to view and compare models with different
parameters and calculation engines. The usual case is to simply plot a
model that you are developing::

    python sascomp path/to/model.py

Once the obvious problems are addressed, check the numerical precision
across a variety of randomly generated inputs::

    python sascomp -engine=single,double path/to/model.py -sets=10

You can compare different parameter values for the same or different models.
For example when looking along the long axis of a cylinder ($\theta=0$),
dispersity in $\theta$ should be equivalent to dispersity in $\phi$
when $\phi=90$::

    python sascomp -2d cylinder theta=0 phi=0,90 theta_pd_type=rectangle \\
    phi_pd_type=rectangle phi_pd=10,1 theta_pd=1,10 length=500 radius=10

It turns out that they are not because the equirectangular map projection
weights the points by $\cos(\theta)$ so $\Delta\theta$ is not identical
to $\Delta\phi$.  Setting *PROJECTION=2* in :mod:`sasmodels.generate` helps
somewhat.  Postel would help even more in this case, though leading
to distortions elsewhere.  See :mod:`sasmodels.compare` for many more details.

*sascomp -ngauss=n* allows you to set the number of quadrature points used
for the 1D integration for any model.  For example, a carbon nanotube with
length 10 $\mu$\ m and radius 1 nm is not computed correctly at high $q$::

    python sascomp cylinder length=100000 radius=10 -ngauss=76,500 -double -highq

Note: ticket 702 gives some forms for long rods and thin disks that may
be used in place of the integration, depending on $q$, radius and length; if
the cylinder model is updated with these corrections then above call will show
no difference.

*explore/check1d.py* uses sasmodels 1D integration and compares that with a
rectangle distribution in $\theta$ and $\phi$, with $\theta$ limits set to
$\pm 90/\sqrt(3)$ and $\phi$ limits set to $\pm 180/\sqrt(3)$ [The rectangle
weight function uses the fact that the distribution width column is labelled
sigma to decide that the 1-$\sigma$ width of a rectangular distribution needs to
be multiplied by $\sqrt(3)$ to get the corresponding gaussian equivalent, or
similar reasoning.] This should rotate the sample through the entire
$\theta$-$\phi$ surface according to the pattern that you see in jitter.py when
you use 'rectangle' rather than 'gaussian' for its distribution without
changing the viewing angle. In order to match the 1-D pattern for an arbitrary
viewing angle on triaxial shapes, we need to integrate
over $\theta$, $\phi$ and $\Psi$.

*sascomp -sphere=n* uses the same rectangular distribution as check1d to
compute the pattern of the $q_x$-$q_y$ grid.  This ought to produce a
spherically symmetric pattern.

*explore/precision.py* investigates the accuracy of individual functions
on the different execution platforms.  This includes the basic special
functions as well as more complex expressions used in scattering.  In many
cases the OpenCL function in sasmodels will use a polynomial approximation
over part of the range to improve accuracy, as shown in::

    python explore/precision.py 3j1/x

*explore/realspace.py* allows you to set up a Monte Carlo simulation of your
model by sampling random points within and computing the 1D and 2D scattering
patterns.  This was used to check the core shell parallelepiped example.  This
is similar to the general sas calculator in sasview, though it uses different
code.

*sasmodels/jitter.py* is for exploring different options for handling
orientation and orientation dispersity.  It uses *sasmodels/guyou.py* to
generate the Guyou projection.

*explore/asymint.py* is a direct implementation of the 1D integration for
a number of different models.  It calculates the integral for a particular
$q$ using several different integration schemes, including mpmath with 100
bits of precision (33 digits), so you can use it to check the target values
for the 1D model tests.  This is not a general purpose tool; you will need to
edit the file to change the model parameters. It does not currently
apply the $u=cos(\theta)$ substitution that many models are using
internally so the 76-point gaussian quadrature may not match the value
produced by the eqivalent model from sasmodels.

*explore/symint.py* is for rotationally symmetric models (just cylinder for
now), so it can compute an entire curve rather than a single $q$ point.  It
includes code to compare the long cylinder approximation to cylinder.

*explore/rpa.py* is for checking different implementations of the RPA model
against calculations for specific blends.  This is a work in (suspended)
progress.

There are a few modules left over in *explore* that can probably be removed.
*angular_pd.py* was an early version of *jitter.py*.  *sc.py* and *sc.c*
was used to explore different calculations for paracrystal models; it has
been absorbed into *asymint.py*. *transform_angles.py* translates orientation
parameters from the SasView 3.x definition to sasmodels.

*explore/angles.py* generates code for the view and jitter transformations.
Keep this around since it may be needed if we add new projections.

Testing
-------

* :mod:`model_test`
* :mod:`compare`
* :mod:`compare_many`
* :mod:`rst2html`
* :mod:`list_pars`

Individual models should all have test values to make sure that the
evaluation is correct.  This is particularly important in the context
of OpenCL since sasmodels doesn't control the compiler or the hardware,
and since GPUs are notorious for preferring speed over precision.  The
tests can be run as a group using :mod:`model_test` as main::

    $ python -m sasmodels.model_test all

Individual models can be listed instead of *all*, which is convenient when
adding new models.

The :mod:`compare` module, usually invoked using *./sascomp* provides a
rich interface for exploring model accuracy, execution speed and parameter
ranges.  It also allows different models to be compared.
The :mod:`compare_many` module does batch comparisons, keeping a list of
the particular random seeds which lead to large differences in output
between different computing platforms.

The :mod:`rst2html` module provides tools for converting model docs to
html and viewing the html.  This is used by :mod:`compare` to display
the model description, such as::

    $ ./sascomp -html sphere

This makes debugging the latex much faster, though this may require
Qt in order for mathjax to work correctly.

When run as main, it can display arbitrary ReStructuredText files. E.g.,

::

    $ python -m sasmodels.rst2html doc/developer/overview.rst

This is handy for sorting out rst and latex syntax.  With some work
the results could be improved so that it recognizes sphinx roles
such as *mod*, *class* and *func*, and so that it uses the style sheets
from the sasmodels docs.

The :mod:`list_pars` module lists all instances of parameters across
all models.  This helps to make sure that similar parameters have
similar names across the different models.  With the verbose flag,
the particular models which use each named parameter are listed.


Model conversion
----------------

* :mod:`convert`
* :mod:`conversion_table`

Model definitions are not static.  As needs change or problems are found,
models may be updated with new model names or may be reparameterized
with new parameter definitions.  For example, in translating the
Teubner-Strey model from SasView 3.x to sasmodels, the definition
in terms of *drho*, *k*, *c1*, *c2*, *a2* and prefactor was replaced
by the defintion in terms of *volfraction_a*, *xi*, *d*, *sld_a* and
*sld_b*.  Within :mod:`convert`, the *_hand_convert_3_1_2_to_4_1*
function must be called when loading a 3.x model definition to update it to
4.1, and then the model should be further updated to 4.2, 5.0, and so on.
The :func:`convert.convert_model` function does this, using the conversion
table in :mod:`conversion_table` (which handled the major renaming from
SasView 3.x to sasmodels), and using the internal *_hand_convert* function
for the more complicated cases.

Other
-----

* :mod:`exception`
* :mod:`alignment`

The :func:`exception.annotate_exception` function annotates the current
exception with a context string, such as "while opening myfile.dat" without
adjusting the traceback.

The :mod:`alignment` module is unused.