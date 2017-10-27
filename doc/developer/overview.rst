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

Orientation and Numerical Integration
-------------------------------------

For 2d data from oriented anisotropic particles, the mean particle orientation is defined by angles $\theta$, $\phi$ and $\Psi$, which
are not in general the same as similarly named angles in many form factors. The wikipedia page on Euler angles 
(https://en.wikipedia.org/wiki/Euler_angles) lists the different conventions available. To quote: "Different authors may use different 
sets of rotation axes to define Euler angles, or different names for the same angles. Therefore, any discussion employing Euler angles
should always be preceded by their definition."

We are using the z-y-z convention with extrinsic rotations $\Psi-\theta-\phi$ for the particle orientation and $x-y-z$ convention with 
extrinsic rotations $\psi-\theta-\phi$ for jitter, with jitter applied before particle orientation.

For numerical integration within form factors etc. sasmodels is mostly using Gaussian quadrature with 20, 76 or 150 points depending on 
the model.  It also makes use of symmetries such as calculating only over one quadrant rather than the whole sphere.  There is often a 
U-substitution replacing $\theta$ with $cos(\theta)$ which changes the limits of integration from 0 to $\pi/2$ to 0 to 1 and also conveniently 
absorbs the $sin(\theta)$ scale factor in the integration.  This can cause confusion if checking equations to say include in a paper or thesis!
Most models use the same core kernel code expressed in terms of the rotated view (qa, qb, qc) for both the 1D and the 2D models, but there 
are also historical quirks such as the parallelepiped model, which has a useless transformation representing j0(a qa) as j0(b qa a/b).

Useful testing routines include -

:mod:`asymint` a direct implementation of the surface integral for certain models to get a more trusted value for the 1D integral using a reimplementation of the 2D kernel in python and mpmath 
(which computes math functions to arbitrary precision). It uses $\theta$ ranging from 0 to $\pi$ and $\phi$ ranging from 0 to $2\pi$.  It perhaps would benefit
from including the U-substitution for theta.

:mod:`check1d` uses sasmodels 1D integration and compares that with a rectangle distribution in $\theta$ and $\phi$, with $\theta$ limits set to
$\pm90/\sqrt(3)$ and $\phi$ limits set to $\pm180/\sqrt(3)$  
[The rectangle weight function uses the fact that the distribution width column is labelled sigma to decide 
that the 1-sigma width of a rectangular distribution needs to be multiplied by $\sqrt(3)$ to get the corresponding gaussian equivalent, 
or similar reasoning.]  This should rotate the sample through the entire $\theta-\phi$ 
surface according to the pattern that you see in jitter.py when you modify it to use 'rectangle' rather than 'gaussian' for its distribution 
without changing the viewing angle. When computing the dispersity integral, weights are scaled by abs(cos(dtheta)) to account for the points in 
phi getting closer together as dtheta increases.  This integrated dispersion is computed at a set of $(qx, qy)$ points $(q cos(\alpha), q sin(\alpha))$ 
at some angle $\alpha$ (currently angle=0) for each q used in the 1-D integration.  The individual q points should be equivalent to asymint rect-n 
when the viewing angle is set to (theta,phi,psi) = (90, 0, 0). Such tests can help to validate that 2d intensity is consistent with 1d models.
 
:mod:`sascomp -sphere=n` uses the identical rectangular distribution to compute the pattern of the qx-qy grid.  You can see from triaxial_ellipsoid
that there may be something wrong conceptually since the pattern is no longer circular when the view (theta,phi,psi) is not (90, phi, 0).  
check1d shows that it is different from the sasmodels 1D integral even when at theta=0, psi=0. Cross checking the values with asymint, 
the sasmodels 1D integral is better at low q, though for very large structures there are not enough points in the integration for sasmodels 1D 
to compute the high q 1D integral correctly. [Some of that may now be fixed?]

The :mod:`sascomp` utility can be used for 2d as well as 1d calculations to compare results for two sets of parameters or processor types, for example
these two oriented cylinders here should be equivalent.

:mod:`\./sascomp -2d cylinder theta=0 phi=0,90 theta_pd_type=rectangle phi_pd_type=rectangle phi_pd=10,1 theta_pd=1,10 length=500 radius=10`


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