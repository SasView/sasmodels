.. py:currentmodule:: sasmodels

.. _Writing_a_Plugin:

Writing a Plugin Model
======================

Overview
^^^^^^^^

In addition to the models provided with the sasmodels package, you are free to
create your own models. This document describes how to create plugin models
from first principles.

If you are using SasView and simply want to combine existing models into a new
plugin then use the "Add/Multiply Models" menu option instead.

Models can be of three types:

- A pure python model : Example -
  `broadpeak.py <https://github.com/SasView/sasmodels/blob/master/sasmodels/models/broad_peak.py>`_

- A python model with embedded C : Example -
  `sphere.py <https://github.com/SasView/sasmodels/blob/master/sasmodels/models/sphere.py>`_

- A python wrapper with separate C code : Example -
  `cylinder.py <https://github.com/SasView/sasmodels/blob/master/sasmodels/models/cylinder.py>`_,
  `cylinder.c <https://github.com/SasView/sasmodels/blob/master/sasmodels/models/cylinder.c>`_


When using SasView, plugin models should be saved to the SasView
*plugin_models* folder *C:\\Users\\{username}\\.sasview\\plugin_models*
(on Windows) or */Users/{username}/.sasview\\plugin_models* (on Mac).
The next time SasView is started it will compile the plugin and add
it to the list of *Plugin Models* in a FitPage.  Scripts can load
the models from anywhere.

The built-in modules are available in the *models* subdirectory
of the sasmodels package.  For SasView on Windows, these will
be found in *C:\\Program Files (x86)\\SasView\\sasmodels-data\\models*.
On Mac OSX, these will be within the application bundle as
*/Applications/SasView.app/Contents/Resources/sasmodels-data/models*.

Other models are available for download from the
`Model Marketplace <http://marketplace.sasview.org/>`_. You can contribute your
own models to the Marketplace as well.

Create New Model Files
^^^^^^^^^^^^^^^^^^^^^^

Copy the appropriate files to your plugin models directory (we recommend
using the examples above as templates) as mymodel.py (and mymodel.c, etc)
as required, where "mymodel" is the name for the model you are creating.

*Please follow these naming rules:*

- No capitalization and thus no CamelCase
- If necessary use underscore to separate words (i.e. barbell not BarBell or
  broad_peak not BroadPeak)
- Do not include "model" in the name (i.e. barbell not BarBellModel)


Edit New Model Files
^^^^^^^^^^^^^^^^^^^^

Model Contents
..............

The model interface definition is in the .py file.  This file contains:

- a **model name**:
   - this is the **name** string in the *.py* file
   - titles should be:

    - all in *lower* case
    - without spaces (use underscores to separate words instead)
    - without any capitalization or CamelCase
    - without incorporating the word "model"
    - examples: *barbell* **not** *BarBell*; *broad_peak* **not** *BroadPeak*;
      *barbell* **not** *BarBellModel*

- a **model title**:
   - this is the **title** string in the *.py* file
   - this is a one or two line description of the model, which will appear
     at the start of the model documentation and as a tooltip in the SasView GUI

- a **short description**:
   - this is the **description** string in the *.py* file
   - this is a medium length description which appears when you click
     *Description* on the model FitPage

- a **parameter table**:
   - this will be auto-generated from the *parameters* in the *.py* file

- a **long description**:
   - this is ReStructuredText enclosed between the r""" and """ delimiters
     at the top of the *.py* file
   - what you write here is abstracted into the SasView help documentation
   - this is what other users will refer to when they want to know what
     your model does; so please be helpful!

- a **definition** of the model:
   - as part of the **long description**

- a **formula** defining the function the model calculates:
   - as part of the **long description**

- an **explanation of the parameters**:
   - as part of the **long description**
   - explaining how the symbols in the formula map to the model parameters

- a **plot** of the function, with a **figure caption**:
   - this is automatically generated from your default parameters

- at least one **reference**:
   - as part of the **long description**
   - specifying where the reader can obtain more information about the model

- the **name of the author**
   - as part of the **long description**
   - the *.py* file should also contain a comment identifying *who*
     converted/created the model file

Models that do not conform to these requirements will *never* be incorporated
into the built-in library.


Composite Models
................

The simplest way to define a new model is to combine existing models.  For
example, *cyl_sphere.py* could contain::

    import os.path
    from sasmodels.core import load_model_info
    model_info = load_model_info('cylinder+sphere@hardsphere')
    model_info.name = os.path.basename(__file__).split('.')[0]

This defines a mixture of :ref:`cylinder` and :ref:`sphere` shapes, with
the spheres clustered closely enough to require an interaction effective
with the :ref:`hardsphere` structure factor.

The magic code at the end extracts the base filename, *cyl_sphere* from the
model file path and assigns it to the model name.

.. _Reparameterized_Models:

Reparameterized Models
......................

You can modify an existing model to use new parameters.  For example,
to create an ellipsoid constrained by volume::

    from numpy import inf
    from sasmodels.core import reparameterize
    parameters = [
        # name, units, default, [min, max], type, description
        ["volume", "Ang^3", 1e5, [0, inf], "volume", "ellipsoid volume"],
        ["eccentricity", "", 1, [0, inf], "volume", "polar:equatorial radius"],
    ]
    translation = """
        Re = cbrt(volume/eccentricity/M_4PI_3)
        radius_polar = eccentricity*Re
        radius_equatorial = Re  # python style comments allowed
        """
    model_info = reparameterize('ellipsoid', parameters, translation, __file__)

Here, *volume* and *eccentricity* are new parameters which replace the
*radius_polar* and *radius_equatorial* parameters in the :ref:`ellipsoid`
model.  The parameter properties are described below.  Since *volume* and
*eccentricity* are "volume" parameters they may be polydisperse.

*translation* is a set of equations to compute the underlying ellipsoid
parameters from the new parameters. *Re* is an intermediate value
introduced to make the equations easier to write.  Helper functions can
be included in a C file as *source=['filename.c', ...]* in the same directory
as the model file.

The new parameters replace *radius_polar* and *radius_equatorial* in the
parameter table.  To have more control over parameter placement, use an
*insert_after={...}* argmument to :func:`.core.reparameterize`.
For each insert location provide a list of new parameter names to insert
at that location.  For example, *{'': 'eccentricity,volume'}* inserts
them both at the beginning (before any parameter), whereas
*{'radius_polar': 'eccentricty', 'radius_equatorial': 'volume'}* will
place them after *radius_polar* and *radius_equatorial* in the final
parameter table, before deleting *radius_polar* and *radius_equatorial*.

Model Documentation
...................

The *.py* file starts with an r (for raw) and three sets of quotes
to start the doc string and ends with a second set of three quotes.
For example::

    r"""
    Definition
    ----------

    The 1D scattering intensity of the sphere is calculated in the following
    way (Guinier, 1955)

    .. math::

        I(q) = \frac{\text{scale}}{V} \cdot \left[
            3V(\Delta\rho) \cdot \frac{\sin(qr) - qr\cos(qr))}{(qr)^3}
            \right]^2 + \text{background}

    where *scale* is a volume fraction, $V$ is the volume of the scatterer,
    $r$ is the radius of the sphere and *background* is the background level.
    *sld* and *sld_solvent* are the scattering length densities (SLDs) of the
    scatterer and the solvent respectively, whose difference is $\Delta\rho$.

    You can included figures in your documentation, as in the following
    figure for the cylinder model.

    .. figure:: img/cylinder_angle_definition.jpg

        Definition of the angles for oriented cylinders.

    References
    ----------

    #. A Guinier, G Fournet, *Small-Angle Scattering of X-Rays*,
       John Wiley and Sons, New York, (1955)
    """

This is where the FULL documentation for the model goes (to be picked up by
the automatic documentation system).  Although it feels odd, you
should start the documentation immediately with the **definition**---the model
name, a brief description and the parameter table are automatically inserted
above the definition, and the a plot of the model is automatically inserted
before the **reference**.

Figures can be included using the *figure* command, with the name
of the *.png* file containing the figure and a caption to appear below the
figure.  Figure numbers will be added automatically.

See this `Sphinx cheat sheet <http://matplotlib.org/sampledoc/cheatsheet.html>`_
for a quick guide to the documentation layout commands, or the
`Sphinx Documentation <http://www.sphinx-doc.org/en/stable/>`_ for
complete details.

The model should include a **formula** written using LaTeX markup.
The example above uses the *math* command to make a displayed equation.  You
can also use *\$formula\$* for an inline formula. This is handy for defining
the relationship between the model parameters and formula variables, such
as the phrase "\$r\$ is the radius" used above.  The live demo MathJax
page `<http://www.mathjax.org/>`_ is handy for checking that the equations
will look like you intend.

Math layout uses the `amsmath <http://www.ams.org/publications/authors/tex/amslatex>`_
package for aligning equations (see amsldoc.pdf on that page for complete
documentation). You will automatically be in an aligned environment, with
blank lines separating the lines of the equation.  Place an ampersand before
the operator on which to align.  For example::

    .. math::

      x + y &= 1 \\
      y &= x - 1

produces

.. math::

      x + y &= 1 \\
      y &= x - 1

If you need more control, use::

    .. math::
        :nowrap:


Model Definition
................

Following the documentation string, there are a series of definitions::

    name = "sphere"  # optional: defaults to the filename without .py

    title = "Spheres with uniform scattering length density"

    description = """\
    P(q)=(scale/V)*[3V(sld-sld_solvent)*(sin(qr)-qr cos(qr))
                    /(qr)^3]^2 + background
        r: radius of sphere
        V: The volume of the scatter
        sld: the SLD of the sphere
        sld_solvent: the SLD of the solvent
    """

    category = "shape:sphere"

    single = True   # optional: defaults to True

    opencl = False  # optional: defaults to False

    structure_factor = False  # optional: defaults to False

**name = "mymodel"** defines the name of the model that is shown to the user.
If it is not provided it will use the name of the model file. The name must
be a valid variable name, starting with a letter and contains only letters,
numbers or underscore.  Spaces, dashes, and other symbols are not permitted.

**title = "short description"** is short description of the model which
is included after the model name in the automatically generated documentation.
The title can also be used for a tooltip.

**description = """doc string"""** is a longer description of the model. It
shows up when you press the "Description" button of the SasView FitPage.
It should give a brief description of the equation and the parameters
without the need to read the entire model documentation. The triple quotes
allow you to write the description over multiple lines. Keep the lines
short since the GUI will wrap each one separately if they are too long.
**Make sure the parameter names in the description match the model definition!**

**category = "shape:sphere"** defines where the model will appear in the
model documentation.  In this example, the model will appear alphabetically
in the list of spheroid models in the *Shape* category.

**single = True** indicates that the model can be run using single
precision floating point values.  Set it to False if the numerical
calculation for the model is unstable, which is the case for about 20 of
the built in models.  It is worthwhile modifying the calculation to support
single precision, allowing models to run up to 10 times faster.  The
section `Test_Your_New_Model`_  describes how to compare model values for
single vs. double precision so you can decide if you need to set
single to False.

**opencl = False** indicates that the model should not be run using OpenCL.
This may be because the model definition includes code that cannot be
compiled for the GPU (for example, goto statements).  It can also be used
for large models which can't run on most GPUs.  This flag has not been
used on any of the built in models; models which were failing were
streamlined so this flag was not necessary.

**structure_factor = True** indicates that the model can be used as a
structure factor to account for interactions between particles.  See
`Form_Factors`_ for more details.

**model_info = ...** lets you define a model directly, for example, by
loading and modifying existing models.  This is done implicitly by
:func:`.core.load_model_info`, which can create a mixture model
from a pair of existing models.  For example::

    from sasmodels.core import load_model_info
    model_info = load_model_info('sphere+cylinder')

See :class:`.modelinfo.ModelInfo` for details about the model
attributes that are defined.

Model Parameters
................

Next comes the parameter table.  For example::

    # pylint: disable=bad-whitespace, line-too-long
    #   ["name",        "units", default, [min, max], "type",    "description"],
    parameters = [
        ["sld",         "1e-6/Ang^2",  1, [-inf, inf], "sld",    "Layer scattering length density"],
        ["sld_solvent", "1e-6/Ang^2",  6, [-inf, inf], "sld",    "Solvent scattering length density"],
        ["radius",      "Ang",        50, [0, inf],    "volume", "Sphere radius"],
    ]
    # pylint: enable=bad-whitespace, line-too-long

**parameters = [["name", "units", default, [min,max], "type", "tooltip"],...]**
defines the parameters that form the model.

**Note: The order of the parameters in the definition will be the order of the
parameters in the user interface and the order of the parameters in Fq(), Iq(),
Iqac(), Iqabc(), radius_effective(), form_volume() and shell_volume().
And** *scale* **and** *background* **parameters are implicit to all models,
so they do not need to be included in the parameter table.**

- **"name"** is the name of the parameter shown on the FitPage.

  - the name must be a valid variable name, starting with a letter and
    containing only letters, numbers and underscore.

  - parameter names should follow the mathematical convention; e.g.,
    *radius_core* not *core_radius*, or *sld_solvent* not *solvent_sld*.

  - model parameter names should be consistent between different models,
    so *sld_solvent*, for example, should have exactly the same name
    in every model.

  - to see all the parameter names currently in use, type the following in the
    python shell/editor under the Tools menu::

       import sasmodels.list_pars
       sasmodels.list_pars.list_pars()

    *re-use* as many as possible!!!

  - use "name[n]" for multiplicity parameters, where *n* is the name of
    the parameter defining the number of shells/layers/segments, etc.

- **"units"** are displayed along with the parameter name

  - every parameter should have units; use "None" if there are no units.

  - **sld's should be given in units of 1e-6/Ang^2, and not simply
    1/Ang^2 to be consistent with the builtin models.  Adjust your formulas
    appropriately.**

  - fancy units markup is available for some units, including::

        Ang, 1/Ang, 1/Ang^2, 1e-6/Ang^2, degrees, 1/cm, Ang/cm, g/cm^3, mg/m^2

  - the list of units is defined in the variable *RST_UNITS* within
    `sasmodels/generate.py <https://github.com/SasView/sasmodels/tree/master/sasmodels/generate.py>`_

    - new units can be added using the macros defined in *doc/rst_prolog*
      in the sasmodels source.
    - units should be properly formatted using sub-/super-scripts
      and using negative exponents instead of the / operator, though
      the unit name should use the / operator for consistency.
    - please post a message to the SasView developers mailing list with your changes.

- **default** is the initial value for the parameter.

  - **the parameter default values are used to auto-generate a plot of
    the model function in the documentation.**

- **[min, max]** are the lower and upper limits on the parameter.

  - lower and upper limits can be any number, or *-inf* or *inf*.

  - the limits will show up as the default limits for the fit making it easy,
    for example, to force the radius to always be greater than zero.

  - these are hard limits defining the valid range of parameter values;
    polydisperity distributions will be truncated at the limits.

- **"type"** can be one of: "", "sld", "volume", or "orientation".

  - "sld" parameters can have magnetic moments when fitting magnetic models;
    depending on the spin polarization of the beam and the $q$ value being
    examined, the effective sld for that material will be used to compute the
    scattered intensity.

  - "volume" parameters are passed to Fq(), Iq(), Iqac(), Iqabc(), form_volume()
    and shell_volume(), and have polydispersity loops generated automatically.

  - "orientation" parameters are not passed, but instead are combined with
    orientation dispersity to translate *qx* and *qy* to *qa*, *qb* and *qc*.
    These parameters should appear at the end of the table with the specific
    names *theta*, *phi* and for asymmetric shapes *psi*, in that order.

Some models will have integer parameters, such as number of pearls in the
pearl necklace model, or number of shells in the multi-layer vesicle model.
The optimizers in BUMPS treat all parameters as floating point numbers which
can take arbitrary values, even for integer parameters, so your model should
round the incoming parameter value to the nearest integer inside your model
you should round to the nearest integer.  In C code, you can do this using:

.. code-block:: c

    static double
    Iq(double q, ..., double fp_n, ...)
    {
        int n = (int)(fp_n + 0.5);
        ...
    }

in python::

    def Iq(q, ..., n, ...):
        n = int(n+0.5)
        ...

Derivative based optimizers such as Levenberg-Marquardt will not work
for integer parameters since the partial derivative is always zero, but
the remaining optimizers (DREAM, differential evolution, Nelder-Mead simplex)
will still function.

Model Computation
.................

Models can be defined as pure python models, or they can be a mixture of
python and C models.  C models are run on the GPU if it is available,
otherwise they are compiled and run on the CPU.

Models are defined by the scattering kernel, which takes a set of parameter
values defining the shape, orientation and material, and returns the
expected scattering. Polydispersity and angular dispersion are defined
by the computational infrastructure.  Any parameters defined as "volume"
parameters are polydisperse, with polydispersity defined in proportion
to their value.  "orientation" parameters use angular dispersion defined
in degrees, and are not relative to the current angle.

Based on a weighting function $G(x)$ and a number of points $n$, the
computed value is

.. math::

     \hat I(q)
     = \frac{\int G(x) I(q, x)\,dx}{\int G(x) V(x)\,dx}
     \approx \frac{\sum_{i=1}^n G(x_i) I(q,x_i)}{\sum_{i=1}^n G(x_i) V(x_i)}

That is, the individual models do not need to include polydispersity
calculations, but instead rely on numerical integration to compute the
appropriately smeared pattern.

Each .py file also contains a function::

	def random():
	...

This function provides a model-specific random parameter set which shows model
features in the USANS to SANS range.  For example, core-shell sphere sets the
outer radius of the sphere logarithmically in `[20, 20,000]`, which sets the Q
value for the transition from flat to falling.  It then uses a beta distribution
to set the percentage of the shape which is shell, giving a preference for very
thin or very thick shells (but never 0% or 100%).  Using `-sets=10` in sascomp
should show a reasonable variety of curves over the default sascomp q range.
The parameter set is returned as a dictionary of `{parameter: value, ...}`.
Any model parameters not included in the dictionary will default according to
the code in the `_randomize_one()` function from sasmodels/compare.py.

Python Models
.............

.. note::

   Pure python models do not yet support direct computation of $<F(Q)^2>$ or
   $<F(Q)>^2$. Neither do they support orientational distributions or magnetism
   (use C models if these are required).

For pure python models, define the *Iq* function::

      import numpy as np
      from numpy import cos, sin, ...

      def Iq(q, par1, par2, ...):
          return I(q, par1, par2, ...)
      Iq.vectorized = True

The parameters *par1, par2, ...* are the list of non-orientation parameters
to the model in the order that they appear in the parameter table.
**Note that the auto-generated model file uses** *x* **rather than** *q*.

The *.py* file should import trigonometric and exponential functions from
numpy rather than from math.  This lets us evaluate the model for the whole
range of $q$ values at once rather than looping over each $q$ separately in
python.  With $q$ as a vector, you cannot use if statements, but must instead
do tricks like

::

     a = x*q*(q>0) + y*q*(q<=0)

or

::

     a = np.empty_like(q)
     index = q>0
     a[index] = x*q[index]
     a[~index] = y*q[~index]

which sets $a$ to $q \cdot x$ if $q$ is positive or $q \cdot y$ if $q$
is zero or negative. If you have not converted your function to use $q$
vectors, you can set the following and it will only receive one $q$
value at a time::

    Iq.vectorized = False

Return np.nan if the parameters are not valid (e.g., cap_radius < radius in
barbell).  If I(q; pars) is NaN for any $q$, then those parameters will be
ignored, and not included in the calculation of the weighted polydispersity.

Models should define *form_volume(par1, par2, ...)* where the parameter
list includes the *volume* parameters in order.  This is used for a weighted
volume normalization so that scattering is on an absolute scale.  For
solid shapes, the *I(q)* function should use *form_volume* squared
as its scale factor.  If *form_volume* is not defined, then the default
*form_volume = 1.0* will be used.

Hollow shapes, where the volume fraction of particle corresponds to the
material in the shell rather than the volume enclosed by the shape, must
also define a *shell_volume(par1, par2, ...)* function.  The parameters
are the same as for *form_volume*.  Here the *I(q)* function should use
*shell_volume* squared instead of *form_volume* squared so that the scale
parameter corresponds to the volume fraction of material in the sample.
The structure factor calculation needs the volume fraction of the filled
shapes for its calculation, so the volume fraction parameter in the model
is automatically scaled by *form_volume/shell_volume* prior to calling the
structure factor.

:ref:`Special_Functions` for scattering such as $3j_1(x)/x$
are available for both C and python.

Embedded C Models
.................

Like pure python models, inline C models need to define an *Iq* function::

    Iq = """
        return I(q, par1, par2, ...);
    """

This expands into the equivalent C code:

.. code-block:: c

    double Iq(double q, double par1, double par2, ...);
    double Iq(double q, double par1, double par2, ...)
    {
        return I(q, par1, par2, ...);
    }

*form_volume* defines the volume of the shape. As in python models, it
includes only the volume parameters.

*shell_volume* defines the volume of the shell for hollow shapes. As in
python models, it includes only the volume parameters.

**source=['fn.c', ...]** includes the listed C source files in the
program before *Iq* and *form_volume* are defined. These can include
files from the :ref:`Special_Functions` library.

*c_code* includes arbitrary C code into your kernel, which can be
handy for defining helper functions for *Iq* and *form_volume*. Note that
you can put the full function definition for *Iq* and *form_volume*
(include function declaration) into *c_code* as well, or put them into an
external C file and add that file to the list of sources.

Models are defined using double precision declarations for the
parameters and return values.  When a model is run using single
precision or long double precision, each variable is converted
to the target type, depending on the precision requested.

**Floating point constants must include the decimal point.**  This allows us
to convert values such as 1.0 (double precision) to 1.0f (single precision)
so that expressions that use these values are not promoted to double precision
expressions.  Some graphics card drivers are confused when functions
that expect floating point values are passed integers, such as 4*atan(1); it
is safest to not use integers in floating point expressions.  Even better,
use the builtin constant M_PI rather than 4*atan(1); it is faster and smaller!

The C model operates on a single $q$ value at a time.  The code will be
run in parallel across different $q$ values, either on the graphics card
or the processor.

Rather than returning NAN from Iq, you must provide a conditional expression
which evaluates to True if the parameters are valid and False if they
are not.  This is provided in the python model file as *valid = "expr"*,
where *expr* is a C expression.  For example::

    valid = "bell_radius >= radius && radius >= 0"

Structure Factors
.................

Structure factor calculations may need the underlying $<F(q)>$ and $<F^2(q)>$
rather than $I(q)$.  This is used to compute $\beta = <F(q)>^2/<F^2(q)>$ in
the decoupling approximation to the structure factor.

Instead of defining the *Iq* function, models can define *Fq* as
something like:

.. code-block:: c

    double Fq(double q, double *F1, double *F2, double par1, double par2, ...);
    double Fq(double q, double *F1, double *F2, double par1, double par2, ...)
    {
        // Polar integration loop over all orientations.
        ...
        *F1 = 1e-2 * total_F1 * contrast * volume;
        *F2 = 1e-4 * total_F2 * square(contrast * volume);
        return I(q, par1, par2, ...);
    }

If the volume fraction scale factor is built into the model (as occurs for
the vesicle model, for example), then scale *F1* by $\surd V_f$ so that
$\beta$ is computed correctly.

Structure factor calculations are not yet supported for oriented shapes.

Note: only available as a separate C file listed in *source*, or within
a *c_code* block within the python model definition file.

Oriented Shapes
...............

If the scattering is dependent on the orientation of the shape, then you
will need to include *orientation* parameters *theta*, *phi* and *psi*
at the end of the parameter table.  As described in the section
:ref:`orientation`, the individual $(q_x, q_y)$ points on the detector will
be rotated into $(q_a, q_b, q_c)$ points relative to the sample in its
canonical orientation with $a$-$b$-$c$ aligned with $x$-$y$-$z$ in the
laboratory frame and beam travelling along $-z$.

The oriented C model (oriented pure Python models are not supported)
is called using *Iqabc(qa, qb, qc, par1, par2, ...)* where
*par1*, etc. are the parameters to the model.  If the shape is rotationally
symmetric about *c* then *psi* is not needed, and the model is called
as *Iqac(qab, qc, par1, par2, ...)*.  In either case, the orientation
parameters are not included in the function call.

For 1D oriented shapes, an integral over all angles is usually needed for
the *Iq* function. Given symmetry and the substitution $u = \cos(\alpha)$,
$du = -\sin(\alpha)\,d\alpha$ this becomes

.. math::

    I(q) &= \frac{1}{4\pi} \int_{-\pi/2}^{pi/2} \int_{-pi}^{pi}
            F(q_a, q_b, q_c)^2 \sin(\alpha)\,d\beta\,d\alpha \\
        &= \frac{8}{4\pi} \int_{0}^{pi/2} \int_{0}^{\pi/2}
            F^2 \sin(\alpha)\,d\beta\,d\alpha \\
        &= \frac{8}{4\pi} \int_1^0 \int_{0}^{\pi/2} - F^2 \,d\beta\,du \\
        &= \frac{8}{4\pi} \int_0^1 \int_{0}^{\pi/2} F^2 \,d\beta\,du

for

.. math::

    q_a &= q \sin(\alpha)\sin(\beta) = q \sqrt{1-u^2} \sin(\beta) \\
    q_b &= q \sin(\alpha)\cos(\beta) = q \sqrt{1-u^2} \cos(\beta) \\
    q_c &= q \cos(\alpha) = q u

Using the $z, w$ values for Gauss-Legendre integration in "lib/gauss76.c", the
numerical integration is then:

.. code-block:: c

    double outer_sum = 0.0;
    for (int i = 0; i < GAUSS_N; i++) {
        const double cos_alpha = 0.5*GAUSS_Z[i] + 0.5;
        const double sin_alpha = sqrt(1.0 - cos_alpha*cos_alpha);
        const double qc = cos_alpha * q;
        double inner_sum = 0.0;
        for (int j = 0; j < GAUSS_N; j++) {
            const double beta = M_PI_4 * GAUSS_Z[j] + M_PI_4;
            double sin_beta, cos_beta;
            SINCOS(beta, sin_beta, cos_beta);
            const double qa = sin_alpha * sin_beta * q;
            const double qb = sin_alpha * cos_beta * q;
            const double form = Fq(qa, qb, qc, ...);
            inner_sum += GAUSS_W[j] * form * form;
        }
        outer_sum += GAUSS_W[i] * inner_sum;
    }
    outer_sum *= 0.25; // = 8/(4 pi) * outer_sum * (pi/2) / 4

The *z* values for the Gauss-Legendre integration extends from -1 to 1, so
the double sum of *w[i]w[j]* explains the factor of 4.  Correcting for the
average *dz[i]dz[j]* gives $(1-0) \cdot (\pi/2-0) = \pi/2$.  The $8/(4 \pi)$
factor comes from the integral over the quadrant.  With less symmetry (eg.,
in the bcc and fcc paracrystal models), then an integral over the entire
sphere may be necessary.

For simpler models which are rotationally symmetric a single integral
suffices:

.. math::

    I(q) &= \frac{1}{\pi}\int_{-\pi/2}^{\pi/2}
            F(q_{ab}, q_c)^2 \sin(\alpha)\,d\alpha/\pi \\
        &= \frac{2}{\pi} \int_0^1 F^2\,du

for

.. math::

    q_{ab} &= q \sin(\alpha) = q \sqrt{1 - u^2} \\
    q_c &= q \cos(\alpha) = q u


with integration loop::

    double sum = 0.0;
    for (int i = 0; i < GAUSS_N; i++) {
        const double cos_alpha = 0.5*GAUSS_Z[i] + 0.5;
        const double sin_alpha = sqrt(1.0 - cos_alpha*cos_alpha);
        const double qab = sin_alpha * q;
        const double qc = cos_alpha * q;
        const double form = Fq(qab, qc, ...);
        sum += GAUSS_W[j] * form * form;
    }
    sum *= 0.5; // = 2/pi * sum * (pi/2) / 2

Magnetism
.........

Magnetism is supported automatically for all shapes by modifying the
effective SLD of particle according to the Halpern-Johnson vector
describing the interaction between neutron spin and magnetic field.  All
parameters marked as type *sld* in the parameter table are treated as
possibly magnetic particles with magnitude *M0* and direction
*mtheta* and *mphi*.  Polarization parameters are also provided
automatically for magnetic models to set the spin state of the measurement.

For more complicated systems where magnetism is not uniform throughout
the individual particles, you will need to write your own models.
You should not mark the nuclear sld as type *sld*, but instead leave
them unmarked and provide your own magnetism and polarization parameters.
For 2D measurements you will need $(q_x, q_y)$ values for the measurement
to compute the proper magnetism and orientation, which you can implement
using *Iqxy(qx, qy, par1, par2, ...)*.

**Note: Magnetism is not supported in pure Python models.**


Problems with C models
......................

The graphics processor (GPU) in your computer is a specialized computer tuned
for certain kinds of problems.  This leads to strange restrictions that you
need to be aware of.  Your code may work fine on some platforms or for some
models, but then return bad values on other platforms.  Some examples of
particular problems:

  **(1) Code is too complex, or uses too much memory.** GPU devices only
  have a limited amount of memory available for each processor. If you run
  programs which take too much memory, then rather than running multiple
  values in parallel as it usually does, the GPU may only run a single
  version of the code at a time, making it slower than running on the CPU.
  It may fail to run on some platforms, or worse, cause the screen to go
  blank or the system to reboot.

  **(2) Code takes too long.** Because GPU devices are used for the computer
  display, the OpenCL drivers are very careful about the amount of time they
  will allow any code to run. For example, on OS X, the model will stop
  running after 5 seconds regardless of whether the computation is complete.
  You may end up with only some of your 2D array defined, with the rest
  containing random data. Or it may cause the screen to go blank or the
  system to reboot.

  **(3) Memory is not aligned**. The GPU hardware is specialized to operate
  on multiple values simultaneously. To keep the GPU simple the values in
  memory must be aligned with the different GPU compute engines. Not
  following these rules can lead to unexpected values being loaded into
  memory, and wrong answers computed. The conclusion from a very long and
  strange debugging session was that any arrays that you declare in your
  model should be a multiple of four. For example:

  .. code-block:: c

      double Iq(q, p1, p2, ...)
      {
          double vector[8];  // Only going to use seven slots, but declare 8
          ...
      }

The first step when your model is behaving strangely is to set
**single=False**. This automatically restricts the model to only run on the
CPU, or on high-end GPU cards. There can still be problems even on high-end
cards, so you can force the model off the GPU by setting **opencl=False**.
This runs the model as a normal C program without any GPU restrictions so
you know that strange results are probably from your code rather than the
environment. Once the code is debugged, you can compare your output to the
output on the GPU.

Although it can be difficult to get your model to work on the GPU, the reward
can be a model that runs 1000x faster on a good card.  Even your laptop may
show a 50x improvement or more over the equivalent pure python model.


.. _Form_Factors:

Form Factors
............

Away from the dilute limit you can estimate scattering including
particle-particle interactions using $I(q) = P(q)*S(q)$ where $P(q)$
is the form factor and $S(q)$ is the structure factor.  The simplest
structure factor is the *hardsphere* interaction, which
uses the effective radius of the form factor as an input to the structure
factor model.  The effective radius is the weighted average over all
values of the shape in polydisperse systems.

There can be many notions of effective radius, depending on the shape.  For
a sphere it is clearly just the radius, but for an ellipsoid of revolution
we provide average curvature, equivalent sphere radius, minimum radius and
maximum radius.  These options are listed as *radius_effective_modes* in
the python model defintion, and must be computed by the *radius_effective*
function which takes the *radius_effective_mode* parameter as an integer,
along with the various model parameters.  Unlike normal C/Python arrays,
the first mode is 1, the second is 2, etc.  Mode 0 indicates that the
effective radius from the shape is to be ignored in favour of the the
effective radius parameter in the structure factor model.


Consider the core-shell sphere, which defines the following effective radius
modes in the python model::

    radius_effective_modes = [
        "outer radius",
        "core radius",
    ]

and the following function in the C-file for the model:

.. code-block:: c

    static double
    radius_effective(int mode, double radius, double thickness)
    {
        switch (mode) {
            case 0: return radius + thickness;
            case 1: return radius;
            default: return 0.;
        }
    }

    static double
    form_volume(double radius, double thickness)
    {
        return M_4PI_3 * cube(radius + thickness);
    }

Given polydispersity over *(r1, r2, ..., rm)* in radius and *(t1, t2, ..., tn)*
in thickness, *radius_effective* is called over a mesh grid covering all
possible combinations of radius and thickness, with a single *(ri, tj)* pair
in each call. The weights each of these results according to the
polydispersity distributions and calls the structure factor with the average
effective radius.  Similarly, for *form_volume*.

Hollow models have an additional volume ratio which is needed to scale the
structure factor.  The structure factor uses the volume fraction of the filled
particles as part of its density estimate, but the scale factor for the
scattering intensity (as non-solvent volume fraction / volume) is determined
by the shell volume only.  Therefore the *shell_volume* function is
needed to compute the form:shell volume ratio, which then scales the
*volfraction* parameter prior to calling the structure factor calculator.
In the case of a hollow sphere, this would be:

.. code-block:: c

    static double
    shell_volume(double radius, double thickness)
    {
        double whole = M_4PI_3 * cube(radius + thickness);
        double core = M_4PI_3 * cube(radius);
        return whole - core;
    }

If *shell_volume* is not present, then *form_volume* and *shell_volume* are
assumed to be equal, and the shape is considered solid.

Unit Tests
..........

THESE ARE VERY IMPORTANT. Include at least one test for each model and
PLEASE make sure that the answer value is correct (i.e. not a random number).

::

    tests = [
        [{}, 0.2, 0.726362],
        [{"scale": 1., "background": 0., "sld": 6., "sld_solvent": 1.,
          "radius": 120., "radius_pd": 0.2, "radius_pd_n":45},
         0.2, 0.228843],
        [{"radius": 120., "radius_pd": 0.2, "radius_pd_n":45},
         0.1, None, None, 120., None, 1.],  # q, F, F^2, R_eff, V, form:shell
        [{"@S": "hardsphere"}, 0.1, None],
        [{"@S": "hardsphere"}, 0.1, None, {"S_eff(Q)": ...}],
    ]


**tests=[[{parameters}, q, Iq], ...]** is a list of lists.
Each list is one test and contains, in order:

- a dictionary of parameter values. This can be *{}* using the default
  parameters, or filled with some parameters that will be different from the
  default, such as *{"radius":10.0, "sld":4}*. Unlisted parameters will
  be given the default values.
- the input $q$ value or tuple of $(q_x, q_y)$ values.
- the output $I(q)$ or $I(q_x,q_y)$ expected of the model for the parameters
  and input value given.
- input and output values can themselves be lists if you have several
  $q$ values to test for the same model parameters.
- for testing effective radius, volume and form:shell volume ratio, use the
  extended form of the tests results, with *None, None, R_eff, V, V_r*
  instead of *Iq*.  This calls the kernel *Fq* function instead of *Iq*.
- for testing F and F^2 (used for beta approximation) do the same as the
  effective radius test, but include values for the first two elements,
  $<F(q)>$ and $<F^2(q)>$.
- for testing interaction between form factor and structure factor, specify
  the structure factor name in the parameters as *{"@S": "name", ...}* with
  the remaining list of parameters defined by the *P@S* interaction model.
- for examining the intermediate results computed by the model, use the
  extended form with *{"key": value}* pairs, one for each intermediate value
  returned from the model.  Starting with an empty *{}* will produce a
  complete list of computed intermediates.

.. _Test_Your_New_Model:

Test Your New Model
^^^^^^^^^^^^^^^^^^^

Minimal Testing
...............

In SasView 5.x, plugin models are loaded 'on-the-fly'. So to test a plugin
go to a FitPage (or create a new one with *Fitting* > *New Fit Page*),
select *Category* > *Plugin Models*, and then select your model from the
*Model Name* dropdown. Then click the *Calculate* button. If the model
compiles succesfully, two things will happen: the *Calculate* button will
change to *Show Plot*, and in the *Log Explorer* window something like
this will appear, possibly followed by a report of the unit test in your
plugin model:

  11:06:37 - INFO: make python model peak_voigt

If compilation of your plugin model fails, then Python traceback information
will appear in the Log Explorer:

  11:32:47 - INFO: make python model dummy_voigt

  11:32:47 - ERROR: Traceback (most recent call last):

If you look closely at traceback it should tell you why compilation is failing.
A very common reason is simply that a variable name has been misspelt. Look for
a line that begins *NameError*.

Recommended Testing
...................

**NB: For now, this more detailed testing is only possible if you have a
SasView build environment available!**

If the model compiles and runs, you can next run the unit tests that
you have added using the **test =** values.

From SasView, switch to the *Shell* tab and type the following::

    from sasmodels.model_test import run_one
    run_one("~/.sasview/plugin_models/model.py")

This should print::

    test_model_python (sasmodels.model_test.ModelTestCase) ... ok

To check whether single precision is good enough, type the following::

    from sasmodels.compare import main as compare
    compare("~/.sasview/plugin_models/model.py")

This will pop up a plot showing the difference between single precision
and double precision on a range of $q$ values.

These commands can also be run directly in the python interpreter:

    $ python -m sasmodels.model_test -v ~/.sasview/plugin_models/model.py
    $ python -m sasmodels.compare ~/.sasview/plugin_models/model.py

The options to compare are quite extensive; type the following for help::

    compare()

Options will need to be passed as separate strings.
For example to run your model with a random set of parameters::

    compare("-random", "-pars", "~/.sasview/plugin_models/model.py")

For the random models,

- *sld* will be in the range (-0.5,10.5),
- angles (*theta, phi, psi*) will be in the range (-180,180),
- angular dispersion will be in the range (0,45),
- polydispersity will be in the range (0,1)
- other values will be in the range (0, 2\ *v*), where *v* is the default
  parameter value.

Dispersion parameters *n*\, *sigma* and *type* will be set to random
values if "-poly" is included on the command line.  be unchanged from
demo so that run times are more predictable (polydispersity calculated
across multiple parameters can be very slow).

If your model has 2D orientation calculation, then you should also
test with::

    compare("-2d", "~/.sasview/plugin_models/model.py")

Check The Docs
^^^^^^^^^^^^^^

You can get a rough idea of how the documentation will look using the
following::

    compare("-help", "~/.sasview/plugin_models/model.py")

This does not use the same styling as the rest of the docs, but it will
allow you to check that your ReStructuredText and LaTeX formatting.
Here are some tools to help with the inevitable syntax errors:

- `Sphinx cheat sheet <http://matplotlib.org/sampledoc/cheatsheet.html>`_
- `Sphinx Documentation <http://www.sphinx-doc.org/en/stable/>`_
- `MathJax <http://www.mathjax.org/>`_
- `amsmath <http://www.ams.org/publications/authors/tex/amslatex>`_

There is also a neat online WYSIWYG ReStructuredText editor at
http://rst.ninjs.org\ .


Clean Lint - (Developer Version Only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**NB: For now we are not providing pylint with the installer version
of SasView; so unless you have a SasView build environment available,
you can ignore this section!**

Run the lint check with::

    python -m pylint --rcfile=extra/pylint.rc ~/.sasview/plugin_models/model.py

We are not aiming for zero lint just yet, only keeping it to a minimum.
For now, don't worry too much about *invalid-name*. If you really want a
variable name *Rg* for example because $R_g$ is the right name for the model
parameter then ignore the lint errors.  Also, ignore *missing-docstring*
for standard model functions *Iq*, *Iqac*, etc.

We will have delinting sessions at the SasView Code Camps, where we can
decide on standards for model files, parameter names, etc.

For now, you can tell pylint to ignore things.  For example, to align your
parameters in blocks::

    # pylint: disable=bad-whitespace,line-too-long
    #   ["name",                  "units", default, [lower, upper], "type", "description"],
    parameters = [
        ["contrast_factor",       "barns",    10.0,  [-inf, inf], "", "Contrast factor of the polymer"],
        ["bjerrum_length",        "Ang",       7.1,  [0, inf],    "", "Bjerrum length"],
        ["virial_param",          "1/Ang^2",  12.0,  [-inf, inf], "", "Virial parameter"],
        ["monomer_length",        "Ang",      10.0,  [0, inf],    "", "Monomer length"],
        ["salt_concentration",    "mol/L",     0.0,  [-inf, inf], "", "Concentration of monovalent salt"],
        ["ionization_degree",     "",          0.05, [0, inf],    "", "Degree of ionization"],
        ["polymer_concentration", "mol/L",     0.7,  [0, inf],    "", "Polymer molar concentration"],
        ]
    # pylint: enable=bad-whitespace,line-too-long

Don't put in too many pylint statements, though, since they make the code ugly.

Share Your Model!
^^^^^^^^^^^^^^^^^

Once compare and the unit test(s) pass properly and everything is done,
consider adding your model to the
`Model Marketplace <http://marketplace.sasview.org/>`_ so that others may use it!

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

*Document History*

| 2016-10-25 Steve King
| 2017-05-07 Paul Kienzle - Moved from sasview to sasmodels docs
| 2019-03-28 Paul Kienzle - Update docs for radius_effective and shell_volume
