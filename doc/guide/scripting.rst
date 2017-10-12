.. _Scripting_Interface:

*******************
Scripting Interface
*******************

Need some basic details here of how to load models and data via script, evaluate
them at given parameter values and run bumps fits.

The key functions are :func:`sasmodels.core.load_model` for loading the
model definition and compiling the kernel and
:func:`sasmodels.data.load_data` for calling sasview to load the data. Need
the data because that defines the resolution function and the q values to
evaluate. If there is no data, then use :func:`sasmodels.data.empty_data1D`
or :func:`sasmodels.data.empty_data2D` to create some data with a given $q$.

Using sasmodels through bumps
=============================

With the data and the model, you can wrap it in a *bumps* model with
class:`sasmodels.bumps_model.Model` and create an
class:`sasmodels.bump_model.Experiment` that you can fit with the *bumps*
interface. Here is an example from the *example* directory such as
*example/model.py*::

    import sys
    from bumps.names import *
    from sasmodels.core import load_model
    from sasmodels.bumps_model import Model, Experiment
    from sasmodels.data import load_data, set_beam_stop, set_top

    """ IMPORT THE DATA USED """
    radial_data = load_data('DEC07267.DAT')
    set_beam_stop(radial_data, 0.00669, outer=0.025)
    set_top(radial_data, -.0185)

    kernel = load_model("ellipsoid")

    model = Model(kernel,
        scale=0.08,
        radius_polar=15, radius_equatorial=800,
        sld=.291, sld_solvent=7.105,
        background=0,
        theta=90, phi=0,
        theta_pd=15, theta_pd_n=40, theta_pd_nsigma=3,
        radius_polar_pd=0.222296, radius_polar_pd_n=1, radius_polar_pd_nsigma=0,
        radius_equatorial_pd=.000128, radius_equatorial_pd_n=1, radius_equatorial_pd_nsigma=0,
        phi_pd=0, phi_pd_n=20, phi_pd_nsigma=3,
        )

    # SET THE FITTING PARAMETERS
    model.radius_polar.range(15, 1000)
    model.radius_equatorial.range(15, 1000)
    model.theta_pd.range(0, 360)
    model.background.range(0,1000)
    model.scale.range(0, 10)

    #cutoff = 0     # no cutoff on polydisperisity loops
    #cutoff = 1e-5  # default cutoff
    cutoff = 1e-3  # low precision cutoff
    M = Experiment(data=radial_data, model=model, cutoff=cutoff)
    problem = FitProblem(M)

Assume that bumps has been installed and the bumps command is available.
Maybe need to set the path to sasmodels/sasview
using *PYTHONPATH=path/to/sasmodels:path/to/sasview/src*.
To run the model use the *bumps* program::

    $ bumps example/model.py --preview

Note that bumps and sasmodels are included as part of the SasView
distribution.  On windows, bumps can be called from the cmd prompt
as follows::

    SasViewCom bumps.cli example/model.py --preview

Using sasmodels directly
========================

Bumps has a notion of parameter boxes in which you can set and retrieve
values.  Instead of using bumps, you can create a directly callable function
with :class:`sasmodels.direct_model.DirectModel`.  The resulting object *f*
will be callable as *f(par=value, ...)*, returning the $I(q)$ for the $q$
values in the data.  Polydisperse parameters use the same naming conventions
as in the bumps model, with e.g., radius_pd being the polydispersity associated
with radius.

Getting a simple function that you can call on a set of q values and return
a result is not so simple.  Since the time critical use case (fitting) involves
calling the function over and over with identical $q$ values, we chose to
optimize the call by only transfering the $q$ values to the GPU once at the
start of the fit.  We do this by creating a :class:`sasmodels.kernel.Kernel`
object from the :class:`sasmodels.kernel.KernelModel` returned from
:func:`sasmodels.core.load_model` using the
:meth:`sasmodels.kernel.KernelModel.make_kernel` method.  What it actual
does depends on whether it is running as a DLL, as OpenCL or as a pure
python kernel.  Once the kernel is in hand, we can then marshal a set of
parameters into a :class:`sasmodels.details.CallDetails` object and ship it to
the kernel using the :func:`sansmodels.direct_model.call_kernel` function.  An
example should help, *example/cylinder_eval.py*::

    from numpy import logspace
    from matplotlib import pyplot as plt
    from sasmodels.core import load_model
    from sasmodels.direct_model import call_kernel

    model = load_model('cylinder')
    q = logspace(-3, -1, 200)
    kernel = model.make_kernel([q])
    Iq = call_kernel(kernel, dict(radius=200.))
    plt.loglog(q, Iq)
    plt.show()

On windows, this can be called from the cmd prompt using sasview as::

    SasViewCom example/cylinder_eval.py
