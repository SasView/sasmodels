.. currentmodule:: sasmodels

.. _Calculator_Interface:

Calculator Interface
====================

This document describes the layer between the form factor kernels and the
model calculator which implements the dispersity and magnetic SLD
calculations.  There are three separate implementations of this layer,
:mod:`kernelcl` for OpenCL, which operates on a single Q value at a time,
:mod:`kerneldll` for the DLL, which loops over a vector of Q values, and
:mod:`kernelpy` for python models which operates on vector Q values.

Each implementation provides three different calls *Iq*, *Iqxy* and *Imagnetic*
for 1-D, 2-D and 2-D magnetic kernels respectively. The C code is defined
in *kernel_iq.c*, with the minor differences between OpenCL and DLL handled
by #ifdef statements.

The kernel call looks as follows::

  kernel void KERNEL_NAME(
      int nq,                  // Number of q values in the q vector
      int pd_start,            // Starting position in the dispersity loop
      int pd_stop,             // Ending position in the dispersity loop
      ProblemDetails *details, // dispersity info
      double *values,          // Value and weights vector
      double *q,               // q or (qx,qy) vector
      double *result,          // returned I(q), with result[nq] = pd_weight
      double cutoff)           // dispersity weight cutoff

The details for OpenCL and the python loop are slightly different, but these
data structures are common.

*nq* indicates the number of q values that will be calculated.

The *pd_start* and *pd_stop* parameters set the range of the dispersity
loop to compute for the current kernel call.   Give a dispersity
calculation with 30 weights for length and 30 weights for radius for example,
there are a total of 900 calls to the form factor required to compute the
kernel.  These might be done in chunks of 100, so the first call would start
at zero and stop after 100 iterations.  The second call would then have to set
the length index to 3 and the radius index to 10 for a position of 3*30+10=100,
and could then proceed to position 200.  This allows us to interrupt the
calculation in the middle of a long dispersity loop without having to
do special tricks with the C code.  More importantly, it stops the OpenCL
kernel in a reasonable time; because the GPU is used by the operating
system to show its windows, if a GPU kernel runs too long then it will be
automatically killed and no results will be returned to the caller.

The *ProblemDetails* structure is a direct map of the
:class:`details.CallDetails` buffer.  This indicates which parameters have
dispersity, and where in the values vector the values and weights can be
found.  For each parameter with dispersity there is a parameter id, the length
of the dispersity loop for that parameter, the offset of the parameter
values in the pd value and pd weight vectors and the 'stride' from one index
to the next, which is used to translate between the position in the
dispersity loop and the particular parameter indices.  The *num_eval*
field is the total size of the dispersity loop.  *num_weights* is the
number of elements in the pd value and pd weight vectors.  *num_active* is
the number of non-trivial pd loops (parameters with dispersity should be ordered
by decreasing pd vector length, with a length of 1 meaning no dispersity).
Oriented objects in 2-D need a cos(theta) spherical correction on the angular
variation in order to preserve the 'surface area' of the weight distribution.
*theta_par* is the id of the polar coordinate parameter if there is one.


The *values* vector consists of the fixed values for the model plus pd value
and pd weight vectors.  There are *NUM_VALUES* fixed values for the model,
which includes the initial two slots for scale and background, the NUM_PARS
values for the kernel parameters, three slots for the applied magnetism, and
three slots for each of the SLD parameters for the sample magnetiziation
*(Mx, My, Mz)*.  Sample magnetization is translated from *(M, theta, phi)*
to *(Mx, My, Mz)* before the kernel is called.   After the fixed values comes
the pd value vector, with the dispersity values for each parameter
stacked one after the other.  The order isn't important since the location
for each parameter is stored in the *pd_offset* field of the *ProblemDetails*
structure, but the values do need to be contiguous.  After *num_weights*
values, the pd weight vector is stored, with the same configuration as the
pd value vector.  Note that the pd vectors can contain values that are not
in the dispersity loop; this is used by :class:`mixture.MixtureKernel`
to make it easier to call the various component kernels.

The *q* vector contains one value for each q for *Iq* kernels, or a pair
of *(qx,qy)* values for each q for *Iqxy* and *Imagnetic* kernels.  OpenCL
pads all vectors to 32 value boundaries just to be safe, including the
*values* vector and the the *results* vector.

The *results* vector contains one slot for each of the *nq* values, plus
one extra slot at the end for the weight normalization accumulated across
all points in the dispersity mesh.  This is required when the dispersity
loop is broken across several kernel calls.

*cutoff* is a importance cutoff so that points which contribute negligibly
to the total scattering can be skipped without calculating them.

:func:`generate.make_source` defines the following C macros:

- USE_OPENCL is defined if running in opencl
- MAX_PD is the maximum depth of the dispersity loop [model specific]
- NUM_PARS is the number of parameter values in the kernel.  This may be
  more than the number of parameters if some of the parameters are vector
  values.
- NUM_VALUES is the number of fixed values, which defines the offset in the
  value list to the dispersity value and weight vectors.
- NUM_MAGNETIC is the number of magnetic SLDs
- MAGNETIC_PARS is a comma separated list of the magnetic SLDs, indicating
  their locations in the values vector.
- KERNEL_NAME is the name of the function being declared
- PARAMETER_TABLE is the declaration of the parameters to the kernel:

    Cylinder::

        #define PARAMETER_TABLE \
        double length; \
        double radius; \
        double sld; \
        double sld_solvent;

    Note: scale and background are never included

    Multi-shell cylinder (10 shell max)::

        #define PARAMETER_TABLE \
        double num_shells; \
        double length; \
        double radius[10]; \
        double sld[10]; \
        double sld_solvent;

- CALL_IQ(q, i, var) is the declaration of a call to the kernel:

    Cylinder::

        #define CALL_IQ(q, i, var) Iq(q[i], \
        var.length, \
        var.radius, \
        var.sld, \
        var.sld_solvent)

    Multi-shell cylinder::

        #define CALL_IQ(q, i, var) Iq(q[i], \
        var.num_shells, \
        var.length, \
        var.radius, \
        var.sld, \
        var.sld_solvent)

    Cylinder2D::

        #define CALL_IQ(q, i, var) Iqxy(qa, qc, \
        var.length, \
        var.radius, \
        var.sld, \
        var.sld_solvent)

- CALL_VOLUME(var) is similar, but for calling the form volume::

        #define CALL_VOLUME(var) \
        form_volume(var.length, var.radius)

There is an additional macro that can be defined within the model.c file:

- INVALID(var) is a test for model parameters in the correct range:

    Cylinder::

        #define INVALID(var) 0

    BarBell::

        #define INVALID(var) (var.bell_radius < var.radius)

    Model with complicated constraints::

        inline bool constrained(p1, p2, p3) { return expression; }
        #define INVALID(var) constrained(var.p1, var.p2, var.p3)

Our design supports a limited number of dispersity loops, wherein
we need to cycle through the values of the dispersity, calculate
the I(q, p) for each combination of parameters, and perform a normalized
weighted sum across all the weights.  Parameters may be passed to the
underlying calculation engine as scalars or vectors, but the dispersity
calculator treats the parameter set as one long vector.

Let's assume we have 8 parameters in the model, two of which allow dispersity.
Since this is a 1-D model the orientation parameters won't be used::

    0: scale        {scl = constant}
    1: background   {bkg = constant}
    2: radius       {r = vector of 10pts}
    3: length       {l = vector of 30pts}
    4: sld          {s1 = constant/(radius**2*length)}
    5: sld_solvent  {s2 = constant}
    6: theta        {not used}
    7: phi          {not used}

This generates the following call to the kernel.  Note that parameters 4 and
5 are treated as having dispersity even though they don't --- this is because
it is harmless to do so and simplifies the looping code::

    MAX_PD = 4
    NUM_PARS = 6          // kernel parameters only
    NUM_VALUES = 17       // all values, including scale and background
    NUM_MAGNETIC = 2      // two parameters might be magnetic
    MAGNETIC_PARS = 4, 5  // they are sld and sld_solvent

    details {
        pd_par = {3, 2, 4, 5}         // parameters *radius* and *length* vary
        pd_length = {30, 10, 1, 1}    // *length* has more, so it is first
        pd_offset = {10, 0, 31, 32}   // *length* starts at index 10 in weights
        pd_stride = {1, 30, 300, 300} // cumulative product of pd length
        num_eval = 300   // 300 values in the dispersity loop
        num_weights = 42 // 42 values in the pd vector
        num_active = 2   // only the first two pd are active
        theta_var =  6   // spherical correction
    }

    values = { scl, bkg,                                  // universal
               r, l, s1, s2, theta, phi,                  // kernel pars
               in spin, out spin, spin angle,             // applied magnetism
               mx s1, my s1, mz s1, mx s2, my s2, mz s2,  // magnetic slds
               r0, .., r9, l0, .., l29, s, s2,            // pd values
               r0, .., r9, l0, .., l29, s, s2}            // pd weights

    nq = 130
    q = { q0, q1, ..., q130, x, x }  # pad to 8 element boundary
    result = {r1, ..., r130, pd_norm, x }

The dispersity parameters are stored in as an array of parameter
indices, one for each parameter, stored in pd_par[n]. Parameters which do
not support dispersity do not appear in this array. Each dispersity
parameter has a weight vector whose length is stored in pd_length[n].
The weights are stored in a contiguous vector of weights for all
parameters, with the starting position for the each parameter stored
in pd_offset[n].  The values corresponding to the weights are stored
together in a separate weights[] vector, with offset stored in
par_offset[pd_par[n]]. Dispersity parameters should be stored in
decreasing order of length for highest efficiency.

We limit the number of dispersity dimensions to MAX_PD (currently 4),
though some models may have fewer if they have fewer dispersity
parameters.  The main reason for the limit is to reduce code size.
Each additional dispersity parameter requires a separate dispersity
loop.  If more than 4 levels of dispersity are needed, then we need to
switch to a monte carlo importance sampling algorithm with better
performance for high-dimensional integrals.

Constraints between parameters are not supported.  Instead users will
have to define a new model with the constraints built in by making a
copy of the existing model.  Mac provides OpenCL and we are supplying
the tinycc compiler for Windows so this won't be a complete limitation,
but it is rather inconvenient.  The process could perhaps be automated
so that there is no code copying required, just an alternate CALL_IQ
macro that implements the constraints.  Think carefully about constraining
theta since the polar coordinates normalization is tied to this parameter.

If there is no dispersity we pretend that we have a disperisty mesh over
a single parameter with a single point in the distribution, giving
pd_start=0 and pd_stop=1.

The problem details structure could be allocated and sent in as an integer
array using the read-only flag.  This would allow us to copy it once per fit
along with the weights vector, since features such as the number of
disperity points per pd parameter won't change between function evaluations.
A new parameter vector must be sent for each I(q) evaluation.  This is
not currently implemented, and would require some resturcturing of
the :class:`sasview_model.SasviewModel` interface.

The results array will be initialized to zero for dispersity loop
entry zero, and preserved between calls to [start, stop] so that the
results accumulate by the time the loop has completed.  Background and
scale will be applied when the loop reaches the end.  This does require
that the results array be allocated read-write, which is less efficient
for the GPU, but it makes the calling sequence much more manageable.

For accuracy we may want to introduce Kahan summation into the integration::

    double accumulated_error = 0.0;
    ...
    #if USE_KAHAN_SUMMATION
        const double y = next - accumulated_error;
        const double t = ret + y;
        accumulated_error = (t - ret) - y;
        ret = t;
    #else
        ret += next;
    #endif

This will require accumulated error for each I(q) value to be preserved
between kernel calls to implement this fully.  The *kernel_iq.c* code, which
loops over q for each parameter set in the dispersity loop, will also need
the accumulation vector.
