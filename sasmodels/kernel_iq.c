
/*
    ##########################################################
    #                                                        #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #   !!                                              !!   #
    #   !!  KEEP THIS CODE CONSISTENT WITH KERNELPY.PY  !!   #
    #   !!                                              !!   #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #                                                        #
    ##########################################################
*/

/*
The environment needs to provide the following #defines:

   USE_OPENCL is defined if running in opencl
   KERNEL declares a function to be available externally
   KERNEL_NAME is the name of the function being declared
   NPARS is the number of parameters in the kernel
   PARAMETER_DECL is the declaration of the parameters to the kernel.

       Cylinder:

           #define PARAMETER_DECL \
           double length; \
           double radius; \
           double sld; \
           double sld_solvent

       Note: scale and background are not included

       Multi-shell cylinder (10 shell max):

           #define PARAMETER_DECL \
           double num_shells; \
           double length; \
           double radius[10]; \
           double sld[10]; \
           double sld_solvent

   PARAMETER_CALL(var) is the declaration of a call to the kernel.

       Cylinder:

           #define PARAMETER_CALL(var) \
           var.length, \
           var.radius, \
           var.sld, \
           var.sld_solvent

       Multi-shell cylinder:
           #define PARAMETER_CALL(var) \
           var.num_shells, \
           var.length, \
           var.radius, \
           var.sld, \
           var.sld_solvent

   INVALID is a test for model parameters in the correct range

       Cylinder:

           #define INVALID(var) 0

       BarBell:

           #define INVALID(var) (var.bell_radius > var.radius)

       Model with complicated constraints:

           inline bool constrained(p1, p2, p3) { return expression; }
           #define INVALID(var) constrained(var.p1, var.p2, var.p3)

   IQ_FUNC could be Iq or Iqxy
   IQ_PARS could be q[i] or q[2*i],q[2*i+1]

Our design supports a limited number of polydispersity loops, wherein
we need to cycle through the values of the polydispersity, calculate
the I(q, p) for each combination of parameters, and perform a normalized
weighted sum across all the weights.  Parameters may be passed to the
underlying calculation engine as scalars or vectors, but the polydispersity
calculator treats the parameter set as one long vector.

Let's assume we have 6 parameters in the model, with two polydisperse::

    0: scale        {scl = constant}
    1: background   {bkg = constant}
    5: length       {l = vector of 30pts}
    4: radius       {r = vector of 10pts}
    3: sld          {s = constant/(radius**2*length)}
    2: sld_solvent  {s2 = constant}

This generates the following call to the kernel (where x stands for an
arbitrary value that is not used by the kernel evaluator):

    NPARS = 4  // scale and background are in all models
    problem {
        pd_par = {5, 4, x, x}         // parameters *radius* and *length* vary
        pd_length = {30, 10, 0, 0}    // *length* has more, so it is first
        pd_offset = {10, 0, x, x}     // *length* starts at index 10 in weights
        pd_stride = {1, 30, 300, 300} // cumulative product of pd length
        pd_isvol = {1, 1, x, x}       // true if weight is a volume weight
        par_offset = {2, 3, 303, 313}  // parameter offsets
        par_coord = {0, 3, 2, 1} // bitmap of parameter dependencies
        fast_coord_count = 2  // two parameters vary with *length* distribution
        fast_coord_index = {5, 3, x, x, x, x}
    }

    weight = { l0, .., l29, r0, .., r9}
    pars = { scl, bkg, l0, ..., l29, r0, r1, ..., r9,
             s[l0,r0], ... s[l0,r9], s[l1,r0], ... s[l29,r9] , s2}

    nq = 130
    q = { q0, q1, ..., q130, x, x }  # pad to 8 element boundary
    result = {r1, ..., r130, norm, vol, vol_norm, x, x, x, x, x, x, x}


The polydisperse parameters are stored in as an array of parameter
indices, one for each polydisperse parameter, stored in pd_par[n].
Non-polydisperse parameters do not appear in this array. Each polydisperse
parameter has a weight vector whose length is stored in pd_length[n],
The weights are stored in a contiguous vector of weights for all
parameters, with the starting position for the each parameter stored
in pd_offset[n].  The values corresponding to the weights are stored
together in a separate weights[] vector, with offset stored in
par_offset[pd_par[n]]. Polydisperse parameters should be stored in
decreasing order of length for highest efficiency.

We limit the number of polydisperse dimensions to MAX_PD (currently 4).
This cuts the size of the structure in half compared to allowing a
separate polydispersity for each parameter.  This will help a little
bit for models with large numbers of parameters, such as the onion model.

Parameters may be coordinated.  That is, we may have the value of one
parameter depend on a set of other parameters, some of which may be
polydisperse.  For example, if sld is inversely proportional to the
volume of a cylinder, and the length and radius are independently
polydisperse, then for each combination of length and radius we need a
separate value for the sld.  The caller must provide a coordination table
for each parameter containing the value for each parameter given the
value of the polydisperse parameters v1, v2, etc.  The tables for each
parameter are arranged contiguously in a vector, with offset[k] giving the
starting location of parameter k in the vector.  Each parameter defines
coord[k] as a bit mask indicating which polydispersity parameters the
parameter depends upon. Usually this is zero, indicating that the parameter
is independent, but for the cylinder example given, the bits for the
radius and length polydispersity parameters would both be set, the result
being a (#radius x #length) table, or maybe a (#length x #radius) table
if length comes first in the polydispersity table.

NB: If we can guarantee that a compiler and OpenCL driver are available,
we could instead create the coordination function on the fly for each
parameter, saving memory and transfer time, but requiring a C compiler
as part of the environment.

In ordering the polydisperse parameters by decreasing length we can
iterate over the longest dispersion weight vector first.  All parameters
coordinated with this weight vector (the 'fast' parameters), can be
updated with a simple increment to the next position in the parameter
value table.  The indices of these parameters is stored in fast_coord_index[],
with fast_coord_count being the number of fast parameters.  A total
of NPARS slots is allocated to allow for the case that all parameters
are coordinated with the fast index, though this will likely be mostly
empty.  When the fast increment count reaches the end of the weight
vector, then the index of the second polydisperse parameter must be
incremented, and all of its coordinated parameters updated.  Because this
operation is not in the inner loop, a slower algorithm can be used.

If there is no polydispersity we pretend that it is polydisperisty with one
parameter, pd_start=0 and pd_stop=1.  We may or may not short circuit the
calculation in this case, depending on how much time it saves.

The problem details structure can be allocated and sent in as an integer
array using the read-only flag.  This allows us to copy it once per fit
along with the weights vector, since features such as the number of
polydisperity elements per pd parameter or the coordinated won't change
between function evaluations.  A new parameter vector is sent for
each I(q) evaluation.

To protect against expensive evaluations taking all the GPU resource
on large fits, the entire polydispersity will not be computed at once.
Instead, a start and stop location will be sent, indicating where in the
polydispersity loop the calculation should start and where it should
stop.  We can do this for arbitrary start/stop points since we have
unwound the nested loop.  Instead, we use the same technique as array
index translation, using div and mod to figure out the i,j,k,...
indices in the virtual nested loop.

The results array will be initialized to zero for polydispersity loop
entry zero, and preserved between calls to [start, stop] so that the
results accumulate by the time the loop has completed.  Background and
scale will be applied when the loop reaches the end.  This does require
that the results array be allocated read-write, which is less efficient
for the GPU, but it makes the calling sequence much more manageable.

Scale and background cannot be coordinated with other polydisperse parameters

Cutoff paramater is basically used to restrict the region where integration
is peformed i.e. polydispersity hypercude is limitted spheres.

*/

#define MAX_PD 4  // MAX_PD is the max number of polydisperse parameters
#define PD_2N 16  // PD_2N is the size of the coordination step table

typedef struct {
    int pd_par[MAX_PD];     // index of the nth polydispersity variable
    int pd_length[MAX_PD];  // length of the nth polydispersity weight vector
    int pd_offset[MAX_PD];  // offset of pd weights in the par & weight vector
    int pd_stride[MAX_PD];  // stride to move to the next index at this level
    int pd_isvol[MAX_PD];   // True if parameter is a volume weighting parameter
    int par_offset[NPARS];  // offset of par values in the par & weight vector
    int par_coord[NPARS];   // polydispersity coordination bitvector
    int fast_coord_count;   // number of parameters coordinated with pd 1
    int fast_coord_index[NPARS]; // index of the fast coordination parameters
} ProblemDetails;

typedef struct {
    PARAMETER_DECL;
} ParameterBlock;

#define KERNEL_NAME test_Iq
#define FULL_KERNEL_NAME test_Iq
#define IQ_FUNC Iq

#define IQ_PARAMETERS ignored
#define IQ_FIXED_PARAMETER_DECLARATIONS const double scale, \
    const double background, \
    const double ignored
#define IQ_PARAMETER_DECLARATIONS double ignored
#define IQXY_KERNEL_NAME bessel_Iqxy
#define IQXY_PARAMETERS ignored
#define IQXY_FIXED_PARAMETER_DECLARATIONS const double scale, \
    const double background, \
    const double ignored
#define IQXY_PARAMETER_DECLARATIONS double ignored


void FULL_KERNEL_NAME(
    int nq,                 // number of q values
    global const ProblemDetails *problem,
    global const double *weights,
    global const double *pars,
    global const double *q, // nq q values, with padding to boundary
    global double *result,  // nq return values, again with padding
    const double cutoff,    // cutoff in the polydispersity weight product
    const int pd_start,     // where we are in the polydispersity loop
    const int pd_stop,      // where we are stopping in the polydispersity loop
    )
{

  // Storage for the current parameter values.  These will be updated as we
  // walk the polydispersity cube.
  local ParameterBlock local_pars;  // current parameter values
  const double *parvec = &local_pars;  // Alias named parameters with a vector

  local int offset[NPARS-2];

#if defined(USE_SHORTCUT_OPTIMIZATION)
  if (pd_length[0] == 1) {
    // Shouldn't need to copy!!
    for (int k=0; k < NPARS; k++) {
      parvec[k] = pars[k+2];  // skip scale and background
    }

    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i < nq; i++) {
    {
      const double scattering = IQ_FUNC(IQ_PARS, IQ_PARAMETERS);
      result[i] += pars[0]*scattering + pars[1];
    }
    return;
  }
#endif


  // Since we are no longer looping over the entire polydispersity hypercube
  // for each q, we need to track the normalization values for each q in a
  // separate work vector.
  double norm;   // contains sum over weights
  double vol; // contains sum over volume
  double norm_vol; // contains weights over volume

  // Initialize the results to zero
  if (pd_start == 0) {
    norm_vol = 0.0;
    norm = 0.0;
    vol = 0.0;

    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i < nq; i++) {
      result[i] = 0.0;
    }
  } else {
    //Pulling values from previous segment
    norm = result[nq];
    vol = result[nq+1];
    norm_vol = results[nq+2];
  }

  // Location in the polydispersity hypercube, one index per dimension.
  local int pd_index[PD_MAX];

  // Trigger the reset behaviour that happens at the end the fast loop
  // by setting the initial index >= weight vector length.
  pd_index[0] = pd_length[0];

  // need product of weights at every Iq calc, so keep product of
  // weights from the outer loops so that weight = partial_weight * fast_weight
  double partial_weight = NAN; // product of weight w4*w3*w2 but not w1
  double partial_volweight = NAN;
  double weight = 1.0;        // set to 1 in case there are no weights
  double vol_weight = 1.0;    // set to 1 in case there are no vol weights

  // Loop over the weights then loop over q, accumulating values
  for (int loop_index=pd_start; loop_index < pd_stop; loop_index++) {
    // check if indices need to be updated
    if (pd_index[0] >= pd_length[0]) {
      pd_index[0] = loop_index%pd_length[0];
      partial_weight = 1.0;
      partial_volweight = 1.0;
      for (int k=1; k < MAX_PD; k++) {
        pd_index[k] = (loop_index%pd_length[k])/pd_stride[k];
        const double wi = weights[pd_offset[0]+pd_index[0]];
        partial_weight *= wi;
        if (pd_isvol[k]) partial_volweight *= wi;
      }
      weight = partial_weight * weights[pd_offset[0]+pd_index[0]]
      for (int k=0; k < NPARS; k++) {
        int coord = par_coord[k];
        int this_offset = par_offset[k];
        int block_size = 1;
        for (int bit=0; bit < MAX_PD && coord != 0; bit++) {
          if (coord&1) {
              this_offset += block_size * pd_index[bit];
              block_size *= pd_length[bit];
          }
          coord /= 2;
        }
        offset[k] = this_offset;
        parvec[k] = pars[this_offset];
      }
    } else {
      pd_index[0] += 1;
      const double wi = weights[pd_offset[0]+pd_index[0]];
      weight = partial_weight*wi;
      if (pd_isvol[0]) vol_weight *= wi;
      for (int k=0; k < problem->fast_coord_count; k++) {
        parvec[ fast_coord_index[k]]
            = pars[offset[fast_coord_index[k]] + pd_index[0]];
      }
    }
    #ifdef INVALID
    if (INVALID(local_pars)) continue;
    #endif
    if (weight > cutoff) {
      norm += weight;
      vol += vol_weight * form_volume(VOLUME_PARAMETERS);
      norm_vol += vol_weight;

      #ifdef USE_OPENMP
      #pragma omp parallel for
      #endif
      for (int i=0; i < nq; i++) {
      {
        const double scattering = IQ_FUNC(IQ_PARS, IQ_PARAMETERS);
        //const double scattering = Iq(q[i], IQ_PARAMETERS);
        result[i] += weight*scattering;
      }
  }
  //Makes a normalization avialable for the next round
  result[nq] = norm;
  result[nq+1] = vol;
  result[nq+2] = norm_vol;

  //End of the PD loop we can normalize
  if (pd_stop == pd_stride[MAX_PD-1]) {}
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i < nq; i++) {
      if (vol*norm_vol != 0.0) {
        result[i] *= norm_vol/vol;
      }
        result[i] = pars[0]*result[i]/norm + pars[1];
    }
  }
}
}
}