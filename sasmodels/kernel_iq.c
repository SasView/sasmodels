
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
The environment needs to provide the following #defines
   USE_OPENCL is defined if running in opencl
   KERNEL declares a function to be available externally
   KERNEL_NAME is the name of the function being declared
   NPARS is the number of parameters in the kernel
   PARAMETER_TABLE is the declaration of the parameters to the kernel.

       Cylinder:

           #define PARAMETER_TABLE \
           double length; \
           double radius; \
           double sld; \
           double sld_solvent

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

Our design supports a limited number of polydispersity loops, wherein
we need to cycle through the values of the polydispersity, calculate
the I(q, p) for each combination of parameters, and perform a normalized
weighted sum across all the weights.  Parameters may be passed to the
underlying calculation engine as scalars or vectors, but the polydispersity
calculator treats the parameter set as one long vector.

The polydisperse parameters are stored in as an array of parameter
indices, one for each polydisperse parameter, stored in pd_par[n].
Non-polydisperse parameters do not appear in this array. Each polydisperse
parameter has a weight vector whose length is stored in pd_length[n].
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
value of the polydisperse parameters v1, v2, etc.   The tables for each
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
*/

#define MAX_PD 4  // MAX_PD is the max number of polydisperse parameters
#define PD_2N 16  // PD_2N is the size of the coordination step table

typedef struct {
    int pd_par[MAX_PD];     // index of the nth polydispersity variable
    int pd_length[MAX_PD];  // length of the nth polydispersity weight vector
    int pd_offset[MAX_PD];  // offset of pd weights in the par & weight vector
    int pd_stride[MAX_PD];  // stride to move to the next index at this level
    int par_offset[NPARS];  // offset of par values in the par & weight vector
    int par_coord[NPARS];   // polydispersity coordination bitvector
    int fast_coord_count;   // number of parameters coordinated with pd 1
    int fast_coord_index[NPARS]; // index of the fast coordination parameters
} ProblemDetails;

typedef struct {
    PARAMETER_DECL;
} ParameterBlock;


KERNEL
void KERNEL_NAME(
    int nq,                 // number of q values
    global const ProblemDetails *problem,
    global const double *weights,
    global const double *pars,
    global const double *q, // nq q values, with padding to boundary
    global double *result,  // nq return values, again with padding
    global double *work,    // 3*(nq+padding) space for normalization
    global int *int_work,   // nq+padding space for current position
    const double cutoff,    // cutoff in the polydispersity weight product
    const int pd_start,     // where we are in the polydispersity loop
    const int pd_stop,      // where we are stopping in the polydispersity loop
    )
{
  // Storage for the current parameter values.  These will be updated as we
  // walk the polydispersity cube.
  local ParameterBlock local_pars;
  const double *parvec = &local_pars;  // Alias named parameters with a vector

  // Since we are no longer looping over the entire polydispersity hypercube
  // for each q, we need to track the normalization values for each q in a
  // separate work vector.
  double *norm = work;   // contains sum over weights
  double *vol = norm + (nq + padding); // contains sum over volume
  double *norm_vol = vol + (nq + padding);

  // Initialize the results to zero
  if (pd_start == 0) {
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i < Nq; i++) {
      norm_vol[i] = 0.0;
      norm[i] = 0.0;
      vol[i] = 0.0;
      result[i] = 0.0;
    }
  }

  // Location in the polydispersity cube, one index per dimension.
  local int pd_index[PD_MAX];

  // Set the initial index greater than its vector length in order to
  // trigger the reset behaviour that happens at the end the fast loop.
  pd_index[0] = pd_length[0];

  // Loop over the weights then loop over q, accumulating values
  // par
  double partial_weight = NaN;
  for (int loop_index=pd_start; loop_index < pd_stop; loop_index++) {
    // check if indices need to be updated
    if (pd_index[0] >= pd_length[0]) {
      pd_index[0] = loop_index%pd_length[0];
      partial_weight = 1.0;
      for (int k=0; k < MAX_PD; k++) {
        pd_index[k] = (loop_index%pd_length[k])/pd_stride[k];
        partial_weight *= weights[pd_offset[k]+pd_index[k]];
      }
      weight = partial_weight * weights[pd_offset[0]+pd_index[0]]
      for (int k=0; k < NPARS; k++) {
        int coord = par_coord[k];
        int this_offset = 0;
        int block_size = 1;
        for (int bit=0; bit < MAX_PD && coord != 0; bit++) {
          if (coord&1) {
              this_offset += block_size * pd_index[bit];
              block_size *= pd_length[bit];
          }
          coord /= 2;
        }
        offset[k] = this_offset;
      }
    } else {
      pd_index[0] += 1;
      weight = partial_weight*weights[pd_offset[0]+pd_index[0]];
      for (int k=0; k < problem->fast_coord_count; k++) {
        parvec[ fast_coord_index[k]]
            = pars[offset[fast_coord_index[k]] + pd_index[0]];
      }
    }
    if (weight > cutoff) {
      const double vol_weight = VOLUME_WEIGHT_PRODUCT;
      const double weighted_vol = vol_weight*form_volume(VOLUME_PARAMTERS);
      #ifdef USE_OPENMP
      #pragma omp parallel for
      #endif
      for (int i=0; i < Nq; i++) {
      {
        const double scattering = Iq(qi, IQ_PARAMETERS);
        // allow kernels to exclude invalid regions by returning NaN
        if (!isnan(scattering)) {
          result[i] += weight*scattering;
          // can almost get away with only having one constant rather than
          // one constant per q.  Maybe want a "is_valid" test?
          norm[i] += weight;
          vol[i] += weighted_vol;
          norm_vol[i] += vol_weight;
        }
    }
  }

  if (pd_stop == pd_stride[MAX_PD-1]) {}
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i < Nq; i++) {
      if (vol[i]*norm_vol[i] != 0.0) {
        result[i] *= norm_vol[i]/vol[i];
      }
      //TODO: Ask Richard if scale and background may be corridanted parameters
      result[i] = scale*result[i]/norm[i]+background;
    }
  }
}
}
