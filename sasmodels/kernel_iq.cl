
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

#ifndef _PAR_BLOCK_ // protected block so we can include this code twice.
#define _PAR_BLOCK_

typedef struct {
#if MAX_PD > 0
    int32_t pd_par[MAX_PD];     // id of the nth polydispersity variable
    int32_t pd_length[MAX_PD];  // length of the nth polydispersity weight vector
    int32_t pd_offset[MAX_PD];  // offset of pd weights in the value & weight vector
    int32_t pd_stride[MAX_PD];  // stride to move to the next index at this level
#endif // MAX_PD > 0
    int32_t par_offset[NPARS];  // offset of par value blocks in the value & weight vector
    int32_t par_coord[NPARS];   // ids of the coordination parameters
    int32_t pd_coord[NPARS];    // polydispersity coordination bitvector
    int32_t num_active;         // number of non-trivial pd loops
    int32_t total_pd;           // total number of voxels in hypercube
    int32_t num_coord;          // number of coordinated parameters
    int32_t theta_par;          // id of spherical correction variable
} ProblemDetails;

typedef struct {
    PARAMETER_TABLE;
} ParameterBlock;
#endif


kernel
void KERNEL_NAME(
    int32_t nq,                 // number of q values
    const int32_t pd_start,     // where we are in the polydispersity loop
    const int32_t pd_stop,      // where we are stopping in the polydispersity loop
    global const ProblemDetails *details,
    global const double *weights,
    global const double *values,
    global const double *q, // nq q values, with padding to boundary
    global double *result,  // nq+3 return values, again with padding
    const double cutoff     // cutoff in the polydispersity weight product
    )
{
  // Storage for the current parameter values.  These will be updated as we
  // walk the polydispersity cube.
  ParameterBlock local_values;  // current parameter values
  double *pvec = (double *)(&local_values);  // Alias named parameters with a vector
  double norm;

  // who we are and what element we are working with
  const int q_index = get_global_id(0);

  // number of active loops
  const int num_active = details->num_active;

  // Fill in the initial variables
  for (int k=0; k < NPARS; k++) {
    pvec[k] = values[details->par_offset[k]];
  }

  // Monodisperse computation
  if (num_active == 0) {
    #ifdef INVALID
    if (INVALID(local_values)) { return; }
    #endif
    norm = CALL_VOLUME(local_values);

    double scale, background;
    scale = values[0];
    background = values[1];

    if (q_index < nq) {
      double scattering = CALL_IQ(q, q_index, local_values);
      result[q_index] = (norm>0. ? scale*scattering/norm + background : background);
    }
    return;
  }

#if MAX_PD > 0

  double this_result;

  //printf("Entering polydispersity from %d to %d\n", pd_start, pd_stop);
  // norm will be shared across all threads.

  // need product of weights at every Iq calc, so keep product of
  // weights from the outer loops so that weight = partial_weight * fast_weight
  double partial_weight; // product of weight w4*w3*w2 but not w1
  double spherical_correction; // cosine correction for latitude variation
  double weight; // product of partial_weight*w1*spherical_correction

  // Location in the polydispersity hypercube, one index per dimension.
  int pd_index[MAX_PD];

  // Location of the coordinated parameters in their own sub-cubes.
  int offset[NPARS];

  // Number of coordinated indices
  const int num_coord = details->num_coord;

  // Number of elements in the longest polydispersity loop
  const int fast_length = details->pd_length[0];

  // Trigger the reset behaviour that happens at the end the fast loop
  // by setting the initial index >= weight vector length.
  pd_index[0] = fast_length;

  // Default the spherical correction to 1.0 in case it is not otherwise set
  spherical_correction = 1.0;

  // Since we are no longer looping over the entire polydispersity hypercube
  // for each q, we need to track the result and normalization values between
  // calls.  This means initializing them to 0 at the start and accumulating
  // them between calls.
  norm = pd_start == 0 ? 0.0 : result[nq];
  if (q_index < nq) {
    this_result = pd_start == 0 ? 0.0 : result[q_index];
  }

  // Loop over the weights then loop over q, accumulating values
  for (int loop_index=pd_start; loop_index < pd_stop; loop_index++) {
    // check if fast loop needs to be reset
    if (pd_index[0] == fast_length) {
      //printf("should be here with %d active\n", num_active);

      // Compute position in polydispersity hypercube
      for (int k=0; k < num_active; k++) {
        pd_index[k] = (loop_index/details->pd_stride[k])%details->pd_length[k];
        //printf("pd_index[%d] = %d\n",k,pd_index[k]);
      }

      // need to compute the product of the weights.  If the vector were really
      // long, we could split the work into groups, with each thread taking
      // every nth weight, but there really is no call for it here.  We could
      // also do some clever pair-wise multiplication similar to parallel
      // prefix, but again simpler is probably faster since n is likely small.
      // Compute partial weights
      partial_weight = 1.0;
      //printf("partial weight %d: ", loop_index);
      for (int k=1; k < num_active; k++) {
        double wi = weights[details->pd_offset[k] + pd_index[k]];
        //printf("pd[%d]=par[%d]=%g ", k, details->pd_par[k], wi);
        partial_weight *= wi;
      }
      //printf("\n");

      // Update parameter offsets in weight vector
      //printf("slow %d: ", loop_index);
      for (int k=0; k < num_coord; k++) {
        int par = details->par_coord[k];
        int coord = details->pd_coord[k];
        int this_offset = details->par_offset[par];
        int block_size = 1;
        for (int bit=0; coord != 0; bit++) {
          if (coord&1) {
              this_offset += block_size * pd_index[bit];
              block_size *= details->pd_length[bit];
          }
          coord >>= 1;
        }
        offset[par] = this_offset;
        pvec[par] = values[this_offset];
        //printf("par[%d]=v[%d]=%g \n", k, offset[k], pvec[k]);
        // if theta is not coordinated with fast index, precompute spherical correction
        if (par == details->theta_par && !(details->par_coord[k]&1)) {
          spherical_correction = fmax(fabs(cos(M_PI_180*pvec[details->theta_par])), 1.e-6);
        }
      }
      //printf("\n");
    }

    // Update fast parameters
    //printf("fast %d: ", loop_index);
    for (int k=0; k < num_coord; k++) {
      if (details->pd_coord[k]&1) {
        const int par = details->par_coord[k];
        pvec[par] = values[offset[par]++];
        //printf("p[%d]=v[%d]=%g ", par, offset[par]-1, pvec[par]);
        // if theta is coordinated with fast index, compute spherical correction each time
        if (par == details->theta_par) {
          spherical_correction = fmax(fabs(cos(M_PI_180*pvec[details->theta_par])), 1.e-6);
        }
      }
    }
    //printf("\n");

    // Increment fast index
    const double wi = weights[details->pd_offset[0] + pd_index[0]];
    weight = partial_weight*wi;
    pd_index[0]++;

    #ifdef INVALID
    if (INVALID(local_values)) continue;
    #endif

    // Accumulate I(q)
    // Note: weight==0 must always be excluded
    if (weight > cutoff) {
      // spherical correction has some nasty effects when theta is +90 or -90
      // where it becomes zero.  If the entirety of the correction
      weight *= spherical_correction;
      norm += weight * CALL_VOLUME(local_values);

      const double scattering = CALL_IQ(q, q_index, local_values);
      this_result += weight*scattering;
    }
  }

  if (q_index < nq) {
    if (pd_stop >= details->total_pd) {
      // End of the PD loop we can normalize
      double scale, background;
      scale = values[0];
      background = values[1];
      result[q_index] = (norm>0. ? scale*this_result/norm + background : background);
    } else {
      // Partial result, so remember it but don't normalize it.
      result[q_index] = this_result;
    }

    // Remember the updated norm.
    if (q_index == 0) result[nq] = norm;
  }

#endif // MAX_PD > 0
}
