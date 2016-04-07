
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
    global const ProblemDetails *problem,
    global const double *weights,
    global const double *values,
    global const double *q, // nq q values, with padding to boundary
    global double *result,  // nq+3 return values, again with padding
    const double cutoff     // cutoff in the polydispersity weight product
    )
{
  // Storage for the current parameter values.  These will be updated as we
  // walk the polydispersity cube.
  local ParameterBlock local_values;  // current parameter values
  double *pvec = (double *)(&local_values);  // Alias named parameters with a vector

  // Fill in the initial variables
  #ifdef USE_OPENMP
  #pragma omp parallel for
  #endif
  for (int k=0; k < NPARS; k++) {
    pvec[k] = values[problem->par_offset[k]];
  }

  // If it is the first round initialize the result to zero, otherwise
  // assume that the previous result has been passed back.
  // Note: doing this even in the monodisperse case in order to handle the
  // rare case where the model parameters are invalid and zero is returned.
  // So slightly increased cost for slightly smaller code size.
  if (pd_start == 0) {
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i < nq+1; i++) {
      result[i] = 0.0;
    }
  }

  // Monodisperse computation
  if (problem->num_active == 0) {
    #ifdef INVALID
    if (INVALID(local_values)) { return; }
    #endif

    const double norm = CALL_VOLUME(local_values);
    const double scale = values[0];
    const double background = values[1];
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    result[nq] = norm; // Total volume normalization
    for (int i=0; i < nq; i++) {
      double scattering = CALL_IQ(q, i, local_values);
      result[i] = (norm>0. ? scale*scattering/norm + background : background);
    }
    return;
  }

#if MAX_PD > 0
  //printf("Entering polydispersity from %d to %d\n", pd_start, pd_stop);
  // Since we are no longer looping over the entire polydispersity hypercube
  // for each q, we need to track the normalization values between calls.
  double norm = 0.0;

  // need product of weights at every Iq calc, so keep product of
  // weights from the outer loops so that weight = partial_weight * fast_weight
  double partial_weight = NAN; // product of weight w4*w3*w2 but not w1
  double spherical_correction = 1.0;  // cosine correction for latitude variation

  // Location in the polydispersity hypercube, one index per dimension.
  local int pd_index[MAX_PD];

  // Location of the coordinated parameters in their own sub-cubes.
  local int offset[NPARS];

  // Trigger the reset behaviour that happens at the end the fast loop
  // by setting the initial index >= weight vector length.
  const int fast_length = problem->pd_length[0];
  pd_index[0] = fast_length;

  // Loop over the weights then loop over q, accumulating values
  for (int loop_index=pd_start; loop_index < pd_stop; loop_index++) {
    // check if indices need to be updated
    if (pd_index[0] == fast_length) {
      //printf("should be here with %d active\n", problem->num_active);

      // Compute position in polydispersity hypercube
      for (int k=0; k < problem->num_active; k++) {
        pd_index[k] = (loop_index/problem->pd_stride[k])%problem->pd_length[k];
        //printf("pd_index[%d] = %d\n",k,pd_index[k]);
      }

      // Compute partial weights
      partial_weight = 1.0;
      //printf("partial weight %d: ", loop_index);
      for (int k=1; k < problem->num_active; k++) {
        double wi = weights[problem->pd_offset[k] + pd_index[k]];
        //printf("pd[%d]=par[%d]=%g ", k, problem->pd_par[k], wi);
        partial_weight *= wi;
      }
      //printf("\n");

      // Update parameter offsets in weight vector
      //printf("slow %d: ", loop_index);
      for (int k=0; k < problem->num_coord; k++) {
        int par = problem->par_coord[k];
        int coord = problem->pd_coord[k];
        int this_offset = problem->par_offset[par];
        int block_size = 1;
        for (int bit=0; coord != 0; bit++) {
          if (coord&1) {
              this_offset += block_size * pd_index[bit];
              block_size *= problem->pd_length[bit];
          }
          coord >>= 1;
        }
        offset[par] = this_offset;
        pvec[par] = values[this_offset];
        //printf("par[%d]=v[%d]=%g \n", k, offset[k], pvec[k]);
        // if theta is not coordinated with fast index, precompute spherical correction
        if (par == problem->theta_par && !(problem->par_coord[k]&1)) {
          spherical_correction = fmax(fabs(cos(M_PI_180*pvec[problem->theta_par])), 1e-6);
        }
      }
      //printf("\n");
    }

    // Increment fast index
    const double wi = weights[problem->pd_offset[0] + pd_index[0]++];
    double weight = partial_weight*wi;
    //printf("fast %d: ", loop_index);
    for (int k=0; k < problem->num_coord; k++) {
      if (problem->pd_coord[k]&1) {
        const int par = problem->par_coord[k];
        pvec[par] = values[offset[par]++];
        //printf("p[%d]=v[%d]=%g ", par, offset[par]-1, pvec[par]);
        // if theta is coordinated with fast index, compute spherical correction each time
        if (par == problem->theta_par) {
          spherical_correction = fmax(fabs(cos(M_PI_180*pvec[problem->theta_par])), 1e-6);
        }
      }
    }
    //printf("\n");

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

      #ifdef USE_OPENMP
      #pragma omp parallel for
      #endif
      for (int i=0; i < nq; i++) {
        const double scattering = CALL_IQ(q, i, local_values);
        result[i] += weight*scattering;
      }
    }
  }

  // Accumulate norm.
  result[nq] += norm;

  // End of the PD loop we can normalize
  if (pd_stop >= problem->total_pd) {
    const double scale = values[0];
    const double background = values[1];
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i < nq; i++) {
      result[i] = (norm>0. ? scale*result[i]/norm + background : background);
    }
  }
#endif // MAX_PD > 0
}
