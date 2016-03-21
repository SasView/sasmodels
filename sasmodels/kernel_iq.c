
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

#define MAX_PD 4  // MAX_PD is the max number of polydisperse parameters

typedef struct {
    int32_t pd_par[MAX_PD];     // index of the nth polydispersity variable
    int32_t pd_length[MAX_PD];  // length of the nth polydispersity weight vector
    int32_t pd_offset[MAX_PD];  // offset of pd weights in the par & weight vector
    int32_t pd_stride[MAX_PD];  // stride to move to the next index at this level
    int32_t pd_isvol[MAX_PD];   // True if parameter is a volume weighting parameter
    int32_t par_offset[NPARS];  // offset of par values in the par & weight vector
    int32_t par_coord[NPARS];   // polydispersity coordination bitvector
    int32_t fast_coord_index[NPARS]; // index of the fast coordination parameters
    int32_t fast_coord_count;   // number of parameters coordinated with pd 1
    int32_t theta_var;          // id of spherical correction variable
    int32_t fast_theta;         // true if spherical correction depends on pd 1
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
    global const double *pars,
    global const double *q, // nq q values, with padding to boundary
    global double *result,  // nq+3 return values, again with padding
    const double cutoff     // cutoff in the polydispersity weight product
    )
{
  // Storage for the current parameter values.  These will be updated as we
  // walk the polydispersity cube.
  local ParameterBlock local_pars;  // current parameter values
  double *pvec = (double *)(&local_pars);  // Alias named parameters with a vector

  local int offset[NPARS-2];

#if 1 // defined(USE_SHORTCUT_OPTIMIZATION)
  if (problem->pd_length[0] == 1) {
    // Shouldn't need to copy!!

    for (int k=0; k < NPARS; k++) {
      pvec[k] = pars[k+2];  // skip scale and background
    }
    printf("calculating\n");

    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i < nq; i++) {
      const double scattering = CALL_IQ(q, i, local_pars);
      result[i] += pars[0]*scattering + pars[1];
    }
    printf("returning\n");
    return;
  }
  printf("falling through\n");
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
    norm_vol = result[nq+2];
  }

  // Location in the polydispersity hypercube, one index per dimension.
  local int pd_index[MAX_PD];

  // Trigger the reset behaviour that happens at the end the fast loop
  // by setting the initial index >= weight vector length.
  pd_index[0] = problem->pd_length[0];


  // need product of weights at every Iq calc, so keep product of
  // weights from the outer loops so that weight = partial_weight * fast_weight
  double partial_weight = NAN; // product of weight w4*w3*w2 but not w1
  double partial_volweight = NAN;
  double weight = 1.0;        // set to 1 in case there are no weights
  double vol_weight = 1.0;    // set to 1 in case there are no vol weights
  double spherical_correction = 1.0;  // correction for latitude variation

  // Loop over the weights then loop over q, accumulating values
  for (int loop_index=pd_start; loop_index < pd_stop; loop_index++) {
    // check if indices need to be updated
    if (pd_index[0] >= problem->pd_length[0]) {

      // RESET INDICES
      pd_index[0] = loop_index%problem->pd_length[0];
      partial_weight = 1.0;
      partial_volweight = 1.0;
      for (int k=1; k < MAX_PD; k++) {
        pd_index[k] = (loop_index%problem->pd_length[k])/problem->pd_stride[k];
        const double wi = weights[problem->pd_offset[0]+pd_index[0]];
        partial_weight *= wi;
        if (problem->pd_isvol[k]) partial_volweight *= wi;
      }
      for (int k=0; k < NPARS; k++) {
        int coord = problem->par_coord[k];
        int this_offset = problem->par_offset[k];
        int block_size = 1;
        for (int bit=0; bit < MAX_PD && coord != 0; bit++) {
          if (coord&1) {
              this_offset += block_size * pd_index[bit];
              block_size *= problem->pd_length[bit];
          }
          coord /= 2;
        }
        offset[k] = this_offset;
        pvec[k] = pars[this_offset];
      }
      weight = partial_weight * weights[problem->pd_offset[0]+pd_index[0]];
      if (problem->theta_var >= 0) {
        spherical_correction = fabs(cos(M_PI_180*pvec[problem->theta_var]));
      }
      if (!problem->fast_theta) {
        weight *= spherical_correction;
      }

    } else {

      // INCREMENT INDICES
      pd_index[0] += 1;
      const double wi = weights[problem->pd_offset[0]+pd_index[0]];
      weight = partial_weight*wi;
      if (problem->pd_isvol[0]) vol_weight *= wi;
      for (int k=0; k < problem->fast_coord_count; k++) {
        pvec[problem->fast_coord_index[k]]
            = pars[offset[problem->fast_coord_index[k]]++];
      }
      if (problem->fast_theta) {
        weight *= fabs(cos(M_PI_180*pvec[problem->theta_var]));
      }
    }

    #ifdef INVALID
    if (INVALID(local_pars)) continue;
    #endif

    // Accumulate I(q)
    // Note: weight==0 must always be excluded
    if (weight > cutoff) {
      norm += weight;
      vol += vol_weight * CALL_VOLUME(local_pars);
      norm_vol += vol_weight;

      #ifdef USE_OPENMP
      #pragma omp parallel for
      #endif
      for (int i=0; i < nq; i++) {
        const double scattering = CALL_IQ(q, i, local_pars);
        result[i] += weight*scattering;
      }
    }
  }
  //Makes a normalization avialable for the next round
  result[nq] = norm;
  result[nq+1] = vol;
  result[nq+2] = norm_vol;

  //End of the PD loop we can normalize
  if (pd_stop >= problem->pd_stride[MAX_PD-1]) {
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
