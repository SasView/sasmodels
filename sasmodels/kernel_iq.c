
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
    int32_t pd_isvol[MAX_PD];   // True if parameter is a volume weighting parameter
#endif // MAX_PD > 0
    int32_t par_offset[NPARS];  // offset of par values in the value & weight vector
    int32_t par_coord[NPARS];   // polydispersity coordination bitvector
    int32_t fast_coord_pars[NPARS]; // ids of the fast coordination parameters
    int32_t fast_coord_count;   // number of parameters coordinated with pd 1
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

  // Monodisperse computation
  if (pd_stop == 1) {
    // Shouldn't need to copy!!
    for (int k=0; k < NPARS; k++) {
      pvec[k] = values[k+2];  // skip scale and background
    }

    const double volume = CALL_VOLUME(local_values);
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i < nq; i++) {
      double scattering = CALL_IQ(q, i, local_values);
      if (volume != 0.0) scattering /= volume;
      result[i] = values[0]*scattering + values[1];
    }
    return;
  }

#if MAX_PD > 0
  //printf("Entering polydispersity\n");

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

  // polydispersity loop index positions
  local int offset[NPARS];  // NPARS excludes scale/background

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
        const double wi = weights[problem->pd_offset[k]+pd_index[k]];
        partial_weight *= wi;
        if (problem->pd_isvol[k]) partial_volweight *= wi;
      }
      //printf("slow %d: ", loop_index);
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
        pvec[k] = values[this_offset];
        //printf("p[%d]=v[%d]=%g ", k, offset[k], pvec[k]);
      }
      //printf("\n");
      weight = partial_weight * weights[problem->pd_offset[0]+pd_index[0]];
      if (problem->theta_par >= 0) {
        spherical_correction = fabs(cos(M_PI_180*pvec[problem->theta_par]));
      }
      if (problem->theta_par == problem->pd_par[0]) {
        weight *= spherical_correction;
      }
      pd_index[0] += 1;

    } else {

      // INCREMENT INDICES
      const double wi = weights[problem->pd_offset[0]+pd_index[0]];
      weight = partial_weight*wi;
      if (problem->pd_isvol[0]) vol_weight *= wi;
      //printf("fast %d: ", loop_index);
      for (int k=0; k < problem->fast_coord_count; k++) {
        const int pindex = problem->fast_coord_pars[k];
        pvec[pindex] = values[++offset[pindex]];
        //printf("p[%d]=v[%d]=%g ", pindex, offset[pindex], pvec[pindex]);
      }
      //printf("\n");
      if (problem->theta_par == problem->pd_par[0]) {
        weight *= fabs(cos(M_PI_180*pvec[problem->theta_par]));
      }
      pd_index[0] += 1;
    }
    #ifdef INVALID
    if (INVALID(local_values)) continue;
    #endif

    // Accumulate I(q)
    // Note: weight==0 must always be excluded
    if (weight > cutoff) {
      norm += weight;
      vol += vol_weight * CALL_VOLUME(local_values);
      norm_vol += vol_weight;

      #ifdef USE_OPENMP
      #pragma omp parallel for
      #endif
      for (int i=0; i < nq; i++) {
        const double scattering = CALL_IQ(q, i, local_values);
        result[i] += weight*scattering;
      }
    }
  }

  // Make normalization available for the next round
  result[nq] = norm;
  result[nq+1] = vol;
  result[nq+2] = norm_vol;

  // End of the PD loop we can normalize
  if (pd_stop >= problem->pd_stride[MAX_PD-1]) {
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i < nq; i++) {
      if (vol*norm_vol != 0.0) {
        result[i] *= norm_vol/vol;
      }
      result[i] = values[0]*result[i]/norm + values[1];
    }
  }
#endif // MAX_PD > 0
}
