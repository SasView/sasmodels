
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
    int fast_coord_index[NPARS]; // index of the fast coordination parameters
    int fast_coord_count;   // number of parameters coordinated with pd 1
    int theta_var;          // id of spherical correction variable
    int fast_theta;         // true if spherical correction depends on pd 1
} ProblemDetails;

typedef struct {
    PARAMETER_DECL;
} ParameterBlock;

#define FULL_KERNEL_NAME KERNEL_NAME ## _ ## IQ_FUNC
KERNEL
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
  double spherical_correction = 1.0;  // correction for latitude variation

  // Loop over the weights then loop over q, accumulating values
  for (int loop_index=pd_start; loop_index < pd_stop; loop_index++) {
    // check if indices need to be updated
    if (pd_index[0] >= pd_length[0]) {

      // RESET INDICES
      pd_index[0] = loop_index%pd_length[0];
      partial_weight = 1.0;
      partial_volweight = 1.0;
      for (int k=1; k < MAX_PD; k++) {
        pd_index[k] = (loop_index%pd_length[k])/pd_stride[k];
        const double wi = weights[pd_offset[0]+pd_index[0]];
        partial_weight *= wi;
        if (pd_isvol[k]) partial_volweight *= wi;
      }
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
      weight = partial_weight * weights[pd_offset[0]+pd_index[0]]
      if (theta_var >= 0) {
        spherical_correction = fabs(cos(M_PI_180*parvec[theta_var]));
      }
      if (!fast_theta) weight *= spherical_correction;

    } else {

      // INCREMENT INDICES
      pd_index[0] += 1;
      const double wi = weights[pd_offset[0]+pd_index[0]];
      weight = partial_weight*wi;
      if (pd_isvol[0]) vol_weight *= wi;
      for (int k=0; k < problem->fast_coord_count; k++) {
        parvec[ fast_coord_index[k]]
            = pars[offset[fast_coord_index[k]] + pd_index[0]];
      }
      if (fast_theta) weight *= fabs(cos(M_PI_180*parvec[theta_var]));

    }

    #ifdef INVALID
    if (INVALID(local_pars)) continue;
    #endif

    if (weight > cutoff) {
      norm += weight;
      vol += vol_weight * CALL_VOLUME(local_pars);
      norm_vol += vol_weight;

      #ifdef USE_OPENMP
      #pragma omp parallel for
      #endif
      for (int i=0; i < nq; i++) {
      {
        const double scattering = CALL_IQ(q, nq, i, local_pars);
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
