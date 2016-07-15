
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
    int32_t pd_prod;            // total number of voxels in hypercube
    int32_t pd_sum;             // total length of the weights vector
    int32_t num_active;         // number of non-trivial pd loops
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
    global const double *values,
    global const double *q, // nq q values, with padding to boundary
    global double *result,  // nq+1 return values, again with padding
    const double cutoff     // cutoff in the polydispersity weight product
    )
{
  // Storage for the current parameter values.  These will be updated as we
  // walk the polydispersity cube.  local_values will be aliased to pvec.
  local ParameterBlock local_values;

  // who we are and what element we are working with
  const int q_index = get_global_id(0);
  const int thread = get_local_id(0);

  // Fill in the initial variables
  event_t e = async_work_group_copy((local double *)&local_values, values+2, NPARS, 0);
  wait_group_events(1, &e);

  // Monodisperse computation
  if (details->num_active == 0) {
    double norm, scale, background;
    // TODO: only needs to be done by one process...
    #ifdef INVALID
    if (INVALID(local_values)) { return; }
    #endif

    norm = CALL_VOLUME(local_values);
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

  // "values" is global and can't be assigned to a local, so even though only
  // the alias is only needed for thread 0 it is allocated in all threads.
  global const double *pd_value = values+2+NPARS;
  global const double *pd_weight = pd_value+details->pd_sum;

  // need product of weights at every Iq calc, so keep product of
  // weights from the outer loops so that weight = partial_weight * fast_weight
  local double pd_norm;
  local double partial_weight; // product of weight w4*w3*w2 but not w1
  local double spherical_correction; // cosine correction for latitude variation
  local double weight; // product of partial_weight*w1*spherical_correction
  local double *pvec;
  local int p0_par;
  local int p0_length;
  local int p0_offset;
  local int p0_is_theta;
  local int p0_index;

  // Number of elements in the longest polydispersity loop
  barrier(CLK_LOCAL_MEM_FENCE);
  if (thread == 0) {
    pvec = (local double *)(&local_values);

    // Number of elements in the longest polydispersity loop
    p0_par = details->pd_par[0];
    p0_length = details->pd_length[0];
    p0_offset = details->pd_offset[0];
    p0_is_theta = (p0_par == details->theta_par);

    // Trigger the reset behaviour that happens at the end the fast loop
    // by setting the initial index >= weight vector length.
    p0_index = p0_length;

    // Default the spherical correction to 1.0 in case it is not otherwise set
    spherical_correction = 1.0;

    // Since we are no longer looping over the entire polydispersity hypercube
    // for each q, we need to track the result and normalization values between
    // calls.  This means initializing them to 0 at the start and accumulating
    // them between calls.
    pd_norm = pd_start == 0 ? 0.0 : result[nq];
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  if (q_index < nq) {
    this_result = pd_start == 0 ? 0.0 : result[q_index];
  }

  // Loop over the weights then loop over q, accumulating values
  for (int loop_index=pd_start; loop_index < pd_stop; loop_index++) {
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread == 0) {
      // check if fast loop needs to be reset
      if (p0_index == p0_length) {
        //printf("should be here with %d active\n", num_active);

        // Compute position in polydispersity hypercube and partial weight
        partial_weight = 1.0;
        for (int k=1; k < details->num_active; k++) {
          int pk = details->pd_par[k];
          int index = details->pd_offset[k] + (loop_index/details->pd_stride[k])%details->pd_length[k];
          pvec[pk] = pd_value[index];
          partial_weight *= pd_weight[index];
          //printf("index[%d] = %d\n",k,index);
          if (pk == details->theta_par) {
            spherical_correction = fmax(fabs(cos(M_PI_180*pvec[pk])), 1.e-6);
          }
        }
        p0_index = loop_index%p0_length;
        //printf("\n");
      }

      // Update parameter p0
      weight = partial_weight*pd_weight[p0_offset + p0_index];
      pvec[p0_par] = pd_value[p0_offset + p0_index];
      if (p0_is_theta) {
        spherical_correction = fmax(fabs(cos(M_PI_180*pvec[p0_par])), 1.e-6);
      }
      p0_index++;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    //printf("\n");

    // Increment fast index

    #ifdef INVALID
    if (INVALID(local_values)) continue;
    #endif

    // Accumulate I(q)
    // Note: weight==0 must always be excluded
    if (weight > cutoff) {
      // spherical correction has some nasty effects when theta is +90 or -90
      // where it becomes zero.  If the entirety of the correction
      weight *= spherical_correction;
      pd_norm += weight * CALL_VOLUME(local_values);

      const double scattering = CALL_IQ(q, q_index, local_values);
      this_result += weight*scattering;
    }
  }

  if (q_index < nq) {
    if (pd_stop >= details->pd_prod) {
      // End of the PD loop we can normalize
      double scale, background;
      scale = values[0];
      background = values[1];
      result[q_index] = (pd_norm>0. ? scale*this_result/pd_norm + background : background);
    } else {
      // Partial result, so remember it but don't normalize it.
      result[q_index] = this_result;
    }

    // Remember the updated norm.
    if (q_index == 0) result[nq] = pd_norm;
  }

#endif // MAX_PD > 0
}
