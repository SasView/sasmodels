
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
    PARAMETER_TABLE
} ParameterBlock;
#endif

#ifdef MAGNETIC
const int32_t magnetic[] = { MAGNETIC_PARS };
#endif

#ifdef MAGNETIC
// Return value restricted between low and high
static double clip(double value, double low, double high)
{
  return (value < low ? low : (value > high ? high : value));
}

// Compute spin cross sections given in_spin and out_spin
// To convert spin cross sections to sld b:
//     uu * (sld - m_sigma_x);
//     dd * (sld + m_sigma_x);
//     ud * (m_sigma_y + 1j*m_sigma_z);
//     du * (m_sigma_y - 1j*m_sigma_z);
static void spins(double in_spin, double out_spin,
    double *uu, double *dd, double *ud, double *du)
{
  in_spin = clip(in_spin, 0.0, 1.0);
  out_spin = clip(out_spin, 0.0, 1.0);
  *uu = sqrt(sqrt(in_spin * out_spin));
  *dd = sqrt(sqrt((1.0-in_spin) * (1.0-out_spin)));
  *ud = sqrt(sqrt(in_spin * (1.0-out_spin)));
  *du = sqrt(sqrt((1.0-in_spin) * out_spin));
}

// Convert polar to rectangular coordinates.
static void polrec(double r, double theta, double phi,
  double *x, double *y, double *z)
{
  double cos_theta, sin_theta, cos_phi, sin_phi;
  SINCOS(theta*M_PI_180, sin_theta, cos_theta);
  SINCOS(phi*M_PI_180, sin_phi, cos_phi);
  *x = r * cos_theta * cos_phi;
  *y = r * sin_theta;
  *z = -r * cos_theta * sin_phi;
}
#endif

kernel
void KERNEL_NAME(
    int32_t nq,                 // number of q values
    const int32_t pd_start,     // where we are in the polydispersity loop
    const int32_t pd_stop,      // where we are stopping in the polydispersity loop
    global const ProblemDetails *details,
   // global const  // TODO: make it const again!
    double *values,
    global const double *q, // nq q values, with padding to boundary
    global double *result,  // nq+3 return values, again with padding
    const double cutoff     // cutoff in the polydispersity weight product
    )
{
  // Storage for the current parameter values.  These will be updated as we
  // walk the polydispersity cube.
  ParameterBlock local_values;  // current parameter values
  double *pvec = (double *)(&local_values);  // Alias named parameters with a vector

  // Fill in the initial variables
  //   values[0] is scale
  //   values[1] is background
  #ifdef USE_OPENMP
  #pragma omp parallel for
  #endif
  for (int k=0; k < NPARS; k++) {
    pvec[k] = values[k+2];
  }
#ifdef MAGNETIC
  const double up_frac_i = values[NPARS+2];
  const double up_frac_f = values[NPARS+3];
  const double up_angle = values[NPARS+4];
  #define MX(_k) (values[NPARS+5+3*_k])
  #define MY(_k) (values[NPARS+6+3*_k])
  #define MZ(_k) (values[NPARS+7+3*_k])

  // TODO: precompute this on the python side
  // Convert polar to rectangular coordinates in place.
  if (pd_start == 0) {  // Update in place; only do this for the first hunk!
//printf("spin: %g %g %g\n", up_frac_i, up_frac_f, up_angle);
    for (int mag=0; mag < NUM_MAGNETIC; mag++) {
//printf("mag %d: %g %g %g\n", mag, MX(mag), MY(mag), MZ(mag));
        polrec(MX(mag), MY(mag), MZ(mag), &MX(mag), &MY(mag), &MZ(mag));
//printf("   ==>: %g %g %g\n", MX(mag), MY(mag), MZ(mag));
    }
  }
  // Interpret polarization cross section.
  double uu, dd, ud, du;
  double cos_mspin, sin_mspin;
  spins(up_frac_i, up_frac_f, &uu, &dd, &ud, &du);
  SINCOS(-up_angle*M_PI_180, sin_mspin, cos_mspin);
#endif

  // Monodisperse computation
  if (details->num_active == 0) {
    double norm, scale, background;
    #ifdef INVALID
    if (INVALID(local_values)) { return; }
    #endif

    norm = CALL_VOLUME(local_values);
    scale = values[0];
    background = values[1];

    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int q_index=0; q_index < nq; q_index++) {
#ifdef MAGNETIC
      const double qx = q[2*q_index];
      const double qy = q[2*q_index+1];
      const double qsq = qx*qx + qy*qy;

      // Constant across orientation, polydispersity for given qx, qy
      double px, py, pz;
      if (qsq > 1e-16) {
        px = (qy*cos_mspin + qx*sin_mspin)/qsq;
        py = (qy*sin_mspin - qx*cos_mspin)/qsq;
        pz = 1.0;
      } else {
        px = py = pz = 0.0;
      }

      double scattering = 0.0;
      if (uu > 1e-8) {
        for (int mag=0; mag<NUM_MAGNETIC; mag++) {
            const double perp = (qy*MX(mag) - qx*MY(mag));
            pvec[magnetic[mag]] = (values[magnetic[mag]+2] - perp*px)*uu;
        }
        scattering += CALL_IQ(q, q_index, local_values);
      }
      if (dd > 1e-8){
        for (int mag=0; mag<NUM_MAGNETIC; mag++) {
            const double perp = (qy*MX(mag) - qx*MY(mag));
            pvec[magnetic[mag]] = (values[magnetic[mag]+2] + perp*px)*dd;
        }
        scattering += CALL_IQ(q, q_index, local_values);
      }
      if (ud > 1e-8){
        for (int mag=0; mag<NUM_MAGNETIC; mag++) {
            const double perp = (qy*MX(mag) - qx*MY(mag));
            pvec[magnetic[mag]] = perp*py*ud;
        }
        scattering += CALL_IQ(q, q_index, local_values);
        for (int mag=0; mag<NUM_MAGNETIC; mag++) {
            pvec[magnetic[mag]] = MZ(mag)*pz*ud;
        }
        scattering += CALL_IQ(q, q_index, local_values);
      }
      if (du > 1e-8) {
        for (int mag=0; mag<NUM_MAGNETIC; mag++) {
            const double perp = (qy*MX(mag) - qx*MY(mag));
            pvec[magnetic[mag]] = perp*py*du;
        }
        scattering += CALL_IQ(q, q_index, local_values);
        for (int mag=0; mag<NUM_MAGNETIC; mag++) {
            pvec[magnetic[mag]] = -MZ(mag)*pz*du;
        }
        scattering += CALL_IQ(q, q_index, local_values);
      }
#else
      double scattering = CALL_IQ(q, q_index, local_values);
#endif
      result[q_index] = (norm>0. ? scale*scattering/norm + background : background);
    }
    return;
  }

#if MAX_PD > 0

#if MAGNETIC
  const double *pd_value = values+2+NPARS+3+3*NUM_MAGNETIC;
#else
  const double *pd_value = values+2+NPARS;
#endif
  const double *pd_weight = pd_value+details->pd_sum;

  // need product of weights at every Iq calc, so keep product of
  // weights from the outer loops so that weight = partial_weight * fast_weight
  double pd_norm;
  double partial_weight; // product of weight w4*w3*w2 but not w1
  double spherical_correction; // cosine correction for latitude variation
  double weight; // product of partial_weight*w1*spherical_correction

  // Number of elements in the longest polydispersity loop
  const int p0_par = details->pd_par[0];
  const int p0_length = details->pd_length[0];
  const int p0_offset = details->pd_offset[0];
  const int p0_is_theta = (p0_par == details->theta_par);
  int p0_index;

  // Trigger the reset behaviour that happens at the end the fast loop
  // by setting the initial index >= weight vector length.
  p0_index = p0_length;

  // Default the spherical correction to 1.0 in case it is not otherwise set
  spherical_correction = 1.0;

  // Since we are no longer looping over the entire polydispersity hypercube
  // for each q, we need to track the result and normalization values between
  // calls.  This means initializing them to 0 at the start and accumulating
  // them between calls.
  pd_norm = (pd_start == 0 ? 0.0 : result[nq]);

  if (pd_start == 0) {
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int q_index=0; q_index < nq; q_index++) {
      result[q_index] = 0.0;
    }
  }

  // Loop over the weights then loop over q, accumulating values
  for (int loop_index=pd_start; loop_index < pd_stop; loop_index++) {
    // check if fast loop needs to be reset
    if (p0_index == p0_length) {

      // Compute position in polydispersity hypercube and partial weight
      partial_weight = 1.0;
      for (int k=1; k < details->num_active; k++) {
        int pk = details->pd_par[k];
        int index = details->pd_offset[k] + (loop_index/details->pd_stride[k])%details->pd_length[k];
        pvec[pk] = pd_value[index];
        partial_weight *= pd_weight[index];
        if (pk == details->theta_par) {
          spherical_correction = fmax(fabs(cos(M_PI_180*pvec[pk])), 1.e-6);
        }
      }
      p0_index = loop_index%p0_length;
    }

    // Update parameter p0
    weight = partial_weight*pd_weight[p0_offset + p0_index];
    pvec[p0_par] = pd_value[p0_offset + p0_index];
    if (p0_is_theta) {
      spherical_correction = fmax(fabs(cos(M_PI_180*pvec[p0_par])), 1.e-6);
    }
    p0_index++;

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

      #ifdef USE_OPENMP
      #pragma omp parallel for
      #endif
      for (int q_index=0; q_index < nq; q_index++) {
#ifdef MAGNETIC
        const double qx = q[2*q_index];
        const double qy = q[2*q_index+1];
        const double qsq = qx*qx + qy*qy;

        // Constant across orientation, polydispersity for given qx, qy
        double px, py, pz;
        if (qsq > 1e-16) {
          px = (qy*cos_mspin + qx*sin_mspin)/qsq;
          py = (qy*sin_mspin - qx*cos_mspin)/qsq;
          pz = 1.0;
        } else {
          px = py = pz = 0.0;
        }

        double scattering = 0.0;
        if (uu > 1e-8) {
          for (int mag=0; mag<NUM_MAGNETIC; mag++) {
              const double perp = (qy*MX(mag) - qx*MY(mag));
              pvec[magnetic[mag]] = (values[magnetic[mag]+2] - perp*px)*uu;
          }
          scattering += CALL_IQ(q, q_index, local_values);
        }
        if (dd > 1e-8){
          for (int mag=0; mag<NUM_MAGNETIC; mag++) {
              const double perp = (qy*MX(mag) - qx*MY(mag));
              pvec[magnetic[mag]] = (values[magnetic[mag]+2] + perp*px)*dd;
          }
          scattering += CALL_IQ(q, q_index, local_values);
        }
        if (ud > 1e-8){
          for (int mag=0; mag<NUM_MAGNETIC; mag++) {
              const double perp = (qy*MX(mag) - qx*MY(mag));
              pvec[magnetic[mag]] = perp*py*ud;
          }
          scattering += CALL_IQ(q, q_index, local_values);
          for (int mag=0; mag<NUM_MAGNETIC; mag++) {
              pvec[magnetic[mag]] = MZ(mag)*pz*ud;
          }
          scattering += CALL_IQ(q, q_index, local_values);
        }
        if (du > 1e-8) {
          for (int mag=0; mag<NUM_MAGNETIC; mag++) {
              const double perp = (qy*MX(mag) - qx*MY(mag));
              pvec[magnetic[mag]] = perp*py*du;
          }
          scattering += CALL_IQ(q, q_index, local_values);
          for (int mag=0; mag<NUM_MAGNETIC; mag++) {
              pvec[magnetic[mag]] = -MZ(mag)*pz*du;
          }
          scattering += CALL_IQ(q, q_index, local_values);
        }
#else
        double scattering = CALL_IQ(q, q_index, local_values);
#endif
        result[q_index] += weight*scattering;
      }
    }
  }

  if (pd_stop >= details->pd_prod) {
    // End of the PD loop we can normalize
    double scale, background;
    scale = values[0];
    background = values[1];
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int q_index=0; q_index < nq; q_index++) {
      result[q_index] = (pd_norm>0. ? scale*result[q_index]/pd_norm + background : background);
    }
  }

  // Remember the updated norm.
  result[nq] = pd_norm;
#endif // MAX_PD > 0
}
