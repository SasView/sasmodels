
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
#endif // _PAR_BLOCK_


#if defined(MAGNETIC) && NUM_MAGNETIC>0

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
static void set_spins(double in_spin, double out_spin, double spins[4])
{
  in_spin = clip(in_spin, 0.0, 1.0);
  out_spin = clip(out_spin, 0.0, 1.0);
  spins[0] = sqrt(sqrt((1.0-in_spin) * (1.0-out_spin))); // dd
  spins[1] = sqrt(sqrt((1.0-in_spin) * out_spin));       // du
  spins[2] = sqrt(sqrt(in_spin * (1.0-out_spin)));       // ud
  spins[3] = sqrt(sqrt(in_spin * out_spin));             // uu
}

static double mag_sld(double qx, double qy, double p,
                       double mx, double my, double sld)
{
    const double perp = qy*mx - qx*my;
    return sld + perp*p;
}

#endif // MAGNETIC

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
  ParameterBlock local_values;
  double *pvec = (double *)&local_values;

#if defined(MAGNETIC) && NUM_MAGNETIC>0
  // Location of the sld parameters in the parameter pvec.
  // These parameters are updated with the effective sld due to magnetism.
  #if NUM_MAGNETIC > 3
  const int32_t slds[] = { MAGNETIC_PARS };
  #endif

  // TODO: could precompute these outside of the kernel.
  // Interpret polarization cross section.
  //     up_frac_i = values[NUM_PARS+2];
  //     up_frac_f = values[NUM_PARS+3];
  //     up_angle = values[NUM_PARS+4];
  double spins[4];
  double cos_mspin, sin_mspin;
  set_spins(values[NUM_PARS+2], values[NUM_PARS+3], spins);
  SINCOS(-values[NUM_PARS+4]*M_PI_180, sin_mspin, cos_mspin);
#endif // MAGNETIC

  // Fill in the initial variables
  //   values[0] is scale
  //   values[1] is background
  #ifdef USE_OPENMP
  #pragma omp parallel for
  #endif
  for (int i=0; i < NUM_PARS; i++) {
    pvec[i] = values[2+i];
//printf("p%d = %g\n",i, pvec[i]);
  }

  double pd_norm;
//printf("start: %d %d\n",pd_start, pd_stop);
  if (pd_start == 0) {
    pd_norm = 0.0;
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int q_index=0; q_index < nq; q_index++) result[q_index] = 0.0;
//printf("initializing %d\n", nq);
  } else {
    pd_norm = result[nq];
  }
//printf("start %d %g %g\n", pd_start, pd_norm, result[0]);

#if MAX_PD>0
  global const double *pd_value = values + NUM_VALUES + 2;
  global const double *pd_weight = pd_value + details->pd_sum;
#endif

  // Jump into the middle of the polydispersity loop
#if MAX_PD>4
  int n4=details->pd_length[4];
  int i4=(pd_start/details->pd_stride[4])%n4;
  const int p4=details->pd_par[4];
  global const double *v4 = pd_value + details->pd_offset[4];
  global const double *w4 = pd_weight + details->pd_offset[4];
#endif
#if MAX_PD>3
  int n3=details->pd_length[3];
  int i3=(pd_start/details->pd_stride[3])%n3;
  const int p3=details->pd_par[3];
  global const double *v3 = pd_value + details->pd_offset[3];
  global const double *w3 = pd_weight + details->pd_offset[3];
//printf("offset %d: %d %d\n", 3, details->pd_offset[3], NUM_VALUES);
#endif
#if MAX_PD>2
  int n2=details->pd_length[2];
  int i2=(pd_start/details->pd_stride[2])%n2;
  const int p2=details->pd_par[2];
  global const double *v2 = pd_value + details->pd_offset[2];
  global const double *w2 = pd_weight + details->pd_offset[2];
#endif
#if MAX_PD>1
  int n1=details->pd_length[1];
  int i1=(pd_start/details->pd_stride[1])%n1;
  const int p1=details->pd_par[1];
  global const double *v1 = pd_value + details->pd_offset[1];
  global const double *w1 = pd_weight + details->pd_offset[1];
#endif
#if MAX_PD>0
  int n0=details->pd_length[0];
  int i0=(pd_start/details->pd_stride[0])%n0;
  const int p0=details->pd_par[0];
  global const double *v0 = pd_value + details->pd_offset[0];
  global const double *w0 = pd_weight + details->pd_offset[0];
//printf("w0:%p, values:%p, diff:%d, %d\n",w0,values,(w0-values),NUM_VALUES);
#endif


  double spherical_correction=1.0;
  const int theta_par = details->theta_par;
#if MAX_PD>0
  const int fast_theta = (theta_par == p0);
  const int slow_theta = (theta_par >= 0 && !fast_theta);
#else
  const int slow_theta = (theta_par >= 0);
#endif

  int step = pd_start;


#if MAX_PD>4
  const double weight5 = 1.0;
  while (i4 < n4) {
    pvec[p4] = v4[i4];
    double weight4 = w4[i4] * weight5;
//printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 4, p4, i4, n4, pvec[p4], weight4);
#elif MAX_PD>3
    const double weight4 = 1.0;
#endif
#if MAX_PD>3
  while (i3 < n3) {
    pvec[p3] = v3[i3];
    double weight3 = w3[i3] * weight4;
//printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 3, p3, i3, n3, pvec[p3], weight3);
#elif MAX_PD>2
    const double weight3 = 1.0;
#endif
#if MAX_PD>2
  while (i2 < n2) {
    pvec[p2] = v2[i2];
    double weight2 = w2[i2] * weight3;
//printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 2, p2, i2, n2, pvec[p2], weight2);
#elif MAX_PD>1
    const double weight2 = 1.0;
#endif
#if MAX_PD>1
  while (i1 < n1) {
    pvec[p1] = v1[i1];
    double weight1 = w1[i1] * weight2;
//printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 1, p1, i1, n1, pvec[p1], weight1);
#elif MAX_PD>0
    const double weight1 = 1.0;
#endif
    if (slow_theta) { // Theta is not in inner loop
      spherical_correction = fmax(fabs(cos(M_PI_180*pvec[theta_par])), 1.e-6);
    }
#if MAX_PD>0
  while(i0 < n0) {
    pvec[p0] = v0[i0];
    double weight0 = w0[i0] * weight1;
//printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 0, p0, i0, n0, pvec[p0], weight0);
    if (fast_theta) { // Theta is in inner loop
      spherical_correction = fmax(fabs(cos(M_PI_180*pvec[p0])), 1.e-6);
    }
#else
    const double weight0 = 1.0;
#endif

//printf("step:%d of %d, pars:",step,pd_stop); for (int i=0; i < NUM_PARS; i++) printf("p%d=%g ",i, pvec[i]); printf("\n");
//printf("sphcor: %g\n", spherical_correction);

    #ifdef INVALID
    if (!INVALID(local_values))
    #endif
    {
      // Accumulate I(q)
      // Note: weight==0 must always be excluded
      if (weight0 > cutoff) {
        // spherical correction is set at a minimum of 1e-6, otherwise there
        // would be problems looking at models with theta=90.
        const double weight = weight0 * spherical_correction;
        pd_norm += weight * CALL_VOLUME(local_values);

        #ifdef USE_OPENMP
        #pragma omp parallel for
        #endif
        for (int q_index=0; q_index<nq; q_index++) {
#if defined(MAGNETIC) && NUM_MAGNETIC > 0
          const double qx = q[2*q_index];
          const double qy = q[2*q_index+1];
          const double qsq = qx*qx + qy*qy;

          // Constant across orientation, polydispersity for given qx, qy
          double scattering = 0.0;
          // TODO: what is the magnetic scattering at q=0
          if (qsq > 1.e-16) {
            double p[4];
            p[0] = (qy*cos_mspin + qx*sin_mspin)/qsq;
            p[3] = -p[0];
            p[1] = p[2] = (qy*sin_mspin - qx*cos_mspin)/qsq;

            for (int index=0; index<4; index++) {
              const double xs = spins[index];
              if (xs > 1.e-8) {
                const int spin_flip = (index==1) || (index==2);
                const double pk = p[index];
                for (int axis=0; axis<=spin_flip; axis++) {
                  #define M1 NUM_PARS+5
                  #define M2 NUM_PARS+8
                  #define M3 NUM_PARS+13
                  #define SLD(_M_offset, _sld_offset) \
                      pvec[_sld_offset] = xs * (axis \
                      ? (index==1 ? -values[_M_offset+2] : values[_M_offset+2]) \
                      : mag_sld(qx, qy, pk, values[_M_offset], values[_M_offset+1], \
                                (spin_flip ? 0.0 : values[_sld_offset+2])))
                  #if NUM_MAGNETIC==1
                      SLD(M1, MAGNETIC_PAR1);
                  #elif NUM_MAGNETIC==2
                      SLD(M1, MAGNETIC_PAR1);
                      SLD(M2, MAGNETIC_PAR2);
                  #elif NUM_MAGNETIC==3
                      SLD(M1, MAGNETIC_PAR1);
                      SLD(M2, MAGNETIC_PAR2);
                      SLD(M3, MAGNETIC_PAR3);
                  #else
                  for (int sk=0; sk<NUM_MAGNETIC; sk++) {
                      SLD(M1+3*sk, slds[sk]);
                  }
                  #endif
                  scattering += CALL_IQ(q, q_index, local_values);
                }
              }
            }
          }
#else  // !MAGNETIC
          const double scattering = CALL_IQ(q, q_index, local_values);
#endif // !MAGNETIC
//printf("q_index:%d %g %g %g %g\n",q_index, scattering, weight, spherical_correction, weight0);
          result[q_index] += weight * scattering;
        }
      }
    }
    ++step;
#if MAX_PD>0
    if (step >= pd_stop) break;
    ++i0;
  }
  i0 = 0;
#endif
#if MAX_PD>1
    if (step >= pd_stop) break;
    ++i1;
  }
  i1 = 0;
#endif
#if MAX_PD>2
    if (step >= pd_stop) break;
    ++i2;
  }
  i2 = 0;
#endif
#if MAX_PD>3
    if (step >= pd_stop) break;
    ++i3;
  }
  i3 = 0;
#endif
#if MAX_PD>4
    if (step >= pd_stop) break;
    ++i4;
  }
  i4 = 0;
#endif

//printf("res: %g/%g\n", result[0], pd_norm);
  // Remember the updated norm.
  result[nq] = pd_norm;
}
