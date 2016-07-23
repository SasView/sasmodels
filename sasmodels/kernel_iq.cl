
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

  // who we are and what element we are working with
  const int q_index = get_global_id(0);
  if (q_index >= nq) return;

  // Storage for the current parameter values.  These will be updated as we
  // walk the polydispersity cube.  local_values will be aliased to pvec.
  ParameterBlock local_values;
  double *pvec = (double *)&local_values;

  // Fill in the initial variables
  for (int i=0; i < NPARS; i++) {
    pvec[i] = values[2+i];
//if (q_index==0) printf("p%d = %g\n",i, pvec[i]);
  }

#ifdef MAGNETIC
  // Location of the sld parameters in the parameter pvec.
  // These parameters are updated with the effective sld due to magnetism.
  const int32_t slds[] = { MAGNETIC_PARS };

  const double up_frac_i = values[NPARS+2];
  const double up_frac_f = values[NPARS+3];
  const double up_angle = values[NPARS+4];
  #define MX(_k) (values[NPARS+5+3*_k])
  #define MY(_k) (values[NPARS+6+3*_k])
  #define MZ(_k) (values[NPARS+7+3*_k])

  // TODO: could precompute these outside of the kernel.
  // Interpret polarization cross section.
  double uu, dd, ud, du;
  double cos_mspin, sin_mspin;
  spins(up_frac_i, up_frac_f, &uu, &dd, &ud, &du);
  SINCOS(-up_angle*M_PI_180, sin_mspin, cos_mspin);
#endif // MAGNETIC

  double pd_norm, this_result;
  if (pd_start == 0) {
    pd_norm = this_result = 0.0;
  } else {
    pd_norm = result[nq];
    this_result = result[q_index];
  }
//if (q_index==0) printf("start %d %g %g\n", pd_start, pd_norm, this_result);

  global const double *pd_value = values + NUM_VALUES + 2;
  global const double *pd_weight = pd_value + details->pd_sum;

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
//if (q_index==0) printf("offset %d: %d %d\n", 3, details->pd_offset[3], NUM_VALUES);
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
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 4, p4, i4, n4, pvec[p4], weight4);
#elif MAX_PD>3
    const double weight4 = 1.0;
#endif
#if MAX_PD>3
  while (i3 < n3) {
    pvec[p3] = v3[i3];
    double weight3 = w3[i3] * weight4;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 3, p3, i3, n3, pvec[p3], weight3);
#elif MAX_PD>2
    const double weight3 = 1.0;
#endif
#if MAX_PD>2
  while (i2 < n2) {
    pvec[p2] = v2[i2];
    double weight2 = w2[i2] * weight3;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 2, p2, i2, n2, pvec[p2], weight2);
#elif MAX_PD>1
    const double weight2 = 1.0;
#endif
#if MAX_PD>1
  while (i1 < n1) {
    pvec[p1] = v1[i1];
    double weight1 = w1[i1] * weight2;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 1, p1, i1, n1, pvec[p1], weight1);
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
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 0, p0, i0, n0, pvec[p0], weight0);
    if (fast_theta) { // Theta is in inner loop
      spherical_correction = fmax(fabs(cos(M_PI_180*pvec[p0])), 1.e-6);
    }
#else
    const double weight0 = 1.0;
#endif

//if (q_index == 0) {printf("step:%d of %d, pars:",step,pd_stop); for (int i=0; i < NPARS; i++) printf("p%d=%g ",i, pvec[i]); printf("\n"); }
//if (q_index == 0) printf("sphcor: %g\n", spherical_correction);

    #ifdef INVALID
    if (!INVALID(local_values))
    #endif
    {
      // Accumulate I(q)
      // Note: weight==0 must always be excluded
      if (weight0 > cutoff) {
        // spherical correction has some nasty effects when theta is +90 or -90
        // where it becomes zero.
        const double weight = weight0 * spherical_correction;
        pd_norm += weight * CALL_VOLUME(local_values);

#ifdef MAGNETIC
          const double qx = q[2*q_index];
          const double qy = q[2*q_index+1];
          const double qsq = qx*qx + qy*qy;

          // Constant across orientation, polydispersity for given qx, qy
          double px, py, pz;
          if (qsq > 1.e-16) {
            px = (qy*cos_mspin + qx*sin_mspin)/qsq;
            py = (qy*sin_mspin - qx*cos_mspin)/qsq;
            pz = 1.0;
          } else {
            px = py = pz = 0.0;
          }

          double scattering = 0.0;
          if (uu > 1.e-8) {
            for (int sk=0; sk<NUM_MAGNETIC; sk++) {
                const double perp = (qy*MX(sk) - qx*MY(sk));
                pvec[slds[sk]] = (values[slds[sk]+2] - perp*px)*uu;
            }
            scattering += CALL_IQ(q, q_index, local_values);
          }

          if (dd > 1.e-8){
            for (int sk=0; sk<NUM_MAGNETIC; sk++) {
                const double perp = (qy*MX(sk) - qx*MY(sk));
                pvec[slds[sk]] = (values[slds[sk]+2] + perp*px)*dd;
            }
            scattering += CALL_IQ(q, q_index, local_values);
          }
          if (ud > 1.e-8){
            for (int sk=0; sk<NUM_MAGNETIC; sk++) {
                const double perp = (qy*MX(sk) - qx*MY(sk));
                pvec[slds[sk]] = perp*py*ud;
            }
            scattering += CALL_IQ(q, q_index, local_values);
            for (int sk=0; sk<NUM_MAGNETIC; sk++) {
                pvec[slds[sk]] = MZ(sk)*pz*ud;
            }
            scattering += CALL_IQ(q, q_index, local_values);
          }
          if (du > 1.e-8) {
            for (int sk=0; sk<NUM_MAGNETIC; sk++) {
                const double perp = (qy*MX(sk) - qx*MY(sk));
                pvec[slds[sk]] = perp*py*du;
            }
            scattering += CALL_IQ(q, q_index, local_values);
            for (int sk=0; sk<NUM_MAGNETIC; sk++) {
                pvec[slds[sk]] = -MZ(sk)*pz*du;
            }
            scattering += CALL_IQ(q, q_index, local_values);
          }
#else  // !MAGNETIC
        const double scattering = CALL_IQ(q, q_index, local_values);
#endif // !MAGNETIC
        this_result += weight * scattering;
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

//if (q_index==0) printf("res: %g/%g\n", this_result, pd_norm);
  // Remember the current result and the updated norm.
  result[q_index] = this_result;
  if (q_index == 0) result[nq] = pd_norm;
}
