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

// NOTE: the following macros are defined in generate.py:
//
//  MAX_PD : the maximum number of dispersity loops allowed for this model,
//      which will be at most modelinfo.MAX_PD.
//  NUM_PARS : the number of parameters in the parameter table
//  NUM_VALUES : the number of values to skip at the start of the
//      values array before you get to the dispersity values.
//  PARAMETER_TABLE : list of parameter declarations used to create the
//      ParameterTable type.
//  KERNEL_NAME : model_Iq, model_Iqxy or model_Imagnetic.  This code is
//      included three times, once for each kernel type.
//  MAGNETIC : defined when the magnetic kernel is being instantiated
//  NUM_MAGNETIC : the number of magnetic parameters
//  MAGNETIC_PARS : a comma-separated list of indices to the sld
//      parameters in the parameter table.
//  TRANSLATION_VARS(table) : series of intermediate expressions used to
//      compute parameter substitions when reparameterizing a model.
//  VALID(table) : test if the current point is feesible to calculate.
//  CALL_VOLUME(form, shell, table) : assign form and shell values.
//  CALL_RADIUS_EFFECTIVE(mode, table) : call the R_eff function.
//  CALL_IQ(q, table) : call the Iq function for 1D calcs.
//  CALL_IQ_A(q, table) : call the Iq function with |q| for 2D data.
//  CALL_FQ(q, F1, F2, table) : call the Fq function for 1D calcs.
//  CALL_FQ_A(q, F1, F2, table) : call the Iq function with |q| for 2D data.
//  CALL_IQ_AC(qa, qc, table) : call the Iqxy function for symmetric shapes
//  CALL_IQ_ABC(qa, qc, table) : call the Iqxy function for asymmetric shapes
//  CALL_IQ_XY(qx, qy, table) : call the Iqxy function for arbitrary models
//  PROJECTION : equirectangular=1, sinusoidal=2
//      see explore/jitter.py for definitions.


#ifndef _PAR_BLOCK_ // protected block so we can include this code twice.
#define _PAR_BLOCK_

typedef struct {
#if MAX_PD > 0
    int32_t pd_par[MAX_PD];     // id of the nth dispersity variable
    int32_t pd_length[MAX_PD];  // length of the nth dispersity weight vector
    int32_t pd_offset[MAX_PD];  // offset of pd weights in the value & weight vector
    int32_t pd_stride[MAX_PD];  // stride to move to the next index at this level
#endif // MAX_PD > 0
    int32_t num_eval;           // total number of voxels in hypercube
    int32_t num_weights;        // total length of the weights vector
    int32_t num_active;         // number of non-trivial pd loops
    int32_t theta_par;          // id of first orientation variable
} ProblemDetails;

// Intel HD 4000 needs private arrays to be a multiple of 4 long
typedef struct {
    PARAMETER_TABLE
} ParameterTable;
typedef union {
    ParameterTable table;
    double vector[4*((NUM_PARS+3)/4)];
} ParameterBlock;
#endif // _PAR_BLOCK_

#if defined(MAGNETIC) && NUM_MAGNETIC > 0
// ===== Helper functions for magnetism =====

// vector algebra
void SET_VEC(double *vector, double v0, double v1, double v2)
{
    vector[0] = v0;
    vector[1] = v1;
    vector[2] = v2;
}

void SCALE_VEC(double *vector, double a)
{
    vector[0] = a*vector[0];
    vector[1] = a*vector[1];
    vector[2] = a*vector[2];
}

void ADD_VEC(double *result_vec, double *vec1, double *vec2)
{
    result_vec[0] = vec1[0] + vec2[0];
    result_vec[1] = vec1[1] + vec2[1];
    result_vec[2] = vec1[2] + vec2[2];
}

static double SCALAR_VEC( double *vec1, double *vec2)
{
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

static double MAG_VEC( double *vec)
{
    return sqrt(SCALAR_VEC(vec,vec));
}

void ORTH_VEC(double *result_vec, double *vec1, double *vec2)
{
    double scale =  SCALAR_VEC(vec1,vec2) / SCALAR_VEC(vec2,vec2);
    result_vec[0] = vec1[0] - scale * vec2[0];
    result_vec[1] = vec1[1] - scale * vec2[1];
    result_vec[2] = vec1[2] - scale * vec2[2];
}


// Return value restricted between low and high
static double clip(double value, double low, double high)
{
  return (value < low ? low : (value > high ? high : value));
}

// Compute spin cross sections given in_spin and out_spin
// To convert spin cross sections to sld b:
//     uu * (sld - m_perp_x);
//     dd * (sld + m_perp_x);
//     ud * (m_perp_y - 1j*m_perp_z);
//     du * (m_perp_y + 1j*m_perp_z);
//(x,y,z) is a local magnetic coordinate system. m_perp_x denotes the magnetic scattering vector along the polarisation and m_perp_y the component along the scattering vector, m_perpz is orthogonal to the others.
// weights for spin crosssections: dd du real, ud real, uu, du imag, ud imag
static void set_spin_weights(double in_spin, double out_spin, double weight[6])
{

  double norm;
  in_spin = clip(in_spin, 0.0, 1.0);
  out_spin = clip(out_spin, 0.0, 1.0);
  // Previous version of this function took the square root of the weights,
  // under the assumption that
  //
  //     w*I(q, rho1, rho2, ...) = I(q, sqrt(w)*rho1, sqrt(w)*rho2, ...)
  //
  // However, since the weights are applied to the final intensity and
  // are not interned inside the I(q) function, we want the full
  // weight and not the square root.  Anyway no function will ever use
  // set_spin_weights as part of calculating an amplitude, as the weights are
  // related to polarisation efficiency of the instrument. The weights serve to
  // construct various magnet scattering cross sections, which are linear combinations
  // of the spin-resolved cross sections. The polarisation efficiency e_in and e_out
  // are parameters ranging from 0.5 (unpolarised) beam to 1 (perfect optics).
  // For in_spin or out_spin <0.5 one assumes a CS, where the spin is reversed/flipped
  // with respect to the initial supermirror polariser. The actual polarisation efficiency
  // in this case is however e_in/out = 1-in/out_spin.

  if (out_spin < 0.5){norm=1-out_spin;}
  else{norm=out_spin;}


// The norm is needed to make sure that the scattering cross sections are
//correctly weighted, such that the sum of spin-resolved measurements adds up to
// the unpolarised or half-polarised scattering cross section. No intensity weighting
// needed on the incoming polariser side (assuming that a user), has normalised
// to the incoming flux with polariser in for SANSPOl and unpolarised beam, respectively.


  weight[0] = (1.0-in_spin) * (1.0-out_spin) / norm; // dd
  weight[1] = (1.0-in_spin) * out_spin / norm;       // du
  weight[2] = in_spin * (1.0-out_spin) / norm;       // ud
  weight[3] = in_spin * out_spin / norm;             // uu
  weight[4] = weight[1]; // du.imag
  weight[5] = weight[2]; // ud.imag
}

// Compute the magnetic sld
static double mag_sld(
  const unsigned int xs, // 0=dd, 1=du.real, 2=ud.real, 3=uu, 4=du.imag, 5=ud.imag
  const double qx, const double qy,
  const double cos_mspin, const double sin_mspin,
  const double cos_mphi, const double sin_mphi,
  const double sld,
  const double mx, const double my, const double mz
)
{
  double Mvector[3];
  double Pvector[3]; 
  double qvector[3];
  double rhom[3]; 
  double perpy[3];
  double perpz[3];  
  double Mperp[3];

  const double qsq = sqrt(qx*qx + qy*qy); 
  SET_VEC(qvector, qx / qsq, qy / qsq, 0);
  SET_VEC(Mvector, mx, my, mz);
  //const double px = sin_mspin * cos_mphi;
  //const double py = sin_mspin * sin_mphi;
  //const double pz = cos_mspin;
  SET_VEC(Pvector, sin_mspin * cos_mphi, sin_mspin * sin_mphi, cos_mspin);
  //two unit vectors spanning up the plane perpendicular to polarisation for SF scattering
  SET_VEC(perpy, -sin_mphi, cos_mphi, 0);
  SET_VEC(perpz, -cos_mspin * cos_mphi, -cos_mspin * sin_mphi, sin_mspin);
  ORTH_VEC(Mperp, Mvector, qvector);



  if (xs < 4) {

    switch (xs) {
      default: // keep compiler happy; condition ensures xs in [0,1,2,3]
      case 0: // dd => sld - D Pvector \cdot Mperp
          return sld - SCALAR_VEC(Pvector, Mperp);
      case 1: // du.real          
          return SCALAR_VEC(perpy, Mperp);
      case 2: // ud.real   
          return SCALAR_VEC(perpy, Mperp);
      case 3: // uu => sld + D Pvector \cdot Mperp
          return sld + SCALAR_VEC(Pvector, Mperp); 
    }
  } else {
    if (xs== 4) {
      return - SCALAR_VEC(perpz, Mperp);  // du.imag => - i  nperp \cdot MperpP
    } else { // index == 5 
      return + SCALAR_VEC(perpz, Mperp);   // du.imag => + i nperp \cdot MperpP
    }
  }
}




#endif

// ===== Helper functions for orientation and jitter =====

// To change the definition of the angles, run explore/angles.py, which
// uses sympy to generate the equations.

#if !defined(_QAC_SECTION) && defined(CALL_IQ_AC)
#define _QAC_SECTION

typedef struct {
    double R31, R32;
} QACRotation;

// Fill in the rotation matrix R from the view angles (theta, phi) and the
// jitter angles (dtheta, dphi).  This matrix can be applied to all of the
// (qx, qy) points in the image to produce R*[qx,qy]' = [qa,qc]'
static void
qac_rotation(
    QACRotation *rotation,
    double theta, double phi,
    double dtheta, double dphi)
{
    double sin_theta, cos_theta;
    double sin_phi, cos_phi;

    // reverse view matrix
    SINCOS(theta*M_PI_180, sin_theta, cos_theta);
    SINCOS(phi*M_PI_180, sin_phi, cos_phi);
    const double V11 = cos_phi*cos_theta;
    const double V12 = sin_phi*cos_theta;
    const double V21 = -sin_phi;
    const double V22 = cos_phi;
    const double V31 = sin_theta*cos_phi;
    const double V32 = sin_phi*sin_theta;

    // reverse jitter matrix
    SINCOS(dtheta*M_PI_180, sin_theta, cos_theta);
    SINCOS(dphi*M_PI_180, sin_phi, cos_phi);
    const double J31 = sin_theta;
    const double J32 = -sin_phi*cos_theta;
    const double J33 = cos_phi*cos_theta;

    // reverse matrix
    rotation->R31 = J31*V11 + J32*V21 + J33*V31;
    rotation->R32 = J31*V12 + J32*V22 + J33*V32;
}

// Apply the rotation matrix returned from qac_rotation to the point (qx,qy),
// returning R*[qx,qy]' = [qa,qc]'
static void
qac_apply(
    QACRotation *rotation,
    double qx, double qy,
    double *qab_out, double *qc_out)
{
    // Indirect calculation of qab, from qab^2 = |q|^2 - qc^2
    const double dqc = rotation->R31*qx + rotation->R32*qy;
    const double dqab_sq = -dqc*dqc + qx*qx + qy*qy;
    //*qab_out = sqrt(fabs(dqab_sq));
    *qab_out = dqab_sq > 0.0 ? sqrt(dqab_sq) : 0.0;
    *qc_out = dqc;
}
#endif // _QAC_SECTION

#if !defined(_QABC_SECTION) && defined(CALL_IQ_ABC)
#define _QABC_SECTION

typedef struct {
    double R11, R12;
    double R21, R22;
    double R31, R32;
} QABCRotation;

// Fill in the rotation matrix R from the view angles (theta, phi, psi) and the
// jitter angles (dtheta, dphi, dpsi).  This matrix can be applied to all of the
// (qx, qy) points in the image to produce R*[qx,qy]' = [qa,qb,qc]'
static void
qabc_rotation(
    QABCRotation *rotation,
    double theta, double phi, double psi,
    double dtheta, double dphi, double dpsi)
{
    double sin_theta, cos_theta;
    double sin_phi, cos_phi;
    double sin_psi, cos_psi;

    // reverse view matrix
    SINCOS(theta*M_PI_180, sin_theta, cos_theta);
    SINCOS(phi*M_PI_180, sin_phi, cos_phi);
    SINCOS(psi*M_PI_180, sin_psi, cos_psi);
    const double V11 = -sin_phi*sin_psi + cos_phi*cos_psi*cos_theta;
    const double V12 = sin_phi*cos_psi*cos_theta + sin_psi*cos_phi;
    const double V21 = -sin_phi*cos_psi - sin_psi*cos_phi*cos_theta;
    const double V22 = -sin_phi*sin_psi*cos_theta + cos_phi*cos_psi;
    const double V31 = sin_theta*cos_phi;
    const double V32 = sin_phi*sin_theta;

    // reverse jitter matrix
    SINCOS(dtheta*M_PI_180, sin_theta, cos_theta);
    SINCOS(dphi*M_PI_180, sin_phi, cos_phi);
    SINCOS(dpsi*M_PI_180, sin_psi, cos_psi);
    const double J11 = cos_psi*cos_theta;
    const double J12 = sin_phi*sin_theta*cos_psi + sin_psi*cos_phi;
    const double J13 = sin_phi*sin_psi - sin_theta*cos_phi*cos_psi;
    const double J21 = -sin_psi*cos_theta;
    const double J22 = -sin_phi*sin_psi*sin_theta + cos_phi*cos_psi;
    const double J23 = sin_phi*cos_psi + sin_psi*sin_theta*cos_phi;
    const double J31 = sin_theta;
    const double J32 = -sin_phi*cos_theta;
    const double J33 = cos_phi*cos_theta;

    // reverse matrix
    rotation->R11 = J11*V11 + J12*V21 + J13*V31;
    rotation->R12 = J11*V12 + J12*V22 + J13*V32;
    rotation->R21 = J21*V11 + J22*V21 + J23*V31;
    rotation->R22 = J21*V12 + J22*V22 + J23*V32;
    rotation->R31 = J31*V11 + J32*V21 + J33*V31;
    rotation->R32 = J31*V12 + J32*V22 + J33*V32;
}

// Apply the rotation matrix returned from qabc_rotation to the point (qx,qy),
// returning R*[qx,qy]' = [qa,qb,qc]'
static void
qabc_apply(
    QABCRotation *rotation,
    double qx, double qy,
    double *qa_out, double *qb_out, double *qc_out)
{
    *qa_out = rotation->R11*qx + rotation->R12*qy;
    *qb_out = rotation->R21*qx + rotation->R22*qy;
    *qc_out = rotation->R31*qx + rotation->R32*qy;
}

#endif // _QABC_SECTION

// ==================== KERNEL CODE ========================
kernel
void KERNEL_NAME(
    int32_t nq,                   // number of q values
    const int32_t pd_start,       // where we are in the dispersity loop
    const int32_t pd_stop,        // where we are stopping in the dispersity loop
    pglobal const ProblemDetails *details,
    pglobal const double *values, // parameter values and distributions
    pglobal const double *q,      // nq q values, with padding to boundary
    pglobal double *result,       // nq+1 return values, again with padding
    const double cutoff,          // cutoff in the dispersity weight product
    int32_t radius_effective_mode // which effective radius to compute
    )
{
#if defined(USE_GPU)
  // who we are and what element we are working with
  #if defined(USE_OPENCL)
  const int q_index = get_global_id(0);
  #else // USE_CUDA
  const int q_index = threadIdx.x + blockIdx.x * blockDim.x;
  #endif
  if (q_index >= nq) return;
#else
  // Define q_index here so that debugging statements can be written to work
  // for both OpenCL and DLL using:
  //    if (q_index == 0) {printf(...);}
  int q_index = 0;
#endif

  // ** Fill in the local values table **
  // Storage for the current parameter values.
  // These will be updated as we walk the dispersity mesh.
  ParameterBlock local_values;
  //   values[0] is scale
  //   values[1] is background
  #ifdef USE_OPENMP
  #pragma omp parallel for
  #endif
  for (int i=0; i < NUM_PARS; i++) {
    local_values.vector[i] = values[2+i];
    //if (q_index==0) printf("p%d = %g\n",i, local_values.vector[i]);
  }
  //if (q_index==0) printf("NUM_VALUES:%d  NUM_PARS:%d  MAX_PD:%d\n", NUM_VALUES, NUM_PARS, MAX_PD);
  //if (q_index==0) printf("start:%d  stop:%d\n", pd_start, pd_stop);

  // ** Precompute magnatism values **
#if defined(MAGNETIC) && NUM_MAGNETIC>0
  // Location of the sld parameters in the parameter vector.
  // These parameters are updated with the effective sld due to magnetism.
  const int32_t slds[] = { MAGNETIC_PARS };

  // Interpret polarization cross section.
  //     up_frac_i = values[NUM_PARS+2];
  //     up_frac_f = values[NUM_PARS+3];
  //     up_angle = values[NUM_PARS+4];
  //     up_phi = values[NUM_PARS+5];
  // TODO: could precompute more magnetism parameters before calling the kernel.
  double xs_weights[8];  // uu, ud real, du real, dd, ud imag, du imag, fill, fill
  double cos_mspin, sin_mspin;
  double cos_mphi, sin_mphi;
  set_spin_weights(values[NUM_PARS+2], values[NUM_PARS+3], xs_weights);
  SINCOS(values[NUM_PARS+4]*M_PI_180, sin_mspin, cos_mspin);
  SINCOS(values[NUM_PARS+5]*M_PI_180, sin_mphi, cos_mphi);
#endif // MAGNETIC

  // ** Fill in the initial results **
  // If pd_start is zero that means that we are starting a new calculation,
  // and must initialize the result to zero.  Otherwise, we are restarting
  // the calculation from somewhere in the middle of the dispersity mesh,
  // and we update the value rather than reset it. Similarly for the
  // normalization factor, which is stored as the final value in the
  // results vector (one past the number of q values).
  //
  // The code differs slightly between opencl and dll since opencl is only
  // seeing one q value (stored in the variable "this_F2") while the dll
  // version must loop over all q.
  #if defined(CALL_FQ)
    double weight_norm = (pd_start == 0 ? 0.0 : result[2*nq]);
    double weighted_form = (pd_start == 0 ? 0.0 : result[2*nq+1]);
    double weighted_shell = (pd_start == 0 ? 0.0 : result[2*nq+2]);
    double weighted_radius = (pd_start == 0 ? 0.0 : result[2*nq+3]);
  #else
    double weight_norm = (pd_start == 0 ? 0.0 : result[nq]);
    double weighted_form = (pd_start == 0 ? 0.0 : result[nq+1]);
    double weighted_shell = (pd_start == 0 ? 0.0 : result[nq+2]);
    double weighted_radius = (pd_start == 0 ? 0.0 : result[nq+3]);
  #endif
  #if defined(USE_GPU)
    #if defined(CALL_FQ)
      double this_F2 = (pd_start == 0 ? 0.0 : result[2*q_index+0]);
      double this_F1 = (pd_start == 0 ? 0.0 : result[2*q_index+1]);
    #else
      double this_F2 = (pd_start == 0 ? 0.0 : result[q_index]);
    #endif
  #else // !USE_GPU
    if (pd_start == 0) {
      #ifdef USE_OPENMP
      #pragma omp parallel for
      #endif
      #if defined(CALL_FQ)
          // 2*nq for F^2,F pairs
          for (int q_index=0; q_index < 2*nq; q_index++) result[q_index] = 0.0;
      #else
          for (int q_index=0; q_index < nq; q_index++) result[q_index] = 0.0;
      #endif
    }
    //if (q_index==0) printf("start %d %g %g\n", pd_start, pd_norm, result[0]);
#endif // !USE_GPU


// ====== macros to set up the parts of the loop =======
/*
Based on the level of the loop, uses C preprocessor magic to construct
level-specific looping variables, including these from loop level 3:
  int n3 : length of loop for mesh level 3
  int i3 : current position in the loop for level 3, which is calculated
       from a combination of pd_start, pd_stride[3] and pd_length[3].
  int p3 : is the index into the parameter table for mesh level 3
  double v3[] : pointer into dispersity array to values for loop 3
  double w3[] : pointer into dispersity array to weights for loop 3
  double weight3 : the product of weights from levels 3 and up, computed
       as weight5*weight4*w3[i3].  Note that we need an outermost
       value weight5 set to 1.0 for this to work properly.
After expansion, the loop struction will look like the following:
  // --- PD_INIT(4) ---
  const int n4 = pd_length[4];
  const int p4 = pd_par[4];
  pglobal const double *v4 = pd_value + pd_offset[4];
  pglobal const double *w4 = pd_weight + pd_offset[4];
  int i4 = (pd_start/pd_stride[4])%n4;  // position in level 4 at pd_start
  // --- PD_INIT(3) ---
  const int n3 = pd_length[3];
  ...
  int i3 = (pd_start/pd_stride[3])%n3;  // position in level 3 at pd_start
  PD_INIT(2)
  PD_INIT(1)
  PD_INIT(0)
  // --- PD_OUTERMOST_WEIGHT(5) ---
  const double weight5 = 1.0;
  // --- PD_OPEN(4,5) ---
  while (i4 < n4) {
    parameter[p4] = v4[i4];  // set the value for pd parameter 4 at this mesh point
    const double weight4 = w4[i4] * weight5;
    // from PD_OPEN(3,4)
    while (i3 < n3) {
      parameter[p3] = v3[i3];  // set the value for pd parameter 3 at this mesh point
      const double weight3 = w3[i3] * weight4;
      PD_OPEN(3,2)
      PD_OPEN(2,1)
      PD_OPEN(0,1)
      // ... main loop body ...
      APPLY_PROJECTION    // convert jitter values to spherical coords
      BUILD_ROTATION      // construct the rotation matrix qxy => qabc
      for each q
          FETCH_Q         // set qx,qy from the q input vector
          APPLY_ROTATION  // convert qx,qy to qa,qb,qc
          CALL_KERNEL     // F2 = Iqxy(qa, qb, qc, p1, p2, ...)
      ++step;  // increment counter representing position in dispersity mesh
      PD_CLOSE(0)
      PD_CLOSE(1)
      PD_CLOSE(2)
      // --- PD_CLOSE(3) ---
      if (step >= pd_stop) break;
      ++i3;
    }
    i3 = 0; // reset loop counter for next round through the loop
    // --- PD_CLOSE(4) ---
    if (step >= pd_stop) break;
    ++i4;
  }
  i4 = 0; // reset loop counter even though no more rounds through the loop
*/


// ** prepare inner loops **

// Depending on the shape type (radial, axial, triaxial), the variables
// and calling parameters in the loop body will be slightly different.
// Macros capture the differences in one spot so the rest of the code
// is easier to read. The code below both declares variables for the
// inner loop and defines the macros that use them.


#if defined(CALL_FQ) // COMPUTE_F1_F2 is true
  // unoriented 1D returning <F> and <F^2>
  // Note that F1 and F2 are returned from CALL_FQ by reference, and the
  // user of the CALL_KERNEL macro below is assuming that F1 and F2 are defined.
  double qk;
  double F1, F2;
  #define FETCH_Q() do { qk = q[q_index]; } while (0)
  #define BUILD_ROTATION() do {} while(0)
  #define APPLY_ROTATION() do {} while(0)
  #define CALL_KERNEL() CALL_FQ(qk,F1,F2,local_values.table)

#elif defined(CALL_FQ_A)
  // unoriented 2D return <F> and <F^2>
  // Note that the CALL_FQ_A macro is computing _F1_slot and _F2_slot by
  // reference then returning _F2_slot.  We are calling them _F1_slot and
  // _F2_slot here so they don't conflict with _F1 and _F2 in the macro
  // expansion, or with the use of F2 = CALL_KERNEL() when it is used below.
  double qx, qy;
  double _F1_slot, _F2_slot;
  #define FETCH_Q() do { qx = q[2*q_index]; qy = q[2*q_index+1]; } while (0)
  #define BUILD_ROTATION() do {} while(0)
  #define APPLY_ROTATION() do {} while(0)
  #define CALL_KERNEL() CALL_FQ_A(sqrt(qx*qx+qy*qy),_F1_slot,_F2_slot,local_values.table)

#elif defined(CALL_IQ)
  // unoriented 1D return <F^2>
  double qk;
  #define FETCH_Q() do { qk = q[q_index]; } while (0)
  #define BUILD_ROTATION() do {} while(0)
  #define APPLY_ROTATION() do {} while(0)
  #define CALL_KERNEL() CALL_IQ(qk,local_values.table)

#elif defined(CALL_IQ_A)
  // unoriented 2D
  double qx, qy;
  #define FETCH_Q() do { qx = q[2*q_index]; qy = q[2*q_index+1]; } while (0)
  #define BUILD_ROTATION() do {} while(0)
  #define APPLY_ROTATION() do {} while(0)
  #define CALL_KERNEL() CALL_IQ_A(sqrt(qx*qx+qy*qy), local_values.table)

#elif defined(CALL_IQ_AC)
  // oriented symmetric 2D
  double qx, qy;
  #define FETCH_Q() do { qx = q[2*q_index]; qy = q[2*q_index+1]; } while (0)
  double qa, qc;
  QACRotation rotation;
  // theta, phi, dtheta, dphi are defined below in projection to avoid repeated code.
  #define BUILD_ROTATION() qac_rotation(&rotation, theta, phi, dtheta, dphi);
  #define APPLY_ROTATION() qac_apply(&rotation, qx, qy, &qa, &qc)
  #define CALL_KERNEL() CALL_IQ_AC(qa, qc, local_values.table)

#elif defined(CALL_IQ_ABC)
  // oriented asymmetric 2D
  double qx, qy;
  #define FETCH_Q() do { qx = q[2*q_index]; qy = q[2*q_index+1]; } while (0)
  double qa, qb, qc;
  QABCRotation rotation;
  // theta, phi, dtheta, dphi are defined below in projection to avoid repeated code.
  // psi and dpsi are only for IQ_ABC, so they are processed here.
  const double psi = values[details->theta_par+4];
  local_values.table.psi = 0.;
  #define BUILD_ROTATION() qabc_rotation(&rotation, theta, phi, psi, dtheta, dphi, local_values.table.psi)
  #define APPLY_ROTATION() qabc_apply(&rotation, qx, qy, &qa, &qb, &qc)
  #define CALL_KERNEL() CALL_IQ_ABC(qa, qb, qc, local_values.table)

#elif defined(CALL_IQ_XY)
  // direct call to qx,qy calculator
  double qx, qy;
  #define FETCH_Q() do { qx = q[2*q_index]; qy = q[2*q_index+1]; } while (0)
  #define BUILD_ROTATION() do {} while(0)
  #define APPLY_ROTATION() do {} while(0)
  #define CALL_KERNEL() CALL_IQ_XY(qx, qy, local_values.table)
#endif

// Define APPLY_PROJECTION depending on model symmetries. We do this outside
// the previous if block so that we don't need to repeat the identical
// logic in the IQ_AC and IQ_ABC branches.  This will become more important
// if we implement more projections, or more complicated projections.
#if defined(CALL_IQ) || defined(CALL_IQ_A) || defined(CALL_FQ) || defined(CALL_FQ_A)
  // no orientation
  #define APPLY_PROJECTION() const double weight=weight0
#elif defined(CALL_IQ_XY) // pass orientation to the model
  // CRUFT: support oriented model which define Iqxy rather than Iqac or Iqabc
  // Need to plug the values for the orientation angles back into parameter
  // table in case they were overridden by the orientation offset.  This
  // means that orientation dispersity will not work for these models, but
  // it was broken anyway, so no matter.  Still want to provide Iqxy in case
  // the user model wants full control of orientation/magnetism.
  #if defined(HAVE_PSI)
    const double theta = values[details->theta_par+2];
    const double phi = values[details->theta_par+3];
    const double psi = values[details->theta_par+4];
    double weight;
    #define APPLY_PROJECTION() do { \
      local_values.table.theta = theta; \
      local_values.table.phi = phi; \
      local_values.table.psi = psi; \
      weight=weight0; \
    } while (0)
  #elif defined(HAVE_THETA)
    const double theta = values[details->theta_par+2];
    const double phi = values[details->theta_par+3];
    double weight;
    #define APPLY_PROJECTION() do { \
      local_values.table.theta = theta; \
      local_values.table.phi = phi; \
      weight=weight0; \
    } while (0)
  #else
    #define APPLY_PROJECTION() const double weight=weight0
  #endif
#else // apply jitter and view before calling the model
  // Grab the "view" angles (theta, phi, psi) from the initial parameter table.
  const double theta = values[details->theta_par+2];
  const double phi = values[details->theta_par+3];
  // Make sure jitter angle defaults to zero if there is no jitter distribution
  local_values.table.theta = 0.;
  local_values.table.phi = 0.;
  // The "jitter" angles (dtheta, dphi, dpsi) are stored with the
  // dispersity values and copied to the local parameter table as
  // we go through the mesh.
  double dtheta, dphi, weight;
  #if PROJECTION == 1 // equirectangular
    #define APPLY_PROJECTION() do { \
      dtheta = local_values.table.theta; \
      dphi = local_values.table.phi; \
      weight = fabs(cos(dtheta*M_PI_180)) * weight0; \
    } while (0)
  #elif PROJECTION == 2 // sinusoidal
    #define APPLY_PROJECTION() do { \
      dtheta = local_values.table.theta; \
      dphi = local_values.table.phi; \
      weight = weight0; \
      if (dtheta != 90.0) dphi /= cos(dtheta*M_PI_180); \
      else if (dphi != 0.0) weight = 0.; \
      if (fabs(dphi) >= 180.) weight = 0.; \
    } while (0)
  #endif
#endif // done defining APPLY_PROJECTION

// ** define looping macros **

// Define looping variables
#define PD_INIT(_LOOP) \
  const int n##_LOOP = details->pd_length[_LOOP]; \
  const int p##_LOOP = details->pd_par[_LOOP]; \
  pglobal const double *v##_LOOP = pd_value + details->pd_offset[_LOOP]; \
  pglobal const double *w##_LOOP = pd_weight + details->pd_offset[_LOOP]; \
  int i##_LOOP = (pd_start/details->pd_stride[_LOOP])%n##_LOOP;

// Jump into the middle of the dispersity loop
#define PD_OPEN(_LOOP,_OUTER) \
  while (i##_LOOP < n##_LOOP) { \
    local_values.vector[p##_LOOP] = v##_LOOP[i##_LOOP]; \
    const double weight##_LOOP = w##_LOOP[i##_LOOP] * weight##_OUTER;

// create the variable "weight#=1.0" where # is the outermost level+1 (=MAX_PD).
#define _PD_OUTERMOST_WEIGHT(_n) const double weight##_n = 1.0;
#define PD_OUTERMOST_WEIGHT(_n) _PD_OUTERMOST_WEIGHT(_n)

// Close out the loop
#define PD_CLOSE(_LOOP) \
    if (step >= pd_stop) break; \
    ++i##_LOOP; \
  } \
  i##_LOOP = 0;

// ====== construct the loops =======

// Pointers to the start of the dispersity and weight vectors, if needed.
#if MAX_PD>0
  pglobal const double *pd_value = values + NUM_VALUES;
  pglobal const double *pd_weight = pd_value + details->num_weights;
#endif

// The variable "step" is the current position in the dispersity loop.
// It will be incremented each time a new point in the mesh is accumulated,
// and used to test whether we have reached pd_stop.
int step = pd_start;

// *** define loops for each of 0, 1, 2, ..., modelinfo.MAX_PD-1 ***

// define looping variables
#if MAX_PD>4
  PD_INIT(4)
#endif
#if MAX_PD>3
  PD_INIT(3)
#endif
#if MAX_PD>2
  PD_INIT(2)
#endif
#if MAX_PD>1
  PD_INIT(1)
#endif
#if MAX_PD>0
  PD_INIT(0)
#endif

// open nested loops
PD_OUTERMOST_WEIGHT(MAX_PD)
#if MAX_PD>4
  PD_OPEN(4,5)
#endif
#if MAX_PD>3
  PD_OPEN(3,4)
#endif
#if MAX_PD>2
  PD_OPEN(2,3)
#endif
#if MAX_PD>1
  PD_OPEN(1,2)
#endif
#if MAX_PD>0
  PD_OPEN(0,1)
#endif

//if (q_index==0) {printf("step:%d of %d, pars:",step,pd_stop); for (int i=0; i < NUM_PARS; i++) printf("p%d=%g ",i, local_values.vector[i]); printf("\n");}

  // ====== loop body =======
  // Note: recalc all intermediates for each loop, because we don't know order
  TRANSLATION_VARS(local_values.table);
  if (VALID(local_values.table)) {
     APPLY_PROJECTION();

    // Accumulate I(q)
    // Note: weight==0 must always be excluded
    if (weight > cutoff) {
      double form, shell;
      CALL_VOLUME(form, shell, local_values.table);
      weight_norm += weight;
      weighted_form += weight * form;
      weighted_shell += weight * shell;
      if (radius_effective_mode != 0) {
        weighted_radius += weight * CALL_RADIUS_EFFECTIVE(radius_effective_mode, local_values.table);
      }
      BUILD_ROTATION();

#if !defined(USE_GPU)
      // DLL needs to explicitly loop over the q values.
      #ifdef USE_OPENMP
      #pragma omp parallel for
      #endif
      for (q_index=0; q_index<nq; q_index++)
#endif // !USE_GPU
      {

        FETCH_Q();
        APPLY_ROTATION();

        // ======= COMPUTE SCATTERING ==========
        #if defined(MAGNETIC) && NUM_MAGNETIC > 0
          // Compute the scattering from the magnetic cross sections.
          double F2 = 0.0;
          const double qsq = qx * qx + qy * qy;
          if (qsq > 1.e-16) {
            // TODO: what is the magnetic scattering at q = 0
            //const double px = sin_mspin * cos_mphi;
            //const double py = sin_mspin * sin_mphi;
            //const double pz = cos_mspin;

            // loop over uu, ud real, du real, dd, ud imag, du imag
            for (unsigned int xs = 0; xs < 6; xs++) {
              const double xs_weight = xs_weights[xs];
              if (xs_weight > 1.e-8) {
                // Since the cross section weight is significant, set the slds
                // to the effective slds for this cross section, call the
                // kernel, and add according to weight.
                for (int sk = 0; sk<NUM_MAGNETIC; sk++) {
                  const int32_t mag_index = NUM_PARS + 6 + 3 * sk;
                  const int32_t sld_index = slds[sk];
                  const double mx = values[mag_index];
                  const double my = values[mag_index + 1];
                  const double mz = values[mag_index + 2];
                  local_values.vector[sld_index] =
                    mag_sld(xs, qx, qy, cos_mspin, sin_mspin, cos_mphi, sin_mphi, values[sld_index + 2], mx, my, mz);
//if (q_index==0) printf("%d: (qx,qy)=(%g,%g) xs=%d sld%d=%g p=(%g,%g) m=(%g,%g,%g)\n",
//  q_index, qx, qy, xs, sk, local_values.vector[sld_index], px, py, mx, my, mz);
                }
                F2 += xs_weight * CALL_KERNEL();
              }
            }
          }
        #else  // !MAGNETIC
          #if defined(CALL_FQ)
            CALL_KERNEL(); // sets F1 and F2 by reference
          #else
            const double F2 = CALL_KERNEL();
          #endif
        #endif // !MAGNETIC
//printf("q_index:%d %g %g %g %g\n", q_index, F2, weight0);

        #if defined(USE_GPU)
          #if defined(CALL_FQ)
            this_F2 += weight * F2;
            this_F1 += weight * F1;
          #else
            this_F2 += weight * F2;
          #endif
        #else // !USE_OPENCL
          #if defined(CALL_FQ)
            result[2*q_index+0] += weight * F2;
            result[2*q_index+1] += weight * F1;
          #else
            result[q_index] += weight * F2;
          #endif
        #endif // !USE_OPENCL
      }
    }
  }
// close nested loops
++step;
#if MAX_PD>0
  PD_CLOSE(0)
#endif
#if MAX_PD>1
  PD_CLOSE(1)
#endif
#if MAX_PD>2
  PD_CLOSE(2)
#endif
#if MAX_PD>3
  PD_CLOSE(3)
#endif
#if MAX_PD>4
  PD_CLOSE(4)
#endif

// Remember the results and the updated norm.
#if defined(USE_GPU)
  #if defined(CALL_FQ)
  result[2*q_index+0] = this_F2;
  result[2*q_index+1] = this_F1;
  #else
  result[q_index] = this_F2;
  #endif
  if (q_index == 0)
#endif
  {
#if defined(CALL_FQ)
    result[2*nq] = weight_norm;
    result[2*nq+1] = weighted_form;
    result[2*nq+2] = weighted_shell;
    result[2*nq+3] = weighted_radius;
#else
    result[nq] = weight_norm;
    result[nq+1] = weighted_form;
    result[nq+2] = weighted_shell;
    result[nq+3] = weighted_radius;
#endif
  }

// ** clear the macros in preparation for the next kernel **
#undef PD_INIT
#undef PD_OPEN
#undef PD_CLOSE
#undef FETCH_Q
#undef APPLY_PROJECTION
#undef BUILD_ROTATION
#undef APPLY_ROTATION
#undef CALL_KERNEL
}
