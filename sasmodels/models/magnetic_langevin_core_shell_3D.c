  static double clipp(double value, double low, double high) //from kernel_iq.c
{
  return (value < low ? low : (value > high ? high : value));
} 

static double langevin(
    double x) {
    // cotanh(x)-1\x

    if (x < 0.00001) {
        // avoid dividing by zero
        return 1.0/3.0*x;
    } else {
        return 1.0/tanh(x)-1/x;
    }
}

static double langevinoverx(
    double x) {
    // cotanh(x)-1\x

    if (x < 0.00001) {
        // avoid dividing by zero
        return 1.0/3.0;
    } else {
        return langevin(x)/x;
    }
}

//from coremultishell.c



static double
outer_radius(double radius, double fp_n, double thickness[])
{
  double r = radius;
  int n = (int)(fp_n+0.5);
  for (int i=0; i < n; i++) {
    r += thickness[i];
  }
  return r;
}

static double
form_volume(double radius, double fp_n, double thickness[])
{
  return M_4PI_3 * cube(outer_radius(radius, fp_n, thickness));
}


effective_radius(int mode, double core_radius, double fp_n, double thickness[])
{
  switch (mode) {
  default:
  case 1: // outer radius
    return outer_radius(core_radius, fp_n, thickness);
  case 2: // core radius
    return core_radius;
  }
}


static double fqnucsq(double q, double sld_core, double radius,
   double sld_solvent, double fp_n, double sld[], double thickness[])
{
      const int n = (int)(fp_n+0.5);
  double f, r, last_sld;
  r = radius;
  last_sld = sld_core;
  f = 0.;
  for (int i=0; i<n; i++) {
    f += M_4PI_3 * cube(r) * (sld[i] - last_sld) * sas_3j1x_x(q*r);
    last_sld = sld[i];
    r += thickness[i];
  }
  f += M_4PI_3 * cube(r) * (sld_solvent - last_sld) * sas_3j1x_x(q*r);
    return 1e-4 *f*f;
}

static double fqMtranssq(double q, double sld_core, double radius,
   double sld_solvent, double fp_n, double sld[], double thickness[],double eta_core,double eta_solvent,double eta[], double delta_solvent,double delta[])
{
      const int n = (int)(fp_n+0.5);
  double f,fsq_perp, r, last_sld;
  r = radius;
  
  last_sld = sld_core*sqrt(langevinoverx(eta_core));
  f = 0.;
  fsq_perp=0;
  for (int i=0; i<n; i++) {
    f += M_4PI_3 * cube(r) *(sld[i]*sqrt(langevinoverx(eta[i])*delta[i]) - last_sld) * sas_3j1x_x(q*r);//transversal magnetisation component coaligned with core transversal magnetisation
    fsq_perp += square(M_4PI_3 * cube(r) *sld[i]*sqrt(langevinoverx(eta[i])*(1-delta[i])) * sas_3j1x_x(q*r));//transversal magnetisation component perpendicular with core transversal magnetisation, adding incoherently (only of squared amplitude)
    last_sld = sld[i]*sqrt(langevinoverx(eta[i])*delta[i]);   
    r += thickness[i];
  }
  f += M_4PI_3 * cube(r) * (sld_solvent*sqrt(langevinoverx(eta_solvent)*delta_solvent) - last_sld) * sas_3j1x_x(q*r);//transversal magnetisation component coaligned with core magnetisation
  fsq_perp += square(M_4PI_3 * cube(r) * sld_solvent*sqrt(langevinoverx(eta_solvent)*(1-delta_solvent))  * sas_3j1x_x(q*r));
    return 1e-4 *(f*f +fsq_perp);
}

static double fqMzsq(double q, double sld_core, double radius,
   double sld_solvent, double fp_n, double sld[], double thickness[],double eta_core,double eta_solvent,double eta[])
{
      const int n = (int)(fp_n+0.5);
  double f, r, last_sld;
  r = radius;
  last_sld = sld_core*sqrt(1-2*langevinoverx(eta_solvent));// sqrt() needed for correct scale to intensity later
  f = 0.;
  for (int i=0; i<n; i++) {
    f += M_4PI_3 * cube(r) * (sld[i]*sqrt(1-2*langevinoverx(eta_solvent)) - last_sld) * sas_3j1x_x(q*r);
    last_sld = sld[i]*sqrt(1-2*langevinoverx(eta[i]));
    r += thickness[i];
  }
  f += M_4PI_3 * cube(r) * (sld_solvent*sqrt(1-2*langevinoverx(eta_solvent)) - last_sld) * sas_3j1x_x(q*r);
    return 1e-4 *f*f;
}

static double fqnuc(double q, double sld_core, double radius,
   double sld_solvent, double fp_n, double sld[], double thickness[])
{
      const int n = (int)(fp_n+0.5);
  double f, r, last_sld;
  r = radius;
  last_sld = sld_core;
  f = 0.;
  for (int i=0; i<n; i++) {
    f += M_4PI_3 * cube(r) * (sld[i] - last_sld) * sas_3j1x_x(q*r);
    last_sld = sld[i];
    r += thickness[i];
  }
  f += M_4PI_3 * cube(r) * (sld_solvent - last_sld) * sas_3j1x_x(q*r);
    return 1e-2 *f;
}

static double fqMz(double q, double sld_core, double radius,
   double sld_solvent, double fp_n, double sld[], double thickness[],double eta_core,double eta_solvent,double eta[])
{
      const int n = (int)(fp_n+0.5);
  double f, r, last_sld;
  r = radius;
  last_sld = sld_core*langevin(eta_core);
  f = 0.;
  for (int i=0; i<n; i++) {
    f += M_4PI_3 * cube(r) * (sld[i]*langevin(eta[i]) - last_sld) * sas_3j1x_x(q*r);
    last_sld = sld[i]*langevin(eta[i]);
    r += thickness[i];
  }
  f += M_4PI_3 * cube(r) * (sld_solvent*langevin(eta_solvent) - last_sld) * sas_3j1x_x(q*r);
    return 1e-2 *f;
}




//Magnetization component wrt to field are fixed and defined by Langevin function.
//M_trans=M_x=M_y=0
//M_trans_quad=langevinoverx(xi)
//M_z= langevin(xi)
//M_zquad=(1-2*langevinoverx(xi))
//z axis oriented always along field/polarisation axis.



  
//some basic vector algebra

void SET_VEC(double *vector, double v0, double v1, double v2) {
    vector[0] = v0;
    vector[1] = v1;
    vector[2] = v2;
}

void SCALE_VEC(double *vector, double a) {
    vector[0] = a*vector[0];
    vector[1] = a*vector[1];
    vector[2] = a*vector[2];
}

void ADD_VEC(double *result_vec, double *vec1, double *vec2) {
    result_vec[0] = vec1[0] + vec2[0];
    result_vec[1] = vec1[1] + vec2[1];
    result_vec[2] = vec1[2] + vec2[2];
}



//Rotation of Detector around the magnetic reference frame
//Magnetic coordinates and scattering coordinates coincides when z/field axis along the neutron beam.
//then azimuthal angle theta=0Deg, rotation beta=arbitrary


static double square_mag_scat(
double cos_theta, double sin_theta,
    double alpha, double beta) {
    double cos_alpha, sin_alpha;
    double cos_beta, sin_beta;

    SINCOS(alpha*M_PI/180, sin_alpha, cos_alpha);
    SINCOS(beta*M_PI/180, sin_beta, cos_beta);
    //field is defined along (0,0,1), orientation of detector
    //is precessing in a cone around B with an inclination of theta
    double scat_vector[3];
    double perp_scat_vector[3];

    SET_VEC(scat_vector, cos_alpha* cos_theta, 
cos_theta* sin_alpha*sin_beta + 
 cos_beta* sin_theta, -cos_beta* cos_theta* sin_alpha + 
 sin_beta* sin_theta);



    return square(scat_vector[2]);

}


//weight of the magnetisation vector for the scattering, done by evaluating the magnetic scattering vector (Halpern Johnson vector) for general orientation of q. Inserting result in the spin-resolved scattering cross section
//and collecting terms 

//Calculation of the 4 spin-resolved scattering cross sections for superparamagnetic spheres, field perpendicular to beam 



//spin-resolved (POLARIS) cross sections
//NSF++ (F_N-Mz*(1-rotated_scat_vector[3]^2))^2+ Mx^2*(1-rotated_scat_vector[3]^2)*rotated_scat_vector[3]^2 
static double Idd(double q, double sld_core,double magnetic_sld_core,double eta_core,double radius,
   double sld_solvent,double magnetic_sld_solvent,double eta_solvent,double delta_solvent,
   double fp_n, double sld[],double magnetic_sld[],double eta[], double delta[], double thickness[], double cos_theta, double sin_theta, double alpha, double beta) { 
        
return fqnucsq(q, sld_core, radius, sld_solvent, fp_n, sld, thickness)-fqnuc(q, sld_core, radius, sld_solvent, fp_n, sld, thickness)*fqMz(q, magnetic_sld_core, radius, magnetic_sld_solvent, fp_n, magnetic_sld, thickness,eta_core,eta_solvent,eta)*(1-square_mag_scat(cos_theta, sin_theta,alpha, beta))+fqMzsq(q, magnetic_sld_core, radius, magnetic_sld_solvent, fp_n, magnetic_sld, thickness,eta_core,eta_solvent,eta)*square(1-square_mag_scat(cos_theta, sin_theta,alpha, beta))+fqMtranssq(q, magnetic_sld_core, radius, magnetic_sld_solvent, fp_n, magnetic_sld, thickness,eta_core,eta_solvent,eta, delta_solvent,delta)*(1-square_mag_scat(cos_theta, sin_theta,alpha, beta))*square_mag_scat(cos_theta, sin_theta,alpha, beta);
	}   
	
//NSF-- (F_N+Mz*(1-rotated_scat_vector[3]^2))^2+ Mx^2*(1-rotated_scat_vector[3]^2)*rotated_scat_vector[3]^2 
static double Iuu(double q, double sld_core,double magnetic_sld_core,double eta_core,double radius,
   double sld_solvent,double magnetic_sld_solvent,double eta_solvent,double delta_solvent,
   double fp_n, double sld[],double magnetic_sld[],double eta[], double delta[], double thickness[], double cos_theta, double sin_theta, double alpha, double beta) {
return fqnucsq(q, sld_core, radius, sld_solvent, fp_n, sld, thickness)+fqnuc(q, sld_core, radius, sld_solvent, fp_n, sld, thickness)*fqMz(q, magnetic_sld_core, radius, magnetic_sld_solvent, fp_n, magnetic_sld, thickness,eta_core,eta_solvent,eta)*(1-square_mag_scat(cos_theta, sin_theta,alpha, beta))+fqMzsq(q, magnetic_sld_core, radius, magnetic_sld_solvent, fp_n, magnetic_sld, thickness,eta_core,eta_solvent,eta)*square(1-square_mag_scat(cos_theta, sin_theta,alpha, beta))+fqMtranssq(q, magnetic_sld_core, radius, magnetic_sld_solvent, fp_n, magnetic_sld, thickness,eta_core,eta_solvent,eta, delta_solvent,delta)*(1-square_mag_scat(cos_theta, sin_theta,alpha, beta))*square_mag_scat(cos_theta, sin_theta,alpha, beta);
	}   

//SF Mz^2*(1-rotated_scat_vector[3]^2)*rotated_scat_vector[3]^2 +Mx^2*(1+rotated_scat_vector[3]^4)
static double Idu(double q, double sld_core,double magnetic_sld_core,double eta_core,double radius,
   double sld_solvent,double magnetic_sld_solvent,double eta_solvent,double delta_solvent,
   double fp_n, double sld[],double magnetic_sld[],double eta[], double delta[], double thickness[], double cos_theta, double sin_theta, double alpha, double beta) {
    return fqMtranssq(q, magnetic_sld_core, radius, magnetic_sld_solvent, fp_n, magnetic_sld,thickness,eta_core,eta_solvent,eta, delta_solvent,delta)*(1+square(square_mag_scat(cos_theta, sin_theta,alpha, beta)))+fqMzsq(q, magnetic_sld_core, radius, magnetic_sld_solvent, fp_n,magnetic_sld, thickness,eta_core,eta_solvent,eta)* (1-square_mag_scat(cos_theta, sin_theta,alpha, beta))*square_mag_scat(cos_theta, sin_theta,alpha, beta);
	}   	

static double Iud(double q, double sld_core,double magnetic_sld_core,double eta_core,double radius,
   double sld_solvent,double magnetic_sld_solvent,double eta_solvent,double delta_solvent,
   double fp_n, double sld[],double magnetic_sld[],double eta[], double delta[], double thickness[], double cos_theta, double sin_theta, double alpha, double beta) {
    return Idu(q, sld_core,magnetic_sld_core, eta_core, radius, sld_solvent,magnetic_sld_solvent, eta_solvent,delta_solvent,
   fp_n, sld,magnetic_sld,eta,delta,  thickness, cos_theta, sin_theta, alpha, beta);
	}  



//weighting of spin resolved cross sections to reconstruct partially polarised beam with imperfect optics using up_i/up_f.
static void set_weights(double in_spin, double out_spin, double weight[4]) //from kernel_iq.c
{
  double norm=out_spin;
   
  
  in_spin = clipp(sqrt(square(in_spin)), 0.0, 1.0);//opencl has ambiguities for abs()
  out_spin = clipp(sqrt(square(out_spin)), 0.0, 1.0);

  if (out_spin < 0.5){norm=1-out_spin;}
  


// The norm is needed to make sure that the scattering cross sections are
//correctly weighted, such that the sum of spin-resolved measurements adds up to
// the unpolarised or half-polarised scattering cross section. No intensity weighting
// needed on the incoming polariser side assuming that a user has normalised
// to the incoming flux with polariser in for SANSPOl and unpolarised beam, respectively.


  weight[0] = (1.0-in_spin) * (1.0-out_spin) / norm; // dd
  weight[1] = (1.0-in_spin) * out_spin / norm;       // du
  weight[2] = in_spin * (1.0-out_spin) / norm;       // ud
  weight[3] = in_spin * out_spin / norm;             // uu
}


//calculate 2D from _fq
static double
Iqxy(double qx, double qy, double sld_core,double magnetic_sld_core,double eta_core,double radius,
   double sld_solvent,double magnetic_sld_solvent,double eta_solvent,double delta_solvent,
   double fp_n, double sld[],double magnetic_sld[],double eta[], double delta[], double thickness[], double up_i, double up_f, double alpha, double beta)
{
 const double q = sqrt(qx*qx + qy*qy);
    double cos_theta=qx/q;
    double sin_theta=qy/q;

    double weights[8];  // uu, ud, du, dd, fill,fill, fill, fill (make memory alloc happy)
    set_weights(up_i, up_f, weights);

    
    const double form=weights[0]*Idd(q, sld_core,magnetic_sld_core,eta_core,radius,sld_solvent,magnetic_sld_solvent,eta_solvent,delta_solvent,
   fp_n, sld,magnetic_sld,eta,delta,thickness,cos_theta ,sin_theta, alpha, beta) + weights[1]*Idu(q, sld_core,magnetic_sld_core,eta_core,radius,sld_solvent,magnetic_sld_solvent,eta_solvent,delta_solvent,
   fp_n, sld,magnetic_sld,eta,delta,thickness,cos_theta ,sin_theta, alpha, beta) + weights[2]*Iud(q, sld_core,magnetic_sld_core,eta_core,radius,sld_solvent,magnetic_sld_solvent,eta_solvent,delta_solvent,
   fp_n, sld,magnetic_sld,eta,delta,thickness,cos_theta ,sin_theta, alpha, beta) + weights[3]*Iuu(q, sld_core,magnetic_sld_core,eta_core,radius,sld_solvent,magnetic_sld_solvent,eta_solvent,delta_solvent,
   fp_n, sld,magnetic_sld,eta,delta,thickness,cos_theta ,sin_theta, alpha, beta);
 

   return form;
}



//calculate 1D by averaging over theta
//TODO: choose orientation and sector width for averaging
// 2D to 1D
static double
Iq(double q, double sld_core,double magnetic_sld_core,double eta_core,double radius,
   double sld_solvent,double magnetic_sld_solvent,double eta_solvent,double delta_solvent,
   double fp_n, double sld[],double magnetic_sld[],double eta[], double delta[], double thickness[], double up_i, double up_f, double alpha, double beta)
{
     double weights[8];  // uu, ud, du, dd, fill,fill, fill, fill (make memory alloc happy) 
    set_weights(up_i, up_f, weights);
   double sin_theta, cos_theta; // slots to hold sincos function output of the orientation on the detector plane
  double total_F2 = 0.0;
    for (int i=0; i<GAUSS_N ;i++) {

        const double theta = M_PI * (GAUSS_Z[i] + 1.0); // 0 .. 2 pi
        SINCOS(theta, sin_theta, cos_theta);
        const double form = weights[0]*Idd(q, sld_core,magnetic_sld_core,eta_core,radius,sld_solvent,magnetic_sld_solvent,eta_solvent,delta_solvent,
   fp_n, sld,magnetic_sld,eta,delta,thickness,cos_theta ,sin_theta, alpha, beta) + weights[1]*Idu(q, sld_core,magnetic_sld_core,eta_core,radius,sld_solvent,magnetic_sld_solvent,eta_solvent,delta_solvent,
   fp_n, sld,magnetic_sld,eta,delta,thickness,cos_theta ,sin_theta, alpha, beta) + weights[2]*Iud(q, sld_core,magnetic_sld_core,eta_core,radius,sld_solvent,magnetic_sld_solvent,eta_solvent,delta_solvent,
   fp_n, sld,magnetic_sld,eta,delta,thickness,cos_theta ,sin_theta, alpha, beta) + weights[3]*Iuu(q, sld_core,magnetic_sld_core,eta_core,radius,sld_solvent,magnetic_sld_solvent,eta_solvent,delta_solvent,
   fp_n, sld,magnetic_sld,eta,delta,thickness,cos_theta ,sin_theta, alpha, beta);
       

        total_F2 += GAUSS_W[i] * form ;
    }

    //convert from [1e-12 A-1] to [cm-1]
    return total_F2;
 }






           
