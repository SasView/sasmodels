double clipp(double value, double low, double high) //from kernel_iq.c
{
  return (value < low ? low : (value > high ? high : value));
} 

//cylinder cross section for nuclear scattering

static double
form_volume(double radius, double length)
{
    return M_PI*radius*radius*length;
}


static double
fq(double x, double y, double z,  double radius, double length, double sld, double solvent_sld)
{


  const double qsq=y*y+z*z; 
  const double q=sqrt(qsq); 
    return (sld-solvent_sld) *sas_2J1x_x(q*radius) * sas_sinx_x(x*0.5*length);
}




//Definition of special integral function on Bessel functions see Supplemental Material of Metlov and Michels Sci. Rep. 6, 25055 (2016)

static double total_F2(double k)
{
	return  0.5 * k * M_PI * (sas_J1(k)* struve(0,k)- sas_J0(k)*struve(1,k));
}

static double total_F6(double k)
{
	return   0.5 *k* (2 *k* sas_J1(k)-M_PI *sas_J1(k)* struve(0,k)+M_PI* sas_J0(k)* struve(1,k));
}

static double total_F7(double k)
{
	return  0.5*k* (M_PI *sas_J1(k)* struve(0,k)+sas_J0(k)*(2-M_PI* struve(1,k)));
}

static double total_F8(double k)
{
	return  0.5 *(sas_J1(k)*(-2+k *M_PI* struve(0,k))+k *sas_J0(k)*(2-M_PI* struve(1,k)));
}

static double total_F9(double k)
{
	return  -0.5*k* (2* k *sas_J1(k)-3*M_PI* sas_J1(k)* struve(0,k)+3*M_PI* sas_J0(k)* struve(1,k));
}

static double total_F10(double k)
{
	return  0.5*(sas_J1(k) *(-4+k*M_PI*struve(0,k))+k *sas_J0(k)* (2-M_PI* struve(1,k)));
}


static double muy1(double cos_alpha,double sin_alpha, double Rm, double k)
{
	return 2*cos_alpha*sin_alpha* (-Rm/2/k*total_F10(k)+1/(2*Rm*k*k*k)*total_F9(k));
}

static double muz1(double cos_alpha,double sin_alpha, double Rm, double k)
{
	return 1/(Rm*k*k*k)*(cos_alpha*cos_alpha*(-total_F6(k)+k*k*Rm*Rm*total_F7(k))+(cos_alpha*cos_alpha-sin_alpha*sin_alpha)*(total_F2(k)-total_F8(k)));
}

static double muperp0( double k)
{
	return total_F2(k)/k/k;
}

static double muz0(double cos_alpha,double sin_alpha, double k)
{
	return muperp0(k)*sin_alpha; //complex quantity does not mix with first order terms
}

static double muy0(double cos_alpha,double sin_alpha, double k)
{
	return -muperp0(k)*cos_alpha; //complex quantity does not mix with first order terms
}


//Mz is defined as the longitudinal magnetisation component along the magnetic field. Mx is zero neglecting contributions for the vortex core



static double fqMyreal( double x, double y, double z, double L, double R, double Ms, double b, double Rm)
{

  const double qsq=y*y+z*z; 
  const double q=sqrt(qsq); 
  const double cos_alpha=z/q;
  const double sin_alpha=y/q;
  const double k=R*q;	

  const double fx=Ms*L*R*R/sqrt(2*M_PI)*sas_J0(L*x/2); 
  const double f = fx*b*muy1(cos_alpha,sin_alpha, Rm, k);
  return f;
}

static double fqMyimag( double x, double y, double z, double L, double R, double Ms)
{
  const double qsq=y*y+z*z; 
  const double q=sqrt(qsq); 
  const double cos_alpha=z/q;
  const double sin_alpha=y/q;	
    const double k=R*q;

  const double fx=Ms*L*R*R/sqrt(2*M_PI)*sas_J0(L*x/2); 
  const double f = fx*muy0(cos_alpha, sin_alpha, k);
  return f;
}

static double fqMzreal( double x, double y, double z, double L,  double R, double Ms, double b, double Rm)
{
  const double qsq=y*y+z*z; 
  const double q=sqrt(qsq); 
  const double cos_alpha=z/q;
  const double sin_alpha=y/q;	
  const double k=R*q;

  const double fx=Ms*L*R*R/sqrt(2*M_PI)*sas_J0(L*x/2); 
  const double f = fx*b*muz1(cos_alpha,sin_alpha, Rm, k);
  return f;
}

static double fqMzimag( double x, double y, double z, double L, double R, double Ms)
{
  const double qsq=y*y+z*z; 
  const double q=sqrt(qsq); 
  const double cos_alpha=z/q;
  const double sin_alpha=y/q;	
  const double k=R*q;

  const double fx=Ms*L*R*R/sqrt(2*M_PI)*sas_J0(L*x/2); 
  const double f = fx*muz0(cos_alpha, sin_alpha, k);
  return f;
}

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

double SCALAR_VEC( double *vec1, double *vec2)
{
 return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

double MAG_VEC( double *vec)
{
 return sqrt(SCALAR_VEC(vec,vec));
}

void ORTH_VEC(double *result_vec, double *vec1, double *vec2)
{
 result_vec[0] = vec1[0] - SCALAR_VEC(vec1,vec2) / SCALAR_VEC(vec2,vec2) * vec2[0];
 result_vec[1] = vec1[1] - SCALAR_VEC(vec1,vec2) / SCALAR_VEC(vec2,vec2) * vec2[1];
 result_vec[2] = vec1[2] - SCALAR_VEC(vec1,vec2) / SCALAR_VEC(vec2,vec2) * vec2[2];
}


//Weighting of spin resolved cross sections to reconstruct partially polarised beam with imperfect optics using up_i/up_f.
void set_weights(double in_spin, double out_spin, double weight[8]) //from kernel_iq.c
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


   weight[0] = (1.0-in_spin) * (1.0-out_spin) / norm; // dd.real
   weight[1] = weight[0]; // dd.imag
   weight[2] = in_spin * out_spin / norm;             // uu.real
   weight[3] = weight[2];             // uu.imag
   weight[4] = (1.0-in_spin) * out_spin / norm;       // du.real
   weight[5] = weight[4];       // du.imag 
   weight[6] = in_spin * (1.0-out_spin) / norm;       // ud.real
   weight[7] = weight[6];       // ud.imag
 }

//transforms scattering vector q in polarisation/magnetisation coordinate system
void set_scatvec(double *qmag,double q,
  double cos_theta, double sin_theta,
  double alpha, double beta) {
  double cos_alpha, sin_alpha;
  double cos_beta, sin_beta;

  SINCOS(alpha*M_PI/180, sin_alpha, cos_alpha);
  SINCOS(beta*M_PI/180, sin_beta, cos_beta);
    //field is defined along (0,0,1), orientation of detector
    //is precessing in a cone around B with an inclination of theta

  qmag[0] = q*(cos_alpha * cos_theta);
  qmag[1] = q*(cos_theta * sin_alpha*sin_beta + 
  cos_beta * sin_theta);
  qmag[2] = q*(-cos_beta * cos_theta* sin_alpha + 
  sin_beta * sin_theta);

}


//Evaluating the magnetic scattering vector (Halpern Johnson vector) for general orientation of q and collecting terms for the spin-resolved (POLARIS) cross sections. Mz is along the applied magnetic field direction, which is also the polarisation direction.
void mag_sld(
   // 0=dd.real, 1=dd.imag, 2=uu.real, 3=uu.imag,  4=du.real, 5=du.imag,  6=ud.real, 7=ud.imag
 double x, double y, double z,
  double mxreal, double mximag, double myreal,  double myimag, double mzreal,double mzimag, double nuc, double sld[8])
{
  double vector[3];
  //The (transversal) magnetisation and hence the magnetic scattering sector is here a complex quantity. The spin-flip (magnetic) scattering amplitude is given with
  //  MperpPperpQ \pm i MperpP  (Moon-Riste-Koehler Phys Rev 181, 920, 1969) with Mperp and MperpPperpQ the
  //magnetisation scattering vector components perpendicular to the polarisation/field direction. Collecting terms in SF that are real (MperpPperpQreal + SCALAR_VEC(MperpPimag,qvector) )  and imaginary (MperpPperpQimag \pm SCALAR_VEC(MperpPreal,qvector) )
  double Mvectorreal[3];
  double Mvectorimag[3];
  double Pvector[3];
  double Mperpreal[3];
  double MperpPreal[3];
  double MperpPperpQreal[3];
  double Mperpimag[3];
  double MperpPimag[3];
  double MperpPperpQimag[3];
  
  const double q = sqrt(x*x + y*y + z*z);
  SET_VEC(vector, x/q, y/q, z/q); 

   //Moon-Riste-Koehler choose z as pointing along field/polarisation axis  
  SET_VEC(Mvectorreal, 0, myreal, mzreal);
  SET_VEC(Mvectorimag, 0, myimag, mzimag);
  SET_VEC(Pvector, 0, 0, 1);
   //Magnetic scattering vector Mperp could be simplified like in Moon-Riste-Koehler
   //leave the generic computation just to check
  ORTH_VEC(Mperpreal, Mvectorreal, vector);
  ORTH_VEC(MperpPreal, Mperpreal, Pvector);
  ORTH_VEC(MperpPperpQreal, MperpPreal, vector);
  ORTH_VEC(Mperpimag, Mvectorimag, vector);
  ORTH_VEC(MperpPimag, Mperpimag, Pvector);
  ORTH_VEC(MperpPperpQimag, MperpPimag, vector);
  
  
   sld[0] = nuc - SCALAR_VEC(Pvector,Mperpreal); // dd.real => sld - D Pvector \cdot Mperp 
   sld[1] = +SCALAR_VEC(Pvector,Mperpimag); //dd.imag  = nuc_img - SCALAR_VEC(Pvector,Mperpimg); nuc_img only exist for noncentrosymmetric nuclear structures; 
   sld[2] = nuc + SCALAR_VEC(Pvector,Mperpreal);              // uu => sld + D Pvector \cdot Mperp
   sld[3] = -SCALAR_VEC(Pvector,Mperpimag); //uu.imag

   sld[4] = MAG_VEC(MperpPperpQreal)+SCALAR_VEC(MperpPimag,vector);       // du.real => length of vector MperpPperpQ:  | MperpP- (MperpP\cdot qvector)  qvector | and Mperp perpendicular to P and along q
   sld[5] = MAG_VEC(MperpPperpQimag)-SCALAR_VEC(MperpPreal,vector);       // du.imag => | MperpPimag- (MperpPimag\cdot qvector)  qvector | - i MperpPreal \cdot qvector
   sld[6] = MAG_VEC(MperpPperpQreal)-SCALAR_VEC(MperpPimag,vector);      // ud.real =>  length of vector MperpPperpQ
   sld[7] = MAG_VEC(MperpPperpQimag)+SCALAR_VEC(MperpPreal,vector);         // du.imag => | MperpPimag- (MperpPimag\cdot qvector)  qvector | + i MperpPreal \cdot qvector

 }


//calculate 2D from _fq
 static double
 Iqxy(double qx, double qy, double radius, double length, double core_nuc, double solvent_nuc, double magnetic_sld_disc,double b, double R_m,  double up_i, double up_f, double alpha, double beta)
 {
   const double q = sqrt(qx*qx + qy*qy);
   if (q > 1.0e-16 ) {
     const double cos_theta=qx/q;
     const double sin_theta=qy/q;

     double qmag[3];
     set_scatvec(qmag,q,cos_theta, sin_theta, alpha, beta);

   double weights[8];  // 0=dd.real, 1=dd.imag, 2=uu.real, 3=uu.imag,  4=du.real, 6=du.imag,  7=ud.real, 5=ud.imag
   set_weights(up_i, up_f, weights);
   
   
   double nuc=fq(qmag[0], qmag[1], qmag[2], radius, length, core_nuc, solvent_nuc);


   double sld[8];




      double myreal=fqMyreal(qmag[0], qmag[1], qmag[2], length, radius, magnetic_sld_disc, b, R_m);
      double myimag=fqMyimag(qmag[0], qmag[1], qmag[2], length, radius, magnetic_sld_disc);
      double mzreal=fqMzreal(qmag[0], qmag[1], qmag[2], length, radius, magnetic_sld_disc, b, R_m);
      double mzimag=fqMzimag(qmag[0], qmag[1], qmag[2], length, radius, magnetic_sld_disc);


 mag_sld(qmag[0], qmag[1], qmag[2], 0, 0, myreal, myimag, mzreal, mzimag, nuc, sld);

      double form = 0.0;
      for (unsigned int xs=0; xs<8; xs++) {
       if (weights[xs] > 1.0e-8) {
                // Since the cross section weight is significant, set the slds
                // to the effective slds for this cross section, call the
                // kernel, and add according to weight.
      // loop over uu, ud real, du real, dd, ud imag, du imag 
        form += weights[xs]*sld[xs]*sld[xs];
      }
    }


  return 0.5*1.0e-4*form;}
}



//calculate 1D by averaging over theta
//TODO: choose orientation and sector width for averaging
// 2D to 1D
static double
Iq(double q, double radius, double length, double core_nuc, double solvent_nuc, double magnetic_sld_disc,double b, double R_m,  double up_i, double up_f, double alpha, double beta)

{
   double sin_theta, cos_theta; // slots to hold sincos function output of the orientation on the detector plane
   double total_F1D = 0.0;
   for (int j=0; j<GAUSS_N ;j++) {

      const double theta = M_PI * (GAUSS_Z[j] + 1.0); // 0 .. 2 pi
      SINCOS(theta, sin_theta, cos_theta);

      double qmag[3];
      set_scatvec(qmag,q,cos_theta, sin_theta, alpha, beta);

   double weights[8];  // 0=dd.real, 1=dd.imag, 2=uu.real, 3=uu.imag,  4=du.real, 5=du.imag,  6=ud.real, 7=ud.imag
   set_weights(up_i, up_f, weights);
   
   double nuc=fq(qmag[0], qmag[1], qmag[2], radius, length, core_nuc, solvent_nuc);

   double sld[8];

        //loop over random anisotropy axis with isotropic orientation gamma for Hkx and Hky
//To be modified for textured material see also Weissmueller et al. PRB 63, 214414 (2001)


      double myreal=fqMyreal(qmag[0], qmag[1], qmag[2], length, radius, magnetic_sld_disc, b, R_m);
      double myimag=fqMyimag(qmag[0], qmag[1], qmag[2], length, radius, magnetic_sld_disc);
      double mzreal=fqMzreal(qmag[0], qmag[1], qmag[2], length, radius, magnetic_sld_disc, b, R_m);
      double mzimag=fqMzimag(qmag[0], qmag[1], qmag[2], length, radius, magnetic_sld_disc);


      
      mag_sld(qmag[0], qmag[1], qmag[2], 0, 0, myreal, myimag, mzreal, mzimag, nuc, sld);
  
      double form = 0.0;
      for (unsigned int xs=0; xs<8; xs++) {
       if (weights[xs] > 1.0e-8 ) {
                // Since the cross section weight is significant, set the slds
                // to the effective slds for this cross section, call the
                // kernel, and add according to weight.
      // loop over uu, ud real, du real, dd, ud imag, du imag 
        form += weights[xs]*sld[xs]*sld[xs];
      }
        total_F1D += GAUSS_W[j] * form ;
    }

   //convert from [1e-12 A-1] to [cm-1]
return 0.25*1.0e-4*total_F1D;
}}






