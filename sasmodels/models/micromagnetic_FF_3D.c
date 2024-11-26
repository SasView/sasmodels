//Core-shell form factor for anisotropy field (Hkx, Hky and Hkz), Nuc and lonitudinal magnetization Mz
static double
form_volume(double radius, double thickness)
{
  if( radius + thickness > radius + thickness)
    return M_4PI_3 * cube(radius + thickness);
  else
    return M_4PI_3 * cube(radius + thickness);

}



static double fq(double q, double radius,
 double thickness, double core_sld, double shell_sld, double solvent_sld)
{
 const double form = core_shell_fq(q,
  radius,
  thickness,
  core_sld,
  shell_sld,
  solvent_sld);
 return form;
}


static double reduced_field(double q, double Ms, double Hi,
 double A)
{
    // q in 10e10 m-1, A in 10e-12 J/m, mu0 in 1e-7 
    return Ms / (fmax(Hi, 1.0e-6) + 2.0 * A * 4.0 * M_PI / Ms * q * q * 10.0);

 }

 static double DMI_length(double Ms, double D, double qval)
 {
   return 2.0 * D * 4.0 * M_PI / Ms / Ms * qval ; //q in 10e10 m-1, A in 10e-3 J/m^2, mu0 in 4 M_PI 1e-7 
 }

//Mz is defined as the longitudinal magnetisation component along the magnetic field.
//In the approach to saturation this component is (almost) constant with magnetic 
//field and simplfy reflects the nanoscale variations in the saturation magnetisation
//value in the sample. The misalignment of the magnetisation due to perturbing
//magnetic anisotropy or dipolar magnetic fields in the sample enter Mx and My,
//the two transversal magnetisation components, reacting to a magnetic field.
//The micromagnetic solution for the magnetisation are from Michels et al. PRB 94, 054424 (2016).

static double fqMxreal( double qx, double qy, double qz, double Mz, double Hkx, double Hky, double Hi, double Ms, double A, double D)
{
  double qsq = qx * qx + qy * qy + qz * qz;
  double q = sqrt(qsq); 
  double Hr = reduced_field(q, Ms, Hi, A);
  double DMI = DMI_length(Ms, D,q);
  double DMIz = DMI_length(Ms, D,qz);
  double denominator = (qsq + Hr * (qx * qx + qy * qy) - square(Hr * DMIz * q)) / Hr;
  double f = (Hkx * (qsq + Hr * qy * qy) - Ms * Mz * qx * qz * (1.0 + Hr * square(DMI)) - Hky * Hr * qx * qy) / denominator;
  return f;
}

static double fqMximag(double qx, double qy, double qz, double Mz, double Hkx, double Hky, double Hi, double Ms, double A, double D)
{
  double qsq = qx * qx + qy * qy + qz * qz;
  double q = sqrt(qsq);   
  double Hr = reduced_field(q, Ms, Hi, A);
  double DMI = DMI_length(Ms, D,q);
  double DMIz = DMI_length(Ms, D,qz);
  double DMIy = DMI_length(Ms, D,qy);
  double denominator = (qsq + Hr * (qx * qx + qy * qy) - square(Hr * DMIz * q)) / Hr;
  double f = - qsq * (Ms * Mz * (1.0 + Hr) * DMIy + Hky * Hr * DMIz) / denominator;
  return f;
}

static double fqMyreal( double qx, double qy, double qz, double Mz, double Hkx, double Hky, double Hi, double Ms, double A, double D)
{
  double qsq = qx * qx + qy * qy + qz * qz;
  double q = sqrt(qsq); 
  double Hr = reduced_field(q, Ms, Hi, A);
  double DMI = DMI_length(Ms, D,q);
  double DMIz = DMI_length(Ms, D,qz);
  double denominator = (qsq + Hr * (qx * qx + qy * qy) - square(Hr * DMIz * q))/Hr;
  double f = (Hky * (qsq + Hr * qx * qx) - Ms * Mz * qy * qz * (1.0 + Hr * square(DMI)) - Hkx * Hr * qx * qy)/denominator;
  return f;
}

static double fqMyimag( double qx, double qy, double qz, double Mz, double Hkx, double Hky, double Hi, double Ms, double A, double D)
{
  double qsq = qx * qx + qy * qy + qz * qz;
  double q = sqrt(qsq);  
  double Hr = reduced_field(q, Ms, Hi, A);
  double DMI = DMI_length(Ms, D,q);
  double DMIx = DMI_length(Ms, D,qx);
  double DMIz = DMI_length(Ms, D,qz);
  double denominator = (qsq + Hr * (qx * qx + qy * qy) - square(Hr * DMIz * q))/Hr;
  double f = qsq * (Ms * Mz * (1.0 + Hr) * DMIx - Hkx * Hr * DMIz)/denominator;
  return f;
}

static double
Calculate_Scattering(double q, double cos_theta, double sin_theta, double radius, double thickness, double nuc_sld_core, double nuc_sld_shell, double nuc_sld_solvent, double mag_sld_core, double mag_sld_shell, double mag_sld_solvent, double hk_sld_core, double Hi, double Ms, double A, double D,  double up_i, double up_f, double alpha, double beta)
{
    double qrot[3];
    set_scatvec(qrot,q,cos_theta, sin_theta, alpha, beta);
    // 0=dd.real, 1=dd.imag, 2=uu.real, 3=uu.imag,  4=du.real, 6=du.imag,  7=ud.real, 5=ud.imag
    double weights[8];
    set_weights(up_i, up_f, weights);
   
    double mz = fq(q, radius, thickness, mag_sld_core, mag_sld_shell, mag_sld_solvent);
    double nuc = fq(q, radius, thickness, nuc_sld_core, nuc_sld_shell, nuc_sld_solvent);
	double hk = fq(q, radius, 0, hk_sld_core, 0, 0);

    double cos_gamma, sin_gamma;
    double sld[8];
    //loop over random anisotropy axis with isotropic orientation gamma for Hkx and Hky
    //To be modified for textured material see also Weissmueller et al. PRB 63, 214414 (2001)
    double total_F2 = 0.0;
    for (int i = 0; i<GAUSS_N ;i++) {
      const double gamma = M_PI * (GAUSS_Z[i] + 1.0); // 0 .. 2 pi
      SINCOS(gamma, sin_gamma, cos_gamma);	
      //Only the core of the defect/particle in the matrix has an effective
      //anisotropy (for simplicity), for the effect of different, more complex
      // spatial profile of the anisotropy see Michels PRB 82, 024433 (2010)
      double Hkx = hk * sin_gamma;
      double Hky = hk * cos_gamma;

      double mxreal = fqMxreal(qrot[0], qrot[1], qrot[2], mz, Hkx, Hky, Hi, Ms, A, D);
      double mximag = fqMximag(qrot[0], qrot[1], qrot[2], mz, Hkx, Hky, Hi, Ms, A, D);
      double myreal = fqMyreal(qrot[0], qrot[1], qrot[2], mz, Hkx, Hky, Hi, Ms, A, D);
      double myimag = fqMyimag(qrot[0], qrot[1], qrot[2], mz, Hkx, Hky, Hi, Ms, A, D);

      mag_sld(qrot[0], qrot[1], qrot[2], mxreal, mximag, myreal, myimag, mz, 0, nuc, sld);
      double form = 0.0;
      for (unsigned int xs = 0; xs<8; xs++) {
        if (weights[xs] > 1.0e-8) {
          // Since the cross section weight is significant, set the slds
          // to the effective slds for this cross section, call the
          // kernel, and add according to weight.
          // loop over uu, ud real, du real, dd, ud imag, du imag 
          form += weights[xs] * sld[xs] * sld[xs];
        }
      }
      total_F2 += GAUSS_W[i] * form ;      
    }
    return total_F2;
}


static double
Iqxy(double qx, double qy, double radius, double thickness, double nuc_sld_core, double nuc_sld_shell, double nuc_sld_solvent, double mag_sld_core, double mag_sld_shell, double mag_sld_solvent, double hk_sld_core, double Hi, double Ms, double A, double D,  double up_i, double up_f, double alpha, double beta)
{
  double sin_theta, cos_theta;
  const double q = sqrt(qx * qx + qy * qy);  
  if (q > 1.0e-16 ) {
    cos_theta = qx/q;
    sin_theta = qy/q;
  } else {
    cos_theta = 0.0;
    sin_theta = 0.0;
  }  

  double total_F2 = Calculate_Scattering(q, cos_theta, sin_theta, radius, thickness, nuc_sld_core, nuc_sld_shell, nuc_sld_solvent, mag_sld_core, mag_sld_shell, mag_sld_solvent, hk_sld_core, Hi, Ms, A, D,  up_i, up_f, alpha, beta);
	
  return 0.5 * 1.0e-4 * total_F2;
    
}



static double
Iq(double q, double radius, double thickness, double nuc_sld_core, double nuc_sld_shell, double nuc_sld_solvent, double mag_sld_core, double mag_sld_shell, double mag_sld_solvent, double hk_sld_core, double Hi, double Ms, double A, double D,  double up_i, double up_f, double alpha, double beta)
{
  // slots to hold sincos function output of the orientation on the detector plane
  double sin_theta, cos_theta; 
  double total_F1D = 0.0;
  for (int j = 0; j<GAUSS_N ;j++) {

    const double theta = M_PI * (GAUSS_Z[j] + 1.0); // 0 .. 2 pi
    SINCOS(theta, sin_theta, cos_theta);

    double total_F2 = Calculate_Scattering(q, cos_theta, sin_theta, radius, thickness, nuc_sld_core, nuc_sld_shell, nuc_sld_solvent, mag_sld_core, mag_sld_shell, mag_sld_solvent, hk_sld_core, Hi, Ms, A, D,  up_i, up_f, alpha, beta);	
    total_F1D += GAUSS_W[j] * total_F2 ;
  }
  //convert from [1e-12 A-1] to [cm-1]
  return 0.25 * 1.0e-4 * total_F1D;
}






