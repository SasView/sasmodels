double form_volume(double length_a, double length_b, double length_c, 
                   double thick_rim_a, double thick_rim_b, double thick_rim_c);
double Iq(double q, double core_sld, double arim_sld, double brim_sld, double crim_sld,
          double solvent_sld, double length_a, double length_b, double length_c,
          double thick_rim_a, double thick_rim_b, double thick_rim_c);
double Iqxy(double qx, double qy, double core_sld, double arim_sld, double brim_sld,
            double crim_sld, double solvent_sld, double length_a, double length_b,
            double length_c, double thick_rim_a, double thick_rim_b,
            double thick_rim_c, double theta, double phi, double psi);

double form_volume(double length_a, double length_b, double length_c, 
                   double thick_rim_a, double thick_rim_b, double thick_rim_c)
{
    //return length_a * length_b * length_c;
    return length_a * length_b * length_c + 
           2.0 * thick_rim_a * length_b * length_c + 
           2.0 * thick_rim_b * length_a * length_c +
           2.0 * thick_rim_c * length_a * length_b;
}

double Iq(double q,
    double core_sld,
    double arim_sld,
    double brim_sld,
    double crim_sld,
    double solvent_sld,
    double length_a,
    double length_b,
    double length_c,
    double thick_rim_a,
    double thick_rim_b,
    double thick_rim_c)
{
    // Code converted from functions CSPPKernel and CSParallelepiped in libCylinder.c_scaled
    // Did not understand the code completely, it should be rechecked (Miguel Gonzalez)
    
    const double mu = 0.5 * q * length_b;
    
    //calculate volume before rescaling (in original code, but not used)
    //double vol = form_volume(length_a, length_b, length_c, thick_rim_a, thick_rim_b, thick_rim_c);        
    //double vol = length_a * length_b * length_c + 
    //       2.0 * thick_rim_a * length_b * length_c + 
    //       2.0 * thick_rim_b * length_a * length_c +
    //       2.0 * thick_rim_c * length_a * length_b;
    
    // Scale sides by B
    const double a_scaled = length_a / length_b;
    const double c_scaled = length_c / length_b;

    // ta and tb correspond to the definitions in CSPPKernel, but they don't make sense to me (MG)
    // the a_scaled in the definition of tb was present in CSPPKernel in libCylinder.c,
    // while in cspkernel in csparallelepiped.cpp (used for the 2D), all the definitions
    // for ta, tb, tc use also A + 2*rim_thickness (but not scaled by B!!!)
    double ta = (a_scaled + 2.0*thick_rim_a)/length_b;
    double tb = (a_scaled + 2.0*thick_rim_b)/length_b;

    double Vin = length_a * length_b * length_c;
    //double Vot = (length_a * length_b * length_c +
    //            2.0 * thick_rim_a * length_b * length_c +
    //            2.0 * length_a * thick_rim_b * length_c +
    //            2.0 * length_a * length_b * thick_rim_c);
    double V1 = (2.0 * thick_rim_a * length_b * length_c);    // incorrect V1 (aa*bb*cc+2*ta*bb*cc)
    double V2 = (2.0 * length_a * thick_rim_b * length_c);    // incorrect V2(aa*bb*cc+2*aa*tb*cc)

    // Scale factors (note that drC is not used later)
    const double drho0 = (core_sld-solvent_sld);
    const double drhoA = (arim_sld-solvent_sld);
    const double drhoB = (brim_sld-solvent_sld);
    //const double drC_Vot = (crim_sld-solvent_sld)*Vot;

    // Precompute scale factors for combining cross terms from the shape
    const double scale23 = drhoA*V1;
    const double scale14 = drhoB*V2;
    const double scale12 = drho0*Vin - scale23 - scale14;

    // outer integral (with gauss points), integration limits = 0, 1
    double outer_total = 0; //initialize integral

    for( int i=0; i<76; i++) {
        double sigma = 0.5 * ( Gauss76Z[i] + 1.0 );
        double mu_proj = mu * sqrt(1.0-sigma*sigma);

        // inner integral (with gauss points), integration limits = 0, 1
        double inner_total = 0.0;
        for(int j=0; j<76; j++) {
            const double uu = 0.5 * ( Gauss76Z[j] + 1.0 );
            double sin_uu, cos_uu;
            SINCOS(M_PI_2*uu, sin_uu, cos_uu);
            const double si1 = sinc(mu_proj * sin_uu * a_scaled);
            const double si2 = sinc(mu_proj * cos_uu);
            const double si3 = sinc(mu_proj * sin_uu * ta);
            const double si4 = sinc(mu_proj * cos_uu * tb);

            // Expression in libCylinder.c (neither drC nor Vot are used)
            const double form = scale12*si1*si2 + scale23*si2*si3 + scale14*si1*si4;

            // To note also that in csparallelepiped.cpp, there is a function called
            // cspkernel, which seems to make more sense and has the following comment:
            //   This expression is different from NIST/IGOR package (I strongly believe the IGOR is wrong!!!). 10/15/2010.
            //   tmp =( dr0*tmp1*tmp2*tmp3*Vin + drA*(tmpt1-tmp1)*tmp2*tmp3*V1+ drB*tmp1*(tmpt2-tmp2)*tmp3*V2 + drC*tmp1*tmp2*(tmpt3-tmp3)*V3)*
            //   ( dr0*tmp1*tmp2*tmp3*Vin + drA*(tmpt1-tmp1)*tmp2*tmp3*V1+ drB*tmp1*(tmpt2-tmp2)*tmp3*V2 + drC*tmp1*tmp2*(tmpt3-tmp3)*V3);   //  correct FF : square of sum of phase factors
            // This is the function called by csparallelepiped_analytical_2D_scaled,
            // while CSParallelepipedModel calls CSParallelepiped in libCylinder.c        
            
            //  correct FF : sum of square of phase factors
            inner_total += Gauss76Wt[j] * form * form;
        }
        inner_total *= 0.5;

        // now sum up the outer integral
        const double si = sinc(mu * c_scaled * sigma);
        outer_total += Gauss76Wt[i] * inner_total * si * si;
    }
    outer_total *= 0.5;

    //convert from [1e-12 A-1] to [cm-1]
    return 1.0e-4 * outer_total;
}

double Iqxy(double qx, double qy,
    double core_sld,
    double arim_sld,
    double brim_sld,
    double crim_sld,
    double solvent_sld,
    double length_a,
    double length_b,
    double length_c,
    double thick_rim_a,
    double thick_rim_b,
    double thick_rim_c,
    double theta,
    double phi,
    double psi)
{
    double q, cos_val_a, cos_val_b, cos_val_c;
    ORIENT_ASYMMETRIC(qx, qy, theta, phi, psi, q, cos_val_c, cos_val_b, cos_val_a);

    // cspkernel in csparallelepiped recoded here
    const double dr0 = core_sld-solvent_sld;
    const double drA = arim_sld-solvent_sld;
    const double drB = brim_sld-solvent_sld;
    const double drC = crim_sld-solvent_sld;

    double Vin = length_a * length_b * length_c;
    double V1 = 2.0 * thick_rim_a * length_b * length_c;    // incorrect V1(aa*bb*cc+2*ta*bb*cc)
    double V2 = 2.0 * length_a * thick_rim_b * length_c;    // incorrect V2(aa*bb*cc+2*aa*tb*cc)
    double V3 = 2.0 * length_a * length_b * thick_rim_c;
    // As for the 1D case, Vot is not used
    //double Vot = (length_a * length_b * length_c +
    //              2.0 * thick_rim_a * length_b * length_c +
    //              2.0 * length_a * thick_rim_b * length_c +
    //              2.0 * length_a * length_b * thick_rim_c);

    // The definitions of ta, tb, tc are not the same as in the 1D case because there is no
    // the scaling by B. The use of length_a for the 3 of them seems clearly a mistake to me,
    // but for the moment I let it like this until understanding better the code.
    double ta = length_a + 2.0*thick_rim_a;
    double tb = length_a + 2.0*thick_rim_b;
    double tc = length_a + 2.0*thick_rim_c;
    //handle arg=0 separately, as sin(t)/t -> 1 as t->0
    double siA = sinc(0.5*q*length_a*cos_val_a);
    double siB = sinc(0.5*q*length_b*cos_val_b);
    double siC = sinc(0.5*q*length_c*cos_val_c);
    double siAt = sinc(0.5*q*ta*cos_val_a);
    double siBt = sinc(0.5*q*tb*cos_val_b);
    double siCt = sinc(0.5*q*tc*cos_val_c);
    

    // f uses Vin, V1, V2, and V3 and it seems to have more sense than the value computed
    // in the 1D code, but should be checked!
    double f = ( dr0*siA*siB*siC*Vin
               + drA*(siAt-siA)*siB*siC*V1
               + drB*siA*(siBt-siB)*siC*V2
               + drC*siA*siB*(siCt*siCt-siC)*V3);
   
    return 1.0e-4 * f * f;
}
