// vd = volume * delta_rho
// besarg = q * R * sin(theta)
// siarg = q * L/2 * cos(theta)
static double _cyl(double vd, double besarg, double siarg)
{
    return vd * sas_sinx_x(siarg) * sas_2J1x_x(besarg);
}

static double
form_volume(double radius_core, double length_core, 
        double n_shells, double thick[], double face[] )
        {      
// this one gets used by rest of sasview, it needs all the "volume" parameters
    double r_out = radius_core;
    double h_out = 0.5*length_core;
    for (int j=0; j < n_shells; j++){
         r_out += thick[j];
         h_out += face[j];
    }
    return M_PI*square(r_out)*2.0*h_out;
}

static double
form_volume_here(double radius, double halflength)
// a quicker local one for these routines
{
    return M_PI*square(radius)*2.0*halflength;
}

static double
radius_from_excluded_volume(double radius_tot, double halflength)

{
    return 0.5*cbrt(0.75*radius_tot*(4.0*radius_tot*halflength + 
    (radius_tot + 2.0*halflength)*(M_PI*radius_tot + 2.0*halflength)));
}

static double
radius_from_volume(double radius, double halflength)
{
    return cbrt(form_volume_here(radius,halflength)/M_4PI_3);
}

static double
radius_from_diagonal(double radius_outer, double halflength)
{
     return sqrt(radius_outer*radius_outer + halflength*halflength);
}

static double
radius_effective(int mode, double radius_core, double length_core, 
        double n_shells, double thick[], double face[] )
{
    double r_out = radius_core;
    double h_out = 0.5*length_core;
    for (int j=0; j < n_shells; j++){
         r_out += thick[j];
         h_out += face[j];
    }  
    switch (mode) {
    default:
    case 1: //cylinder excluded volume
        return radius_from_excluded_volume(r_out, h_out);
    case 2: // equivalent volume sphere
        return radius_from_volume(r_out, h_out);
    case 3: // outer radius
        return r_out;
    case 4: // half outer length
        return h_out;
    case 5: // half min outer length
        return (r_out < h_out ? r_out : h_out );
    case 6: // half max outer length
        return (r_out > h_out ? r_out : h_out );
    case 7: // half outer diagonal
        return radius_from_diagonal(r_out, h_out);
    }
}

static void
Fq(double q, double *F1, double *F2, double sld_core, double radius_core, double length_core, 
    double sld_solvent, double sigma,
    double n_shells, double sld[], double thick[], double face[] )
{
#define MAX_SHELLS_PLUS_2 12
    // the sum over shells in the contrast has to be done inside the integration over angle
    // this is messier than onion.c where for spheres the form factor is analytic
    // NOTE assume at least one shell
    double sin_theta, cos_theta;
    double total_F1 = 0.0;
    double total_F2 = 0.0;
    // avoid repeated calculations and create single loop
    // might be able to do double r_out[n_shells]    but Paul K wants to keep multiples of 4 
    // for some possible gpu issues. Since the model's .py file defines the max limit as 10 
    // shells, need 2 more here (for core & solvent) which nicely makes 12
    double r_out[MAX_SHELLS_PLUS_2], h_out[MAX_SHELLS_PLUS_2], vd[MAX_SHELLS_PLUS_2];
    int n = (int)(n_shells+0.5);
    r_out[0] = radius_core;
    // h_out is half length of core or current shell, r_out is its radius
    h_out[0] = length_core*0.5;
    // first term, for core, needs refer to sld_core and first shell (hence restict model to n>= 1)
    vd[0] = form_volume_here(r_out[0], h_out[0]) * ( sld_core - sld[0] );
    // loop over first to next to last shells, only do this for n >= 2 shells
    // n shells, has n+1 terms, indexed as 0 to n 
    for (int j=0; j < n-1; j++){
         r_out[j+1] = r_out[j] + thick[j];
         h_out[j+1] = h_out[j] + face[j];
         vd[j+1] = form_volume_here(r_out[j+1], h_out[j+1]) * ( sld[j] - sld[j+1] );
         }
    // now do the last shell, which needs refer to solvent_sld
    r_out[n] = r_out[n-1] + thick[n-1];
    h_out[n] = h_out[n-1] + face[n-1];
    vd[n] = form_volume_here(r_out[n], h_out[n]) * ( sld[n-1] - sld_solvent );
    //
    // now we are ready to do the longer loop of the quadrature integration
    //printf("n r h  = %d %g %g\n", n, r_out[n], h_out[n]);
    for (int i=0; i<GAUSS_N ;i++) {
       // translate a point in [-1,1] to a point in [0, pi/2]
       //const double theta = ( GAUSS_Z[i]*(upper-lower) + upper + lower )/2.0;
       const double theta = GAUSS_Z[i]*M_PI_4 + M_PI_4;
       SINCOS(theta, sin_theta,  cos_theta);
       const double qab = q*sin_theta;
       const double qc = q*cos_theta;
       double f =0.0;
       // sum over core and shells, terms  j = 0 to n
       for (int j=0; j < n+1; j++){
         f +=  _cyl( vd[j], r_out[j]*qab, h_out[j]*qc );
         }
       total_F1 += GAUSS_W[i] * f * sin_theta;
       total_F2 += GAUSS_W[i] * f * f * sin_theta;
     }
    // end of loops 
    //  approximate interfacial roughness damping    
    const double atten = exp(-0.25*square(q*sigma));
    *F1 = 1.0e-2 * total_F1 * M_PI_4*atten;
    *F2 = 1.0e-4 * total_F2 * M_PI_4*square(atten);

}

static double
Iqac(double qab, double qc,
    double sld_core, double radius_core, double length_core, 
    double sld_solvent, double sigma,
    double n_shells, double sld[], double thick[], double face[] )
                                         
{
       int n = (int)(n_shells+0.5);
       double r_out = radius_core;
       // h_out is half length of core or current shell, r_out is its radius
       double h_out = length_core*0.5;
       // first term, for core, needs refer to core_sld and first shell (hence restict model to n>= 1)
       double vd = form_volume_here(r_out, h_out) * ( sld_core - sld[0]);
       double f = _cyl(vd, r_out*qab, h_out*qc);
       // loop over first to next to last shells, for n >= 2
       for (int j=0; j < n-1; j++){
         r_out += thick[j];
         h_out += face[j];
         vd = form_volume_here(r_out, h_out) * ( sld[j] - sld[j+1] );
         f +=  _cyl(vd, r_out*qab, h_out*qc);
       }
       // now do the last shell, needs refer to solvent_sld
       r_out += thick[n-1];
       h_out += face[n-1];
       vd = form_volume_here(r_out, h_out) * ( sld[n-1]-sld_solvent );
       f +=  _cyl(vd, r_out*qab, h_out*qc);
     //  approximate interfacial roughness damping
    const double atten = exp(-0.5*(qab*qab + qc*qc)*(sigma*sigma));
    return 1.0e-4 * f*f * atten;        
}

