// NOTE that "length" here is the full height of the core!
static double
form_volume(double r_minor,
        double x_core,
        double thick_rim,
        double thick_face,
        double length)
{
    return M_PI*(r_minor+thick_rim)*(r_minor*x_core+thick_rim)*(length+2.0*thick_face);
}

static double
Iq(double q,
        double r_minor,
        double x_core,
        double thick_rim,
        double thick_face,
        double length,
        double rhoc,
        double rhoh,
        double rhor,
        double rhosolv)
{
    double si1,si2,be1,be2;
     // core_shell_bicelle_elliptical, RKH Dec 2016, based on elliptical_cylinder and core_shell_bicelle
     // tested against limiting cases of cylinder, elliptical_cylinder, stacked_discs, and core_shell_bicelle
     //    const double uplim = M_PI_4;
    const double halfheight = 0.5*length;
    //const double va = 0.0;
    //const double vb = 1.0;
    // inner integral limits
    //const double vaj=0.0;
    //const double vbj=M_PI;

    const double r_major = r_minor * x_core;
    const double r2A = 0.5*(square(r_major) + square(r_minor));
    const double r2B = 0.5*(square(r_major) - square(r_minor));
    const double dr1 = (rhoc-rhoh)   *M_PI*r_minor*r_major*(2.0*halfheight);;
    const double dr2 = (rhor-rhosolv)*M_PI*(r_minor+thick_rim)*(r_major+thick_rim)*2.0*(halfheight+thick_face);
    const double dr3 = (rhoh-rhor)   *M_PI*r_minor*r_major*2.0*(halfheight+thick_face);
    //const double vol1 = M_PI*r_minor*r_major*(2.0*halfheight);
    //const double vol2 = M_PI*(r_minor+thick_rim)*(r_major+thick_rim)*2.0*(halfheight+thick_face);
    //const double vol3 = M_PI*r_minor*r_major*2.0*(halfheight+thick_face);

    //initialize integral
    double outer_sum = 0.0;
    for(int i=0;i<76;i++) {
        //setup inner integral over the ellipsoidal cross-section
        // since we generate these lots of times, why not store them somewhere?
        //const double cos_alpha = ( Gauss76Z[i]*(vb-va) + va + vb )/2.0;
        const double cos_alpha = ( Gauss76Z[i] + 1.0 )/2.0;
        const double sin_alpha = sqrt(1.0 - cos_alpha*cos_alpha);
        double inner_sum=0;
        double sinarg1 = q*halfheight*cos_alpha;
        double sinarg2 = q*(halfheight+thick_face)*cos_alpha;
        si1 = sas_sinx_x(sinarg1);
        si2 = sas_sinx_x(sinarg2);
        for(int j=0;j<76;j++) {
            //76 gauss points for the inner integral (WAS 20 points,so this may make unecessarily slow, but playing safe)
            //const double beta = ( Gauss76Z[j]*(vbj-vaj) + vaj + vbj )/2.0;
            const double beta = ( Gauss76Z[j] +1.0)*M_PI_2;
            const double rr = sqrt(r2A - r2B*cos(beta));
            double besarg1 = q*rr*sin_alpha;
            double besarg2 = q*(rr+thick_rim)*sin_alpha;
            be1 = sas_2J1x_x(besarg1);
            be2 = sas_2J1x_x(besarg2);
            inner_sum += Gauss76Wt[j] *square(dr1*si1*be1 +
                                              dr2*si2*be2 +
                                              dr3*si2*be1);
        }
        //now calculate outer integral
        outer_sum += Gauss76Wt[i] * inner_sum;
    }

    return outer_sum*2.5e-05;
}

static double
Iqxy(double qx, double qy,
          double r_minor,
          double x_core,
          double thick_rim,
          double thick_face,
          double length,
          double rhoc,
          double rhoh,
          double rhor,
          double rhosolv,
          double theta,
          double phi,
          double psi)
{
       // THIS NEEDS TESTING
    double q, xhat, yhat, zhat;
    ORIENT_ASYMMETRIC(qx, qy, theta, phi, psi, q, xhat, yhat, zhat);
    const double dr1 = rhoc-rhoh;
    const double dr2 = rhor-rhosolv;
    const double dr3 = rhoh-rhor;
    const double r_major = r_minor*x_core;
    const double halfheight = 0.5*length;
    const double vol1 = M_PI*r_minor*r_major*length;
    const double vol2 = M_PI*(r_minor+thick_rim)*(r_major+thick_rim)*2.0*(halfheight+thick_face);
    const double vol3 = M_PI*r_minor*r_major*2.0*(halfheight+thick_face);

    // Compute effective radius in rotated coordinates
    const double r_hat = sqrt(square(r_major*xhat) + square(r_minor*yhat));
    const double rshell_hat = sqrt(square((r_major+thick_rim)*xhat)
                                   + square((r_minor+thick_rim)*yhat));
    const double be1 = sas_2J1x_x( q*r_hat );
    const double be2 = sas_2J1x_x( q*rshell_hat );
    const double si1 = sas_sinx_x( q*halfheight*zhat );
    const double si2 = sas_sinx_x( q*(halfheight + thick_face)*zhat );
    const double Aq = square( vol1*dr1*si1*be1 + vol2*dr2*si2*be2 +  vol3*dr3*si2*be1);
    return 1.0e-4 * Aq;
}

