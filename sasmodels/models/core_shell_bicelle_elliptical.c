double form_volume(double radius, double x_core, double thick_rim, double thick_face, double length);
double Iq(double q,
          double radius,
          double x_core,
          double thick_rim,
          double thick_face,
          double length,
          double core_sld,
          double face_sld,
          double rim_sld,
          double solvent_sld);


double Iqxy(double qx, double qy,
          double radius,
          double x_core,
          double thick_rim,
          double thick_face,
          double length,
          double core_sld,
          double face_sld,
          double rim_sld,
          double solvent_sld,
          double theta,
          double phi,
          double psi);

// NOTE that "length" here is the full height of the core!
double form_volume(double radius, double x_core, double thick_rim, double thick_face, double length)
{
    return M_PI*(radius+thick_rim)*(radius*x_core+thick_rim)*(length+2.0*thick_face);
}

double Iq(double q,
          double rad,
          double x_core,
          double radthick,
          double facthick,
          double length,
          double rhoc,
          double rhoh,
          double rhor,
          double rhosolv)
{
    double si1,si2,be1,be2;
     // core_shell_bicelle_elliptical, RKH Dec 2016, based on elliptical_cylinder and core_shell_bicelle
     // tested against limiting cases of cylinder, elliptical_cylinder and core_shell_bicelle
     //    const double uplim = M_PI_4;
    const double halfheight = 0.5*length;
    //const double va = 0.0;
    //const double vb = 1.0;
    // inner integral limits
    //const double vaj=0.0;
    //const double vbj=M_PI;

    const double radius_major = rad * x_core;
    const double rA = 0.5*(square(radius_major) + square(rad));
    const double rB = 0.5*(square(radius_major) - square(rad));
    const double dr1 = (rhoc-rhoh)   *M_PI*rad*radius_major*(2.0*halfheight);;
    const double dr2 = (rhor-rhosolv)*M_PI*(rad+radthick)*(radius_major+radthick)*2.0*(halfheight+facthick);
    const double dr3 = (rhoh-rhor)   *M_PI*rad*radius_major*2.0*(halfheight+facthick);
    //const double vol1 = M_PI*rad*radius_major*(2.0*halfheight);
    //const double vol2 = M_PI*(rad+radthick)*(radius_major+radthick)*2.0*(halfheight+facthick);
    //const double vol3 = M_PI*rad*radius_major*2.0*(halfheight+facthick);

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
        double sinarg2 = q*(halfheight+facthick)*cos_alpha;
        si1 = sas_sinx_x(sinarg1);
        si2 = sas_sinx_x(sinarg2);
        for(int j=0;j<76;j++) {
            //76 gauss points for the inner integral (WAS 20 points,so this may make unecessarily slow, but playing safe)
            //const double beta = ( Gauss76Z[j]*(vbj-vaj) + vaj + vbj )/2.0;
            const double beta = ( Gauss76Z[j] +1.0)*M_PI_2;
            const double rr = sqrt(rA - rB*cos(beta));
            double besarg1 = q*rr*sin_alpha;
            double besarg2 = q*(rr+radthick)*sin_alpha;
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

double 
Iqxy(double qx, double qy,
          double rad,
          double x_core,
          double radthick,
          double facthick,
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
    const double radius_major = rad*x_core;
    const double halfheight = 0.5*length;
    const double vol1 = M_PI*rad*radius_major*length;
    const double vol2 = M_PI*(rad+radthick)*(radius_major+radthick)*2.0*(halfheight+facthick);
    const double vol3 = M_PI*rad*radius_major*2.0*(halfheight+facthick);

    // Compute:  r = sqrt((radius_major*zhat)^2 + (radius_minor*yhat)^2)
    // Given:    radius_major = r_ratio * radius_minor  
    // ASSUME the sin_alpha is included in the separate integration over orientation of rod angle
    const double rad_minor = rad;
    const double rad_major = rad*x_core;
    const double r_hat = sqrt(square(rad_major*xhat) + square(rad_minor*yhat));
    const double rshell_hat = sqrt(square((rad_major+radthick)*xhat)
                                   + square((rad_minor+radthick)*yhat));
    const double be1 = sas_2J1x_x( q*r_hat );
    const double be2 = sas_2J1x_x( q*rshell_hat );
    const double si1 = sas_sinx_x( q*halfheight*zhat );
    const double si2 = sas_sinx_x( q*(halfheight + facthick)*zhat );
    const double Aq = square( vol1*dr1*si1*be1 + vol2*dr2*si2*be2 +  vol3*dr3*si2*be1);
    //const double vol = form_volume(radius_minor, r_ratio, length);
    return 1.0e-4 * Aq;
}

