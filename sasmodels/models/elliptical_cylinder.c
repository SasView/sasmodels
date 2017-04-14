static double
form_volume(double radius_minor, double r_ratio, double length)
{
    return M_PI * radius_minor * radius_minor * r_ratio * length;
}

static double
Iq(double q, double radius_minor, double r_ratio, double length,
   double sld, double solvent_sld)
{
    // orientational average limits
    const double va = 0.0;
    const double vb = 1.0;
    // inner integral limits
    const double vaj=0.0;
    const double vbj=M_PI;

    const double radius_major = r_ratio * radius_minor;
    const double rA = 0.5*(square(radius_major) + square(radius_minor));
    const double rB = 0.5*(square(radius_major) - square(radius_minor));

    //initialize integral
    double outer_sum = 0.0;
    for(int i=0;i<76;i++) {
        //setup inner integral over the ellipsoidal cross-section
        const double cos_val = ( Gauss76Z[i]*(vb-va) + va + vb )/2.0;
        const double sin_val = sqrt(1.0 - cos_val*cos_val);
        //const double arg = radius_minor*sin_val;
        double inner_sum=0;
        for(int j=0;j<20;j++) {
            //20 gauss points for the inner integral
            const double theta = ( Gauss20Z[j]*(vbj-vaj) + vaj + vbj )/2.0;
            const double r = sin_val*sqrt(rA - rB*cos(theta));
            const double be = sas_2J1x_x(q*r);
            inner_sum += Gauss20Wt[j] * be * be;
        }
        //now calculate the value of the inner integral
        inner_sum *= 0.5*(vbj-vaj);

        //now calculate outer integral
        const double si = sas_sinx_x(q*0.5*length*cos_val);
        outer_sum += Gauss76Wt[i] * inner_sum * si * si;
    }
    outer_sum *= 0.5*(vb-va);

    //divide integral by Pi
    const double form = outer_sum/M_PI;

    // scale by contrast and volume, and convert to to 1/cm units
    const double vol = form_volume(radius_minor, r_ratio, length);
    const double delrho = sld - solvent_sld;
    return 1.0e-4*square(delrho*vol)*form;
}


static double
Iqxy(double qx, double qy,
     double radius_minor, double r_ratio, double length,
     double sld, double solvent_sld,
     double theta, double phi, double psi)
{
    double q, xhat, yhat, zhat;
    ORIENT_ASYMMETRIC(qx, qy, theta, phi, psi, q, xhat, yhat, zhat);
    const double qa = q*xhat;
    const double qb = q*yhat;
    const double qc = q*zhat;

    // Compute:  r = sqrt((radius_major*cos_nu)^2 + (radius_minor*cos_mu)^2)
    // Given:    radius_major = r_ratio * radius_minor
    const double qr = radius_minor*sqrt(square(r_ratio*qa) + square(qb));
    const double be = sas_2J1x_x(qr);
    const double si = sas_sinx_x(qc*0.5*length);
    const double Aq = be * si;
    const double delrho = sld - solvent_sld;
    const double vol = form_volume(radius_minor, r_ratio, length);
    return 1.0e-4 * square(delrho * vol * Aq);
}
