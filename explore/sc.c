static double
sc_Zq(double qa, double qb, double qc, double dnn, double d_factor)
{
    // Rewriting equations for efficiency, accuracy and readability, and so
    // code is reusable between 1D and 2D models.
    #if 1  // SC
    const double a1 = qa;
    const double a2 = qb;
    const double a3 = qc;
    #elif 1 // BCC
    const double a1 = (+qa + qb + qc)/2.;
    const double a2 = (-qa - qb + qc)/2.;
    const double a3 = (-qa + qb - qc)/2.;
    #elif 1 // FCC
    const double a1 = ( qa + qb)/2.0;
    const double a2 = (-qa + qc)/2.0;
    const double a3 = (-qa + qb)/2.0;
    #endif

    const double arg = -0.5*square(dnn*d_factor)*(a1*a1 + a2*a2 + a3*a3);

    // Numerator: (1 - exp(a)^2)^3
    //         => (-(exp(2a) - 1))^3
    //         => -expm1(2a)^3
    // Denominator: prod(1 - 2 cos(xk) exp(a) + exp(a)^2)
    //         => exp(a)^2 - 2 cos(xk) exp(a) + 1
    //         => (exp(a) - 2 cos(xk)) * exp(a) + 1
    const double exp_arg = exp(arg);
    const double Zq = -cube(expm1(2.0*arg))
        / ( ((exp_arg - 2.0*cos(dnn*a1))*exp_arg + 1.0)
          * ((exp_arg - 2.0*cos(dnn*a2))*exp_arg + 1.0)
          * ((exp_arg - 2.0*cos(dnn*a3))*exp_arg + 1.0));

    return Zq;
}

// occupied volume fraction calculated from lattice symmetry and sphere radius
static double
sc_volume_fraction(double radius, double dnn)
{
    return sphere_volume(radius/dnn);
}

static double
form_volume(double radius)
{
    return sphere_volume(radius);
}


static double Iq(double q, double dnn,
  double d_factor, double radius,
  double sld, double solvent_sld,
  double n, double sym)
{
double phi_m, phi_b, theta_m, theta_b;
if (sym>0.) {
    // translate a point in [-1,1] to a point in [0, 2 pi]
    phi_m = M_PI_4;
    phi_b = M_PI_4;
    // translate a point in [-1,1] to a point in [0, pi]
    theta_m = M_PI_4;
    theta_b = M_PI_4;
} else {
    // translate a point in [-1,1] to a point in [0, 2 pi]
    phi_m = M_PI;
    phi_b = M_PI;
    // translate a point in [-1,1] to a point in [0, pi]
    theta_m = M_PI_2;
    theta_b = M_PI_2;
}

#if 0
    double outer_sum = 0.0;
    for(int i=0; i<150; i++) {
        double inner_sum = 0.0;
        const double theta = Gauss150Z[i]*theta_m + theta_b;
        double sin_theta, cos_theta;
        SINCOS(theta, sin_theta, cos_theta);
        const double qc = q*cos_theta;
        const double qab = q*sin_theta;
        for(int j=0;j<150;j++) {
            const double phi = Gauss150Z[j]*phi_m + phi_b;
            double sin_phi, cos_phi;
            SINCOS(phi, sin_phi, cos_phi);
            const double qa = qab*cos_phi;
            const double qb = qab*sin_phi;
            const double fq = _sq_sc(qa, qb, qc, dnn, d_factor);
            inner_sum += Gauss150Wt[j] * fq;
        }
        inner_sum *= phi_m;  // sum(f(x)dx) = sum(f(x)) dx
        outer_sum += Gauss150Wt[i] * inner_sum * sin_theta;
    }
    outer_sum *= theta_m;
#else
    double outer_sum = 0.0;
    for(int i=0; i<(int)n; i++) {
        double inner_sum = 0.0;
        const double theta = (i*2./n-1.)*theta_m + theta_b;
        double sin_theta, cos_theta;
        SINCOS(theta, sin_theta, cos_theta);
        const double qc = q*cos_theta;
        const double qab = q*sin_theta;
        for(int j=0;j<(int)n;j++) {
            const double phi = (j*2./n-1.)*phi_m + phi_b;
            double sin_phi, cos_phi;
            SINCOS(phi, sin_phi, cos_phi);
            const double qa = qab*cos_phi;
            const double qb = qab*sin_phi;
            const double form = sc_Zq(qa, qb, qc, dnn, d_factor);
            inner_sum += form;
        }
        inner_sum *= phi_m;  // sum(f(x)dx) = sum(f(x)) dx
        outer_sum += inner_sum * sin_theta;
    }
    outer_sum *= theta_m/(n*n);
#endif
double Zq;
if (sym > 0.) {
    Zq = outer_sum/M_PI_2;
} else {
    Zq = outer_sum/(4.0*M_PI);
}

    //return Zq;
    const double Pq = sphere_form(q, radius, sld, solvent_sld);
    return sc_volume_fraction(radius, dnn) * Pq * Zq;
}


static double Iqxy(double qx, double qy,
    double dnn, double d_factor, double radius,
    double sld, double solvent_sld,
    double n, double sym,
    double theta, double phi, double psi)
{
    double q, zhat, yhat, xhat;
    ORIENT_ASYMMETRIC(qx, qy, theta, phi, psi, q, xhat, yhat, zhat);
    const double qa = q*xhat;
    const double qb = q*yhat;
    const double qc = q*zhat;

    q = sqrt(qa*qa + qb*qb + qc*qc);
    const double Pq = sphere_form(q, radius, sld, solvent_sld);
    const double Zq = sc_Zq(qa, qb, qc, dnn, d_factor);
    return Zq;
    return sc_volume_fraction(radius, dnn) * Pq * Zq;
}
