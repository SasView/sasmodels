static double
bcc_Zq(double qa, double qb, double qc, double dnn, double d_factor)
{
#if 0  // Equations as written in Matsuoka
    const double a1 = (+qa + qb + qc)/2.0;
    const double a2 = (-qa - qb + qc)/2.0;
    const double a3 = (-qa + qb - qc)/2.0;
#else
    const double a1 = (+qa + qb - qc)/2.0;
    const double a2 = (+qa - qb + qc)/2.0;
    const double a3 = (-qa + qb + qc)/2.0;
#endif

#if 1
    // Numerator: (1 - exp(a)^2)^3
    //         => (-(exp(2a) - 1))^3
    //         => -expm1(2a)^3
    // Denominator: prod(1 - 2 cos(d ak) exp(a) + exp(2a))
    //         => prod(exp(a)^2 - 2 cos(d ak) exp(a) + 1)
    //         => prod((exp(a) - 2 cos(d ak)) * exp(a) + 1)
    const double arg = -0.5*square(dnn*d_factor)*(a1*a1 + a2*a2 + a3*a3);
    const double exp_arg = exp(arg);
    const double Zq = -cube(expm1(2.0*arg))
        / ( ((exp_arg - 2.0*cos(dnn*a1))*exp_arg + 1.0)
          * ((exp_arg - 2.0*cos(dnn*a2))*exp_arg + 1.0)
          * ((exp_arg - 2.0*cos(dnn*a3))*exp_arg + 1.0));
#else
    // Alternate form, which perhaps is more approachable
    const double arg = -0.5*square(dnn*d_factor)*(a1*a1 + a2*a2 + a3*a3);
    const double sinh_qd = sinh(arg);
    const double cosh_qd = cosh(arg);
    const double Zq = sinh_qd/(cosh_qd - cos(dnn*a1))
                    * sinh_qd/(cosh_qd - cos(dnn*a2))
                    * sinh_qd/(cosh_qd - cos(dnn*a3));
#endif

    return Zq;
}


// occupied volume fraction calculated from lattice symmetry and sphere radius
static double
bcc_volume_fraction(double radius, double dnn)
{
    return 2.0*sphere_volume(sqrt(0.75)*radius/dnn);
}

static double
form_volume(double radius)
{
    return sphere_volume(radius);
}


static double Iq(double q, double dnn,
    double d_factor, double radius,
    double sld, double solvent_sld)
{
    // translate a point in [-1,1] to a point in [0, 2 pi]
    const double phi_m = M_PI;
    const double phi_b = M_PI;
    // translate a point in [-1,1] to a point in [0, pi]
    const double theta_m = M_PI_2;
    const double theta_b = M_PI_2;

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
            const double form = bcc_Zq(qa, qb, qc, dnn, d_factor);
            inner_sum += Gauss150Wt[j] * form;
        }
        inner_sum *= phi_m;  // sum(f(x)dx) = sum(f(x)) dx
        outer_sum += Gauss150Wt[i] * inner_sum * sin_theta;
    }
    outer_sum *= theta_m;
    const double Zq = outer_sum/(4.0*M_PI);
    const double Pq = sphere_form(q, radius, sld, solvent_sld);
    return bcc_volume_fraction(radius, dnn) * Pq * Zq;
}


static double Iqxy(double qa, double qb, double qc,
    double dnn, double d_factor, double radius,
    double sld, double solvent_sld)
{
    const double q = sqrt(qa*qa + qb*qb + qc*qc);
    const double Zq = bcc_Zq(qa, qb, qc, dnn, d_factor);
    const double Pq = sphere_form(q, radius, sld, solvent_sld);
    return bcc_volume_fraction(radius, dnn) * Pq * Zq;
}