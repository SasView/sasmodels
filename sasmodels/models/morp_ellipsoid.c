static double
form_volume(double radius_polar, double radius_equatorial)
{
    return M_4PI_3*radius_polar*radius_equatorial*radius_equatorial;
}

static double
radius_from_volume(double radius_polar, double radius_equatorial)
{
    return cbrt(radius_polar*radius_equatorial*radius_equatorial);
}

static double
radius_from_curvature(double radius_polar, double radius_equatorial)
{
    // Trivial cases
    
    if (radius_polar * radius_equatorial == 0.)  return 0.;

    // see equation (26) in A.Isihara, J.Chem.Phys. 18(1950)1446-1449
    const double ratio = (radius_polar < radius_equatorial
                          ? radius_polar / radius_equatorial
                          : radius_equatorial / radius_polar);
    const double e1 = sqrt(1.0 - ratio*ratio);
    const double b1 = 1.0 + asin(e1) / (e1 * ratio);
    const double bL = (1.0 + e1) / (1.0 - e1);
    const double b2 = 1.0 + 0.5 * ratio * ratio / e1 * log(bL);
    const double delta = 0.75 * b1 * b2;
    const double ddd = 2.0 * (delta + 1.0) * radius_polar * radius_equatorial * radius_equatorial;
    return 0.5 * cbrt(ddd);
}

static double
radius_effective(int mode, double radius_polar, double radius_equatorial)
{
    switch (mode) {
    default:
    case 1: // average curvature
        return radius_from_curvature(radius_polar, radius_equatorial);
    case 2: // equivalent volume sphere
        return radius_from_volume(radius_polar, radius_equatorial);
    case 3: // min radius
        return (radius_polar < radius_equatorial ? radius_polar : radius_equatorial);
    case 4: // max radius
        return (radius_polar > radius_equatorial ? radius_polar : radius_equatorial);
    }
}

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

static void rotated_long_axis(
    double gamma_1, double gamma_2,
    double psi, 
    double *result_rotated_axis) {
    double cos_psi, sin_psi;
    double cos_gamma_1, sin_gamma_1;
    double cos_gamma_2, sin_gamma_2;

    SINCOS(psi, sin_psi, cos_psi);
    SINCOS(gamma_1, sin_gamma_1, cos_gamma_1);
    SINCOS(gamma_2, sin_gamma_2, cos_gamma_2);

    double vector_long[3];
    double easy_cross_long[3];
    double vector_long_rotated[3];
    // B field is along (1,0,0), orientation of easy axis
    //                           is precessing in a cone around B
    // vector_long: long axis direction,
    //               vector lies in plane perpendicular to easy axis
    SET_VEC(vector_long, -sin_psi, cos_psi*sin_gamma_1, cos_psi*cos_gamma_1);

    // to obtain all possible vectors of vector_long, integrate gamma_2
    // use for this purpose Rodrigues rotation formula, need cross product:
    SET_VEC(easy_cross_long, 0.0, cos_gamma_1, -sin_gamma_1);
    // rotation formula:
    // vector_long_rotated = vector_long*cos_gamma_2 + easy_cross_long*sin_gamma_2
    SCALE_VEC(vector_long, cos_gamma_2);
    SCALE_VEC(easy_cross_long, sin_gamma_2);
    ADD_VEC(
        result_rotated_axis,
        vector_long,
        easy_cross_long
    );
    // last part rotation formula not necessary as r_easy, r_long
    // are perpendicular by definition
    // see: r_easy*dot_product(r_easy, r_long)*(1d0-cos(gamma_2))

    // rotate coordinate frame in case the magnetic field is not
    // pointing along (1,0,0)
//    SET_VEC(result_rotated_axis,
//        cos_beta*vector_long_rotated[0] - sin_beta*vector_long_rotated[2],
//        vector_long_rotated[1],
//        sin_beta*vector_long_rotated[0] + cos_beta*vector_long_rotated[2]
//    );
}

static double oriented_spindle_amplitude(
    double qx, double qy, double radius_polar, double radius_equatorial,
    double gamma_1, double gamma_2, double psi) {

    // get orientation of long axis for given gamma_1, psi, gamma_2 angles:
    double vec_long_rotated[3];
    rotated_long_axis(gamma_1, gamma_2, psi, vec_long_rotated);

    // what is angle between q vector and long axis?
    // acos results in NaN if argument is numerically slightly larger than 1
    const double q = sqrt(square(qx) + square(qy));
    double cos_alpha = 1.0;
    if (q > 0) {
        cos_alpha = (qx*vec_long_rotated[0] + qy*vec_long_rotated[1])/q;
    }

    // calculate the effective radius
    const double req2 = square(radius_equatorial);
    const double r_eff = sqrt( req2 + (square(radius_polar) - req2)*square(cos_alpha) );

    // effective radius for ellipsoidal model:
    // http://gisaxs.com/index.php/Form_Factor:Ellipsoid_of_revolution
    return sas_3j1x_x(q*r_eff);
}

static double boltzmann_statistics(
    double xi, double psi) {
    // xi * exp(xi * (cos(psi) - 1) ) / ( 1 - exp(-2xi) )

    if (xi < 0.0001) {
        // avoid dividing by zero
        // use first order taylor approximation of exp(x) in this case
        // which is accurate to 2e-8
        return 0.5*(1.0 + xi * (cos(psi) - 1.0));
    } else {
        return xi * exp(xi*(cos(psi) - 1.0)) / (1.0 - exp(-2*xi));
    }

}

static void
Fq(double q,
    double *F1,
    double *F2,
    double sld,
    double sld_solvent,
    double radius_polar,
    double radius_equatorial,
    double xi)
{
    // perform integration over gamma_1=0..2pi and gamma_2=0..2pi
    // to account for the precession motion of the easy axis around B with
    // cone angle psi
    // as well as the motion of the long axis of the spindle 
    // translate a point in [-1,1] to a point in [0, 1]
    // const double u = GAUSS_Z[i]*(upper-lower)/2 + (upper+lower)/2;
    
    
    double total_F1 = 0.0;
    double total_F2 = 0.0;
    for (int i_detangle=0; i_detangle < GAUSS_N; i_detangle++) {
        double detangle = M_PI_2 * (GAUSS_Z[i_detangle] + 1.0); // 0 .. pi
        double const qx = q*cos(detangle);
        double const qy = q*sin(detangle);

        double integral_psi_F1 = 0.0;
        double integral_psi_F2 = 0.0;
        for (int i_psi=0; i_psi < GAUSS_N; i_psi++) {
            const double psi = M_PI_2 * (GAUSS_Z[i_psi] + 1.0); // 0 .. pi
            double integral_g1_F1 = 0.0;
            double integral_g1_F2 = 0.0;
            for (int i_g1=0; i_g1 < GAUSS_N; i_g1++) {
                const double gamma_1 = M_PI * (GAUSS_Z[i_g1] + 1.0); // 0 .. 2 pi
                double integral_g2_F1 = 0.0;
                double integral_g2_F2 = 0.0;
                for (int i_g2=0; i_g2 < GAUSS_N; i_g2++) {
                    const double gamma_2 = M_PI * (GAUSS_Z[i_g2] + 1.0); // 0 .. 2 pi
                    const double f = oriented_spindle_amplitude(
                        qx, qy, radius_polar, radius_equatorial,
                        gamma_1, gamma_2, psi);
                    integral_g2_F1 += GAUSS_W[i_g2] * f;
                    integral_g2_F2 += GAUSS_W[i_g2] * f * f;
                }
                integral_g1_F1 += GAUSS_W[i_g1] * integral_g2_F1;
                integral_g1_F2 += GAUSS_W[i_g1] * integral_g2_F2;
            }
            const double sin_psi = sin(psi);
            const double boltzmann_weight = boltzmann_statistics(xi, psi);
            integral_psi_F1 += GAUSS_W[i_psi] * sin_psi * boltzmann_weight * integral_g1_F1;
            integral_psi_F2 += GAUSS_W[i_psi] * sin_psi * boltzmann_weight * integral_g1_F2;
        }
        total_F1 += GAUSS_W[i_detangle] * integral_psi_F1;
        total_F2 += GAUSS_W[i_detangle] * integral_psi_F2;
    }
    // factors for the translation of the [-1,1] range to
    // 0 .. 2pi, 0 .. 2pi, 0 .. pi -> pi^3 / 2
    // divide by 4 pi^2 for the integration over the whole space
    // divide by pi for integration over detector
    // -> factor of 1/8
    total_F1 *= 0.125;
    total_F2 *= 0.125;
    const double s = (sld - sld_solvent) * form_volume(radius_polar, radius_equatorial);
    *F1 = 1e-2 * s * total_F1;
    *F2 = 1e-4 * s * s * total_F2;
}

static double
Iqxy(double qx, double qy,
    double sld,
    double sld_solvent,
    double radius_polar,
    double radius_equatorial,
    double xi)
{
    // mu_bohr / k_Boltzmann = 0.6717140430498562
//    const double xi = 0.6717140430498562 * magnetic_moment * magnetic_field / temperature;
    double integral_psi_F2 = 0.0;
    for (int i_psi=0; i_psi < GAUSS_N; i_psi++) {
        const double psi = M_PI_2 * (GAUSS_Z[i_psi] + 1.0); // 0 .. pi
        double integral_g1_F2 = 0.0;
        for (int i_g1=0; i_g1 < GAUSS_N; i_g1++) {
            const double gamma_1 = M_PI * (GAUSS_Z[i_g1] + 1.0); // 0 .. 2 pi
            double integral_g2_F2 = 0.0;
            for (int i_g2=0; i_g2 < GAUSS_N; i_g2++) {
                const double gamma_2 = M_PI * (GAUSS_Z[i_g2] + 1.0); // 0 .. 2 pi
                const double f = oriented_spindle_amplitude(
                    qx, qy, radius_polar, radius_equatorial,
                    gamma_1, gamma_2, psi);
                integral_g2_F2 += GAUSS_W[i_g2] * f * f;
            }
            integral_g1_F2 += GAUSS_W[i_g1] * integral_g2_F2;
        }
        const double sin_psi = sin(psi);
        const double boltzmann_weight = boltzmann_statistics(xi, psi);
        integral_psi_F2 += GAUSS_W[i_psi] * sin_psi * boltzmann_weight * integral_g1_F2;
    }
    // factors for the translation of the [-1,1] range to
    // 0 .. 2pi, 0 .. 2pi, 0 .. pi -> pi^3 / 2
    // divide by 4 pi^2 for the integration over the whole space
    // -> factor of pi/8
    integral_psi_F2 *= 0.125 * M_PI;
    const double s = (sld - sld_solvent) * form_volume(radius_polar, radius_equatorial);

    return 1.0e-4 * square(s) * integral_psi_F2;
}

