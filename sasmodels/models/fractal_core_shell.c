double form_volume(double radius, double thickness);

double Iq(double q,
          double radius,
          double thickness,
          double core_sld,
          double shell_sld,
          double solvent_sld,
          double volfraction,
          double fractal_dim,
          double cor_length);

double form_volume(double radius, double thickness)
{
    return M_4PI_3 * cube(radius + thickness);
}

double Iq(double q,
          double radius,
          double thickness,
          double core_sld,
          double shell_sld,
          double solvent_sld,
          double volfraction,
          double fractal_dim,
          double cor_length) {


    const double pq = core_shell_kernel(q, radius, thickness,
                                        core_sld, shell_sld, solvent_sld);


    //calculate S(q)
    double sq;
    if (q > 0. && fractal_dim > 1.) {
        // q>0, D>0
        const double D = fractal_dim;
        const double Dm1 = fractal_dim - 1.0;
        const double t1 = D*sas_gamma(Dm1)*sin((Dm1)*atan(q*cor_length));
        const double t2 = pow(q*radius, -D);
        const double t3 = pow(1.0 + 1.0/square(q*cor_length), -0.5*Dm1);
        sq = 1.0 + t1 * t2 * t3;
    } else if (q > 0.) {
        // q>0, D=1
        sq = 1.0 + atan(q*cor_length) / (q*radius);
    } else if (fractal_dim > 1.) {
        // q=0, D>1
        const double D = fractal_dim;
        sq = 1.0 + pow(cor_length/radius, D)*sas_gamma(D+1.0);
    } else {
        // q=0, D=1
        sq = 1.0 + cor_length/radius;
    }

    // Note: core_shell_kernel already performs the 1e-4 unit conversion
    return volfraction * sq * pq;
}

