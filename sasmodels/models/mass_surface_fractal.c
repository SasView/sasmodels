double form_volume(double radius);

double Iq(double q,
          double mass_dim,
          double surface_dim,
          double cluster_rg,
          double primary_rg);

static double _mass_surface_fractal_kernel(double q,
          double mass_dim,
          double surface_dim,
          double cluster_rg,
          double primary_rg)
{
     //computation
    double tot_dim = 6.0 - surface_dim - mass_dim;
    mass_dim /= 2.0;
    tot_dim /= 2.0;

    double rc_norm = cluster_rg * cluster_rg / (3.0 * mass_dim);
    double rp_norm = primary_rg * primary_rg / (3.0 * tot_dim);

    //x for P
    double x_val1 = 1.0 +  q * q * rc_norm;
    double x_val2 = 1.0 +  q * q * rp_norm;

    double inv_form = pow(x_val1, mass_dim) * pow(x_val2, tot_dim);

    //another singular
    if (inv_form == 0.0) return 0.0;

    double form_factor = 1.0;
    form_factor /= inv_form;

    return (form_factor);
}
double form_volume(double radius){

    return 1.333333333333333*M_PI*radius*radius*radius;
}

double Iq(double q,
          double mass_dim,
          double surface_dim,
          double cluster_rg,
          double primary_rg)
{
    return _mass_surface_fractal_kernel(q,
            mass_dim,
            surface_dim,
            cluster_rg,
            primary_rg);
}
