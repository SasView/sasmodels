static double
form_volume(double radius, double thickness)
{
    return M_4PI_3 * cube(radius + thickness);
}

static double
Iq(double q,
   double radius,
   double thickness,
   double core_sld,
   double shell_sld,
   double solvent_sld,
   double volfraction,
   double fractal_dim,
   double cor_length)
{
    //The radius for the building block of the core shell particle that is 
    //needed by the Teixeira fractal S(q) is the radius of the whole particle.
    const double cs_radius = radius + thickness;
    const double sq = fractal_sq(q, cs_radius, fractal_dim, cor_length);
    const double pq = core_shell_kernel(q, radius, thickness,
                                        core_sld, shell_sld, solvent_sld);

    // Note: core_shell_kernel already performs the 1e-4 unit conversion
    return volfraction * sq * pq;
}

