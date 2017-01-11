#define INVALID(p) (p.fractal_dim < 0.0)

static double
Iq(double q,
   double volfraction,
   double radius,
   double fractal_dim,
   double cor_length,
   double sld_block,
   double sld_solvent)
{
    const double sq = fractal_sq(q, radius, fractal_dim, cor_length);

    //calculate P(q) for the spherical subunits
    const double V = M_4PI_3*cube(radius);
    const double pq = V * square((sld_block-sld_solvent)*sas_3j1x_x(q*radius));

    // scale to units cm-1 sr-1 (assuming data on absolute scale)
    //    convert I(1/A) to (1/cm)  => 1e8 * I(q)
    //    convert rho^2 in 10^-6 1/A to 1/A  => 1e-12 * I(q)
    //    combined: 1e-4 * I(q)

    return 1.e-4 * volfraction * sq * pq;
}

