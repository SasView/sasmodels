#define INVALID(p) (p.fractal_dim <= 0.0)

static double
Iq(double q,
   double volfraction,
   double radius,
   double fractal_dim,
   double cor_length,
   double sld_block,
   double sld_solvent)
{
    //calculate P(q) for the spherical subunits
    const double pq = M_4PI_3*cube(radius) * square((sld_block-sld_solvent)*sph_j1c(q*radius));
    
    //calculate S(q),  using Teixeira, Eq(15)
    double sq;
    if (q > 0. && fractal_dim > 1.) {
        // q>0, D>0
        const double D = fractal_dim;
        const double Dm1 = fractal_dim - 1.0;
        // Note: for large Dm1, sin(Dm1*atan(q*cor_length) can go negative
        const double t1 = D*sas_gamma(Dm1)*sin(Dm1*atan(q*cor_length));
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

    // scale to units cm-1 sr-1 (assuming data on absolute scale)
    //    convert I(1/A) to (1/cm)  => 1e8 * I(q)
    //    convert rho^2 in 10^-6 1/A to 1/A  => 1e-12 * I(q)
    //    combined: 1e-4 * I(q)

    return 1.e-4 * volfraction * pq * sq;
}

