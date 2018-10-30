// Set OVERLAPPING to 1 in order to fill in the edges of the box, with
// c endcaps and b overlapping a.  With the proper choice of parameters,
// (setting rim slds to sld, core sld to solvent, rim thickness to thickness
// and subtracting 2*thickness from length, this should match the hollow
// rectangular prism.)  Set it to 0 for the documented behaviour.
#define OVERLAPPING 0
static double
form_volume(double length_a, double length_b, double length_c,
    double thick_rim_a, double thick_rim_b, double thick_rim_c)
{
    return
#if OVERLAPPING
        // Hollow rectangular prism only includes the volume of the shell
        // so uncomment the next line when comparing.  Solid rectangular
        // prism, or parallelepiped want filled cores, so comment when
        // comparing.
        //-length_a * length_b * length_c +
        (length_a + 2.0*thick_rim_a) *
        (length_b + 2.0*thick_rim_b) *
        (length_c + 2.0*thick_rim_c);
#else
        length_a * length_b * length_c +
        2.0 * thick_rim_a * length_b * length_c +
        2.0 * length_a * thick_rim_b * length_c +
        2.0 * length_a * length_b * thick_rim_c;
#endif
}

static double
radius_from_volume(double length_a, double length_b, double length_c,
                   double thick_rim_a, double thick_rim_b, double thick_rim_c)
{
    const double volume = form_volume(length_a, length_b, length_c, thick_rim_a, thick_rim_b, thick_rim_c);
    return cbrt(volume/M_4PI_3);
}

static double
radius_from_crosssection(double length_a, double length_b, double thick_rim_a, double thick_rim_b)
{
    const double area_xsec_paral = length_a*length_b + 2.0*thick_rim_a*length_b + 2.0*thick_rim_b*length_a;
    return sqrt(area_xsec_paral/M_PI);
}

static double
effective_radius(int mode, double length_a, double length_b, double length_c,
                 double thick_rim_a, double thick_rim_b, double thick_rim_c)
{
    switch (mode) {
    default:
    case 1: // equivalent sphere
        return radius_from_volume(length_a, length_b, length_c, thick_rim_a, thick_rim_b, thick_rim_c);
    case 2: // half outer length a
        return 0.5 * length_a + thick_rim_a;
    case 3: // half outer length b
        return 0.5 * length_b + thick_rim_b;
    case 4: // half outer length c
        return 0.5 * length_c + thick_rim_c;
    case 5: // equivalent circular cross-section
        return radius_from_crosssection(length_a, length_b, thick_rim_a, thick_rim_b);
    case 6: // half outer ab diagonal
        return 0.5*sqrt(square(length_a+ 2.0*thick_rim_a) + square(length_b+ 2.0*thick_rim_b));
    case 7: // half outer diagonal
        return 0.5*sqrt(square(length_a+ 2.0*thick_rim_a) + square(length_b+ 2.0*thick_rim_b) + square(length_c+ 2.0*thick_rim_c));
    }
}

static void
Fq(double q,
    double *F1,
    double *F2,
    double core_sld,
    double arim_sld,
    double brim_sld,
    double crim_sld,
    double solvent_sld,
    double length_a,
    double length_b,
    double length_c,
    double thick_rim_a,
    double thick_rim_b,
    double thick_rim_c)
{
    // Code converted from functions CSPPKernel and CSParallelepiped in libCylinder.c
    // Did not understand the code completely, it should be rechecked (Miguel Gonzalez)
    // Code is rewritten, the code is compliant with Diva Singh's thesis now (Dirk Honecker)
    // Code rewritten; cross checked against hollow rectangular prism and realspace (PAK)

    const double half_q = 0.5*q;

    const double tA = length_a + 2.0*thick_rim_a;
    const double tB = length_b + 2.0*thick_rim_b;
    const double tC = length_c + 2.0*thick_rim_c;

    // Scale factors
    const double dr0 = (core_sld-solvent_sld);
    const double drA = (arim_sld-solvent_sld);
    const double drB = (brim_sld-solvent_sld);
    const double drC = (crim_sld-solvent_sld);

    // outer integral (with gauss points), integration limits = 0, 1
    // substitute d_cos_alpha for sin_alpha d_alpha
    double outer_sum_F1 = 0; //initialize integral
    double outer_sum_F2 = 0; //initialize integral
    for( int i=0; i<GAUSS_N; i++) {
        const double cos_alpha = 0.5 * ( GAUSS_Z[i] + 1.0 );
        const double mu = half_q * sqrt(1.0-cos_alpha*cos_alpha);
        const double siC = length_c * sas_sinx_x(length_c * cos_alpha * half_q);
        const double siCt = tC * sas_sinx_x(tC * cos_alpha * half_q);

        // inner integral (with gauss points), integration limits = 0, 1
        // substitute beta = PI/2 u (so 2/PI * d_(PI/2 * beta) = d_beta)
        double inner_sum_F1 = 0.0;
        double inner_sum_F2 = 0.0;
        for(int j=0; j<GAUSS_N; j++) {
            const double u = 0.5 * ( GAUSS_Z[j] + 1.0 );
            double sin_beta, cos_beta;
            SINCOS(M_PI_2*u, sin_beta, cos_beta);
            const double siA = length_a * sas_sinx_x(length_a * mu * sin_beta);
            const double siB = length_b * sas_sinx_x(length_b * mu * cos_beta);
            const double siAt = tA * sas_sinx_x(tA * mu * sin_beta);
            const double siBt = tB * sas_sinx_x(tB * mu * cos_beta);

#if OVERLAPPING
            const double f = dr0*siA*siB*siC
                + drA*(siAt-siA)*siB*siC
                + drB*siAt*(siBt-siB)*siC
                + drC*siAt*siBt*(siCt-siC);
#else
            const double f = dr0*siA*siB*siC
                + drA*(siAt-siA)*siB*siC
                + drB*siA*(siBt-siB)*siC
                + drC*siA*siB*(siCt-siC);
#endif

            inner_sum_F1 += GAUSS_W[j] * f;
            inner_sum_F2 += GAUSS_W[j] * f * f;
        }
        // now complete change of inner integration variable (1-0)/(1-(-1))= 0.5
        // and sum up the outer integral
        outer_sum_F1 += GAUSS_W[i] * inner_sum_F1 * 0.5;
        outer_sum_F2 += GAUSS_W[i] * inner_sum_F2 * 0.5;
    }
    // now complete change of outer integration variable (1-0)/(1-(-1))= 0.5
    outer_sum_F1 *= 0.5;
    outer_sum_F2 *= 0.5;

    //convert from [1e-12 A-1] to [cm-1]
    *F1 = 1.0e-2 * outer_sum_F1;
    *F2 = 1.0e-4 * outer_sum_F2;
}

static double
Iqabc(double qa, double qb, double qc,
    double core_sld,
    double arim_sld,
    double brim_sld,
    double crim_sld,
    double solvent_sld,
    double length_a,
    double length_b,
    double length_c,
    double thick_rim_a,
    double thick_rim_b,
    double thick_rim_c)
{
    // cspkernel in csparallelepiped recoded here
    const double dr0 = core_sld-solvent_sld;
    const double drA = arim_sld-solvent_sld;
    const double drB = brim_sld-solvent_sld;
    const double drC = crim_sld-solvent_sld;

    const double tA = length_a + 2.0*thick_rim_a;
    const double tB = length_b + 2.0*thick_rim_b;
    const double tC = length_c + 2.0*thick_rim_c;
    const double siA = length_a*sas_sinx_x(0.5*length_a*qa);
    const double siB = length_b*sas_sinx_x(0.5*length_b*qb);
    const double siC = length_c*sas_sinx_x(0.5*length_c*qc);
    const double siAt = tA*sas_sinx_x(0.5*tA*qa);
    const double siBt = tB*sas_sinx_x(0.5*tB*qb);
    const double siCt = tC*sas_sinx_x(0.5*tC*qc);

#if OVERLAPPING
    const double f = dr0*siA*siB*siC
        + drA*(siAt-siA)*siB*siC
        + drB*siAt*(siBt-siB)*siC
        + drC*siAt*siBt*(siCt-siC);
#else
    const double f = dr0*siA*siB*siC
        + drA*(siAt-siA)*siB*siC
        + drB*siA*(siBt-siB)*siC
        + drC*siA*siB*(siCt-siC);
#endif

    return 1.0e-4 * f * f;
}
