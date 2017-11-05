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
Iq(double q,
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
    //Code is rewritten,the code is compliant with Diva Singhs thesis now (Dirk Honecker)

    const double mu = 0.5 * q * length_b;

    // Scale sides by B
    const double a_over_b = length_a / length_b;
    const double c_over_b = length_c / length_b;

    double tA_over_b = a_over_b + 2.0*thick_rim_a/length_b;
    double tB_over_b = 1+ 2.0*thick_rim_b/length_b;
    double tC_over_b = c_over_b + 2.0*thick_rim_c/length_b;

    double Vin = length_a * length_b * length_c;
#if OVERLAPPING
    const double capA_area = length_b*length_c;
    const double capB_area = (length_a+2.*thick_rim_a)*length_c;
    const double capC_area = (length_a+2.*thick_rim_a)*(length_b+2.*thick_rim_b);
#else
    const double capA_area = length_b*length_c;
    const double capB_area = length_a*length_c;
    const double capC_area = length_a*length_b;
#endif
    const double Va = length_a * capA_area;
    const double Vb = length_b * capB_area;
    const double Vc = length_c * capC_area;
    const double Vat = Va + 2.0 * thick_rim_a * capA_area;
    const double Vbt = Vb + 2.0 * thick_rim_b * capB_area;
    const double Vct = Vc + 2.0 * thick_rim_c * capC_area;

    // Scale factors (note that drC is not used later)
    const double dr0 = (core_sld-solvent_sld);
    const double drA = (arim_sld-solvent_sld);
    const double drB = (brim_sld-solvent_sld);
    const double drC = (crim_sld-solvent_sld);

    // outer integral (with gauss points), integration limits = 0, 1
    double outer_sum = 0; //initialize integral
    for( int i=0; i<76; i++) {
        double sigma = 0.5 * ( Gauss76Z[i] + 1.0 );
        double mu_proj = mu * sqrt(1.0-sigma*sigma);

        // inner integral (with gauss points), integration limits = 0, pi/2
        const double siC = sas_sinx_x(mu * sigma * c_over_b);
        const double siCt = sas_sinx_x(mu * sigma * tC_over_b);
        double inner_sum = 0.0;
        for(int j=0; j<76; j++) {
            const double uu = 0.5 * ( Gauss76Z[j] + 1.0 );
            double sin_uu, cos_uu;
            SINCOS(M_PI_2*uu, sin_uu, cos_uu);
            const double siA = sas_sinx_x(mu_proj * sin_uu * a_over_b);
            const double siB = sas_sinx_x(mu_proj * cos_uu );
            const double siAt = sas_sinx_x(mu_proj * sin_uu * tA_over_b);
            const double siBt = sas_sinx_x(mu_proj * cos_uu * tB_over_b);

#if OVERLAPPING
            const double f = dr0*Vin*siA*siB*siC
                + drA*(Vat*siAt-Va*siA)*siB*siC
                + drB*siAt*(Vbt*siBt-Vb*siB)*siC
                + drC*siAt*siBt*(Vct*siCt-Vc*siC);
#else
            const double f = dr0*Vin*siA*siB*siC
                + drA*(Vat*siAt-Va*siA)*siB*siC
                + drB*siA*(Vbt*siBt-Vb*siB)*siC
                + drC*siA*siB*(Vct*siCt-Vc*siC);
#endif

            inner_sum += Gauss76Wt[j] * f * f;
        }
        inner_sum *= 0.5;
        // now sum up the outer integral
        outer_sum += Gauss76Wt[i] * inner_sum;
    }
    outer_sum *= 0.5;

    //convert from [1e-12 A-1] to [cm-1]
    return 1.0e-4 * outer_sum;
}

static double
Iqxy(double qa, double qb, double qc,
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

    double Vin = length_a * length_b * length_c;
#if OVERLAPPING
    const double capA_area = length_b*length_c;
    const double capB_area = (length_a+2.*thick_rim_a)*length_c;
    const double capC_area = (length_a+2.*thick_rim_a)*(length_b+2.*thick_rim_b);
#else
    const double capA_area = length_b*length_c;
    const double capB_area = length_a*length_c;
    const double capC_area = length_a*length_b;
#endif
    const double Va = length_a * capA_area;
    const double Vb = length_b * capB_area;
    const double Vc = length_c * capC_area;
    const double Vat = Va + 2.0 * thick_rim_a * capA_area;
    const double Vbt = Vb + 2.0 * thick_rim_b * capB_area;
    const double Vct = Vc + 2.0 * thick_rim_c * capC_area;

    // The definitions of ta, tb, tc are not the same as in the 1D case because there is no
    // the scaling by B.
    const double tA = length_a + 2.0*thick_rim_a;
    const double tB = length_b + 2.0*thick_rim_b;
    const double tC = length_c + 2.0*thick_rim_c;
    const double siA = sas_sinx_x(0.5*length_a*qa);
    const double siB = sas_sinx_x(0.5*length_b*qb);
    const double siC = sas_sinx_x(0.5*length_c*qc);
    const double siAt = sas_sinx_x(0.5*tA*qa);
    const double siBt = sas_sinx_x(0.5*tB*qb);
    const double siCt = sas_sinx_x(0.5*tC*qc);

#if OVERLAPPING
    const double f = dr0*Vin*siA*siB*siC
        + drA*(Vat*siAt-Va*siA)*siB*siC
        + drB*siAt*(Vbt*siBt-Vb*siB)*siC
        + drC*siAt*siBt*(Vct*siCt-Vc*siC);
#else
    const double f = dr0*Vin*siA*siB*siC
        + drA*(Vat*siAt-Va*siA)*siB*siC
        + drB*siA*(Vbt*siBt-Vb*siB)*siC
        + drC*siA*siB*(Vct*siCt-Vc*siC);
#endif

    return 1.0e-4 * f * f;
}
