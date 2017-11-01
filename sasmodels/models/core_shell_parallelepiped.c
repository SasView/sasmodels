static double
form_volume(double length_a, double length_b, double length_c,
    double thick_rim_a, double thick_rim_b, double thick_rim_c)
{
    return length_a * length_b * length_c +
           2.0 * thick_rim_a * length_b * length_c +
           2.0 * thick_rim_b * length_a * length_c +
           2.0 * thick_rim_c * length_a * length_b;
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
    double VtA = (2.0 * thick_rim_a * length_b * length_c);
    double VtB = (2.0 * length_a * thick_rim_b * length_c);
    double VtC = (2.0 * length_a * length_b * thick_rim_c);

    // Scale factors (note that drC is not used later)
    const double dr0 = (core_sld-solvent_sld);
    const double drA = (arim_sld-solvent_sld);
    const double drB = (brim_sld-solvent_sld);
    const double drC = (crim_sld-solvent_sld);

    // Precompute scale factors for combining cross terms from the shape
    const double dr0_Vin = dr0*Vin;
    const double drA_VtA = drA*VtA;
    const double drB_VtB = drB*VtB;
    const double drC_VtC = drC*VtC;
    const double drV_delta = dr0_Vin - drA_VtA - drB_VtB - drC_VtC;

    /*  *************** algorithm description ******************

    // Rewrite f as x*siC + y*siCt to move the siC/siCt calculation out
    // of the inner loop.  That is:

    f = (di-da-db-dc) sa sb sc + da sa' sb sc + db sa sb' sc + dc sa sb sc'
      =  [ (di-da-db-dc) sa sb + da sa' sb + db sa sb' ] sc  + [dc sa sb] sc'
      = x sc + y sc'

    // where:
    di = delta rho_core V_core
    da = delta rho_rimA V_rimA
    db = delta rho_rimB V_rimB
    dc = delta rho_rimC V_rimC
    sa = j_0 (q_a a/2)    // siA the code
    sb = j_0 (q_b b/2)
    sc = j_0 (q_c c/2)
    sa' = j_0(q_a a_rim/2)  // siAt the code
    sb' = j_0(q_b b_rim/2)
    sc' = j_0(q_c c_rim/2)

    // qa, qb, and qc are generated using polar coordinates, with the
    // outer loop integrating over [0,1] after the u-substitution
    //    sigma = cos(theta), sqrt(1-sigma^2) = sin(theta)
    // and inner loop integrating over [0, pi/2] as
    //    uu = phi

    ************************************************************  */

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

            const double x = drV_delta*siA*siB + drA_VtA*siB*siAt + drB_VtB*siA*siBt;
            const double form = x*siC + drC_VtC*siA*siB*siCt;

            inner_sum += Gauss76Wt[j] * form * form;
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
    double VtA = 2.0 * thick_rim_a * length_b * length_c;
    double VtB = 2.0 * length_a * thick_rim_b * length_c;
    double VtC = 2.0 * length_a * length_b * thick_rim_c;

    // The definitions of ta, tb, tc are not the same as in the 1D case because there is no
    // the scaling by B.
    double tA = length_a + 2.0*thick_rim_a;
    double tB = length_b + 2.0*thick_rim_b;
    double tC = length_c + 2.0*thick_rim_c;
    //handle arg=0 separately, as sin(t)/t -> 1 as t->0
    double siA = sas_sinx_x(0.5*length_a*qa);
    double siB = sas_sinx_x(0.5*length_b*qb);
    double siC = sas_sinx_x(0.5*length_c*qc);
    double siAt = sas_sinx_x(0.5*tA*qa);
    double siBt = sas_sinx_x(0.5*tB*qb);
    double siCt = sas_sinx_x(0.5*tC*qc);


    // f uses Vin, V1, V2, and V3 and it seems to have more sense than the value computed
    // in the 1D code, but should be checked!
    double f = ( dr0*Vin*siA*siB*siC
               + drA*VtA*(siAt-siA)*siB*siC
               + drB*VtB*siA*(siBt-siB)*siC
               + drC*VtC*siA*siB*(siCt-siC));

    return 1.0e-4 * f * f;
}
