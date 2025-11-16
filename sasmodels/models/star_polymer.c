double form_volume(void);

double Iq(double q, double rg_sq, double arms);

static double star_polymer_kernel(double q, double rg_sq, int n_arms)
{
    // Given the kernel function:
    //    I = 2(v + (e^-v - 1) + (a-1)/2 (e^-v - 1)^2) / av^2
    // At q = 0, u_2 and v are 0 so term1 = term2 = 0, hence 0/0 = nan
    // We can clean up the values at low q using Taylor expansions for
    //    T1 = (e^-v - 1)/v = [-1 +1/2 -1/6 +1/24 -1/120 +1/720 ...]
    //    T2 = [(e^-v - 1)/v]^2 = [+1 -1 +7/12 -1/4 +31/360 -1/40 ...]
    // Rewriting
    //    I = 2/av + 2/av (e^-v - 1)/v + (a-1)/a ((e^v - 1)/v)^2
    // Since the first term of T1 is -1 this cancels the leading 2/av.
    // Muliplying the remaining terms by 2/v yields T1' = 2*T1[1:]
    //    I = 1/a T1'(v) + (1 - 1/a) T2(v) = (T1'(v) - T2(v))/a + T2(v)
    //    T1' = [+1 -1/3 +1/12 -1/60 +1/360 -1/2520]
    // Substituting v=0 we get lim[q->0] I = (1 - 1)/a + 1 = 1 as expected.
    // After testing it turns out that term2 is numerically stable, so we
    // only use the Taylor series for T1', and protect T2 against v = 0.

#if FLOAT_SIZE>4
#  define _SPLIM 0.03  // cutoff for (Q Rg)^2
#else
#  define _SPLIM 1.0
#endif

    if (q == 0.) { // IEEE 754 has -0 == +0
        return 1.0;
    }

    // Note: arms should be one or more
    const double u_sq = rg_sq * q * q;
    const double v = (u_sq * n_arms) / (3 * n_arms - 2);
    double term1;
    if (u_sq < _SPLIM) {
        term1 = 1. + v*(-1./3. + v*(1./12. + v*(-1./60. + v*(1./360. + v*(-1./2520)))));
    } else {
        term1 = 2.0 * (v + expm1(-v)) / (v * v);
    }
    double term2 = square(expm1(-v)/v); // we trap q == 0 above, so v > 0 here
    return (term1 - term2)/n_arms + term2;
}

double form_volume(void)
{
    return 1.0;
}

double Iq(double q, double rg_sq, double arms)
{
    // TODO: round or truncate?
    // rounding: stacked_disk pearl_necklace multilayer_vesicle linear_pearls
    //      lamellar_stack_caille lamellar_hg_stack_caille
    // truncation: lamellar_stack_paracrystal
    return star_polymer_kernel(q, rg_sq, (int)(arms + 0.5));
}
