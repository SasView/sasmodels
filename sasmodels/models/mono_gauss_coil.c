static double form_volume(double rg)
{
    return 1.0;
}

static double
effective_radius(int mode, double rg)
{
    if (mode == 1) {
        return rg;
    } else if (mode == 2) {
        return 2.0*rg;
    } else if (mode == 3) {
        return 3.0*rg;
    } else {
        return sqrt(5.0/3.0)*rg;
    }
}

double Iq(double q, double i_zero, double rg)
{
    const double uarg = square(q*rg);
    const double inten;
    if (q == 0) {
        inten = i_zero;
    } else {
        inten = 2.0*i_zero * (exp(-uarg) + uarg - 1.0)/square(uarg);
    }

    return inten;
}
