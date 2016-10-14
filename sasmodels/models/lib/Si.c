// integral of sin(x)/x Taylor series approximated to w/i 0.1%
double Si(double x);
double Si(double x)
{
    if (x >= M_PI*6.2/4.0){
        const double z = 1./(x*x);
        // Explicitly writing factorial values triples the speed of the calculation
        const double out_cos = (((-720.*z + 24.)*z - 2.)*z + 1.)/x;
        const double out_sin = (((-5040.*z + 120.)*z - 6.)*z + 1)*z;

        double cos_x, sin_x;
        SINCOS(x, cos_x, sin_x);
        return M_PI_2 - cos_x*out_cos - sin_x*out_sin;
    } else {
        const double z = x*x;
        // Explicitly writing factorial values triples the speed of the calculation
        return (((((-1./439084800.*z
            + 1./3265920.)*z
            - 1./35280.)*z
            + 1./600.)*z
            - 1./18.)*z
            + 1.)*x;
    }
}