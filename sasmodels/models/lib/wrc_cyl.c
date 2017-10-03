/*
    Functions for WRC implementation of flexible cylinders
*/

static double
Rgsquare(double L, double b)
{
    const double x = L/b;
    // Use Horner's method to evaluate:
    //     pow(1.0+square(x/3.12)+cube(x/8.67), 0.176/3.0)
    // Too many digits on the coefficients, but necessary for consistency
    const double alphasq =
        pow(1.0 + square(x)*(1.534414548417740e-03*x + 1.027284681130835e-01),
            5.866666666666667e-02);
    return alphasq*L*b/6.0;
}

static double
Rgsquareshort(double L, double b)
{
    const double r = b/L;
    return Rgsquare(L, b) * (1.0 + r*(-1.5 + r*(1.5 + r*0.75*expm1(-2.0/r))));
}

static double
a_long(double qp, double L, double b/*, double p1, double p2, double q0*/)
{
    // Note: caller sends p1, p2, q0 in as constants.
    const double p1 = 4.12;
    const double p2 = 4.42;
    const double q0 = 3.1;

    const double C1 = 1.22;
    const double C2 = 0.4288;
    const double C3 = -1.651;
    const double C4 = 1.523;
    const double C5 = 0.1477;
    const double miu = 0.585;


    const double C = (L/b>10.0 ? 3.06*pow(L/b, -0.44) : 1.0);
    const double r2 = Rgsquare(L,b);
    const double r = sqrt(r2);
    const double qr_b = q0*r/b;
    const double qr_b_sq = qr_b*qr_b;
    const double qr_b_4 = qr_b_sq*qr_b_sq;
    const double qr_b_miu = pow(qr_b, -1.0/miu);
    const double em1_qr_b_sq = expm1(-qr_b_sq);
    const double sech2 = 1.0/square(cosh((qr_b-C4)/C5));
    const double tanh1m = 1.0 - tanh(-C4 + qr_b/C5);

    const double t1 = pow(q0, 1.0 + p1 + p2)/(b*(p1-p2));
    const double t2 = C/(15.0*L) * (
        + 14.0*b*b*em1_qr_b_sq/(q0*qr_b_sq)
        + 2.0*q0*r2*exp(-qr_b_sq)*(11.0 + 7.0/qr_b_sq));
    const double t11 = ((C3*qr_b_miu + C2)*qr_b_miu + C1)*qr_b_miu;
    const double t3 = r*sech2/(2.*C5)*t11;
    const double t4 = r*(em1_qr_b_sq + qr_b_sq)*sech2 / (C5*qr_b_4);
    const double t5 = -2.0 * r*qr_b*em1_qr_b_sq * tanh1m / qr_b_4;
    const double t10 = (em1_qr_b_sq + qr_b_sq) * tanh1m / qr_b_4;
    const double t6 = 4.0*b/q0 * t10;
    const double t7 = r*((-3.0*C3*qr_b_miu -2.0*C2)*qr_b_miu -1.0*C1)*qr_b_miu/(miu*qr_b);
    const double t8 = 2.0 - tanh1m;
    const double t9 = b*C/(15.0*L) * (4.0 - exp(-qr_b_sq) * (11.0 + 7.0/qr_b_sq) + 7.0/qr_b_sq);
    const double t12 = b*b*M_PI/(L*q0*q0) + t2 + t3 - t4 + t5 - t6 + 0.5*t7*t8;
    const double t13 = -b*M_PI/(L*q0) + t9 + t10 + 0.5*t11*t8;

    const double a1 = pow(q0,p1)*t13 - t1*pow(q0,-p2)*(t12 + b*p1/q0*t13);
    const double a2 = t1*pow(q0,-p1)*(t12 + b*p1/q0*t13);

    const double ans = a1*pow(qp*b, -p1) + a2*pow(qp*b, -p2) + M_PI/(qp*L);
    return ans;
}

static double
_short(double r2, double exp_qr_b, double L, double b, double p1short,
        double p2short, double q0)
{
    const double qr2 = q0*q0 * r2;
    const double b3 = b*b*b;
    const double q0p = pow(q0, -4.0 + p1short);

    double yy = 1.0/(L*r2*r2) * b/exp_qr_b*q0p
        * (8.0*b3*L
           - 8.0*b3*exp_qr_b*L
           + 2.0*b3*exp_qr_b*L*p2short
           - 2.0*b*exp_qr_b*L*p2short*qr2
           + 4.0*b*exp_qr_b*L*qr2
           - 2.0*b3*L*p2short
           + 4.0*b*L*qr2
           - M_PI*exp_qr_b*qr2*q0*r2
           + M_PI*exp_qr_b*p2short*qr2*q0*r2);

    return yy;
}
static double
a_short(double qp, double L, double b
        /*double p1short, double p2short*/, double q0)
{
    const double p1short = 5.36;
    const double p2short = 5.62;

    const double r2 = Rgsquareshort(L,b);
    const double exp_qr_b = exp(r2*square(q0/b));
    const double pdiff = p1short - p2short;
    const double a1 = _short(r2,exp_qr_b,L,b,p1short,p2short,q0)/pdiff;
    const double a2= -_short(r2,exp_qr_b,L,b,p2short,p1short,q0)/pdiff;
    const double ans = a1*pow(qp*b, -p1short) + a2*pow(qp*b, -p2short) + M_PI/(qp*L);
    return ans;
}

//WR named this w (too generic)
static double
w_WR(double x)
{
    return 0.5*(1 + tanh((x - 1.523)/0.1477));
}

static double
Sdebye(double qsq)
{
#if FLOAT_SIZE>4
#define DEBYE_CUTOFF 0.1  // 1e-14 error
#else
#define DEBYE_CUTOFF 0.9  // 4e-7 error
#endif

/* For double precision, the following gets 1e-15 error rather than 1e-14
    if (qsq < 9./16.) {
        // PadeApproximant[2*Exp[-x^2] + x^2-1)/x^4, {x, 0, 8}]
        const double A1=1./12., A2=2./99., A3=1./2640., A4=1./23760., A5=-1./1995840.;
        const double B1=5./12., B2=5./66., B3=1./132., B4=1./2376., B5=1./95040.;
        const double x = qsq;
        return (((((A5*x + A4)*x + A3)*x + A2)*x + A1)*x + 1.)
                /(((((B5*x + B4)*x + B3)*x + B2)*x + B1)*x + 1.);
    }
*/

    if (qsq < DEBYE_CUTOFF) {
        const double x = qsq;
        const double C0 = +1.;
        const double C1 = -1./3.;
        const double C2 = +1./12.;
        const double C3 = -1./60.;
        const double C4 = +1./360.;
        const double C5 = -1./2520.;
        const double C6 = +1./20160.;
        const double C7 = -1./181440.;
        return ((((((C7*x + C6)*x + C5)*x + C4)*x + C3)*x + C2)*x + C1)*x + C0;
    } else {
        return 2.*(exp(-qsq) + qsq - 1.)/(qsq*qsq);
    }
}

static double
Sexv(double q, double L, double b)
{
    const double C1=1.22;
    const double C2=0.4288;
    const double C3=-1.651;
    const double miu = 0.585;
    const double qr = q*sqrt(Rgsquare(L,b));
    const double x = pow(qr, -1.0/miu);
    const double w = w_WR(qr);

    const double base = (1.0 - w)*Sdebye(qr*qr);
    const double correction = w*x*(C1 + x*(C2 + x*C3));
    return base + correction;
}


static double
Sexv_new(double q, double L, double b)
{
    // Modified by Yun on Oct. 15,

    // Correction factor to apply to the returned Sexv
    const double qr = q*sqrt(Rgsquare(L,b));
    const double qr2 = qr*qr;
    const double t1 = (4.0/15.0 + 7.0/(15.0*qr2) - (11.0/15.0 + 7.0/(15.0*qr2))*exp(-qr2))*(b/L);
    const double C = (L/b > 10.0) ? 3.06*pow(L/b, -0.44)*t1 : t1;

    const double Sexv_orig = Sexv(q, L, b);

    // calculating the derivative to decide on the correction (cutoff) term?
    // Note: this is modified from WRs original code
    const double del=1.05;
    const double qdel = (Sexv(q*del,L,b) - Sexv_orig)/(q*(del - 1.0));

    if (qdel < 0) {
        // return Sexv with the additional correction
        return Sexv_orig + C;
    } else {
        // recalculate Sexv base and return it with the additional correction
        const double w = w_WR(qr);
        return (1.0 - w)*Sdebye(qr2) + C;
    }
}


static double
Sk_WR(double q, double L, double b)
{
    const double Rg_short = sqrt(Rgsquareshort(L, b));
    double q0short = fmax(1.9/Rg_short, 3.0);
    double ans;

    if( L > 4*b ) { // L > 4*b : Longer Chains
        if (q*b <= 3.1) {
            ans = Sexv_new(q, L, b);
        } else { //q(i)*b > 3.1
            ans = a_long(q, L, b /*, p1, p2, q0*/);
        }
    } else { // L <= 4*b : Shorter Chains
        if (q*b <= q0short) { // q*b <= fmax(1.9/Rg_short, 3)
            if (q*b <= 0.01) {
                ans = 1.0 - square(q*Rg_short)/3.0;
            } else {
                ans = Sdebye(square(q*Rg_short));
            }
        } else {  // q*b > max(1.9/Rg_short, 3)
            ans = a_short(q, L, b /*, p1short, p2short*/, q0short);
        }
    }

    return ans;
}
