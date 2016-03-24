/*
    Functions for WRC implementation of flexible cylinders
*/
double Sk_WR(double q, double L, double b);


static double
AlphaSquare(double x)
{
    // Potentially faster. Needs proper benchmarking.
    // add native_powr to kernel_template
    //double t = native_powr( (1.0 + (x/3.12)*(x/3.12) +
    //     (x/8.67)*(x/8.67)*(x/8.67)),(0.176/3.0) );
    //return t;

    return pow( (1.0 + (x/3.12)*(x/3.12) +
         (x/8.67)*(x/8.67)*(x/8.67)),(0.176/3.0) );
}

//
static double
Rgsquarezero(double q, double L, double b)
{
    const double r = b/L;
    return (L*b/6.0) *
           (1.0 - r*1.5  + 1.5*r*r - 0.75*r*r*r*(1.0 - exp(-2.0/r)));
}

//
static double
Rgsquareshort(double q, double L, double b)
{
    return AlphaSquare(L/b) * Rgsquarezero(q,L,b);
}

//
static double
Rgsquare(double q, double L, double b)
{
    return AlphaSquare(L/b)*L*b/6.0;
}

static inline double
sech_WR(double x)
{
    return(1/cosh(x));
}

static double
a1long(double q, double L, double b, double p1, double p2, double q0)
{
    double C;
    const double onehalf = 1.0/2.0;

    if( L/b > 10.0) {
        C = 3.06/pow((L/b),0.44);
    } else {
        C = 1.0;
    }

    const double C1 = 1.22;
    const double C2 = 0.4288;
    const double C3 = -1.651;
    const double C4 = 1.523;
    const double C5 = 0.1477;
    const double miu = 0.585;

    const double Rg2 = Rgsquare(q,L,b);
    const double Rg22 = Rg2*Rg2;
    const double Rg = sqrt(Rg2);
    const double Rgb = Rg*q0/b;

    const double b2 = b*b;
    const double b3 = b*b*b;
    const double b4 = b3*b;
    const double q02 = q0*q0;
    const double q03 = q0*q0*q0;
    const double q04 = q03*q0;
    const double q05 = q04*q0;

    const double Rg02 = Rg2*q02;

    const double t1 = (b*C*((4.0/15.0 - pow((double)M_E,(-(Rg02/b2))) *
         ((11.0/15.0 + (7.0*b2)/(15.0*Rg02))) +
         (7.0*b2)/(15.0*Rg02))));

    const double t2 = (2.0*b4*(((-1.0) + pow((double)M_E,(-(Rg02/b2))) +
         Rg02/b2))*((1.0 + onehalf*(((-1.0) -
         tanh((-C4 + Rgb/C5)))))));

    const double t3 = ((C3*pow(Rgb,((-3.0)/miu)) +
         C2*pow(Rgb,((-2.0)/miu)) +
         C1*pow(Rgb,((-1.0)/miu))));

    const double t4 = ((1.0 + tanh(((-C4) + Rgb)/C5)));

    const double t5 = (1.0/(b*p1*pow(q0,((-1.0) - p1 - p2)) -
         b*p2*pow(q0,((-1.0) - p1 - p2))));

    const double t6 = (b*C*(((-((14.0*b3)/(15.0*q03*Rg2))) +
         (14.0*b3*pow((double)M_E,(-(Rg02/b2))))/(15.0*q03*Rg2) +
         (2.0*pow((double)M_E,(-(Rg02/b2)))*q0*((11.0/15.0 +
         (7.0*b2)/(15.0*Rg02)))*Rg2)/b)));

    const double t7 = (Rg*((C3*pow(((Rgb)),((-3.0)/miu)) +
         C2*pow(((Rgb)),((-2.0)/miu)) +
         C1*pow(((Rgb)),((-1.0)/miu))))*pow(sech_WR(((-C4) +
         Rgb)/C5),2));

    const double t8 = (b4*Rg*(((-1.0) + pow((double)M_E,(-(Rg02/b2))) +
         Rg02/b2))*pow(sech_WR(((-C4) + Rgb)/C5),2));

    const double t9 = (2.0*b4*(((2.0*q0*Rg2)/b -
         (2.0*pow((double)M_E,(-(Rg02/b2)))*q0*Rg2)/b))*((1.0 + onehalf*(((-1.0) -
         tanh(((-C4) + Rgb)/C5))))));

    const double t10 = (8.0*b4*b*(((-1.0) + pow((double)M_E,(-(Rg02/b2))) +
         Rg02/b2))*((1.0 + onehalf*(((-1.0) - tanh(((-C4) +
         Rgb)/C5))))));

    const double t11 = (((-((3.0*C3*Rg*pow(((Rgb)),((-1.0) -
          3.0/miu)))/miu)) - (2.0*C2*Rg*pow(((Rgb)),((-1.0) -
          2.0/miu)))/miu - (C1*Rg*pow(((Rgb)),((-1.0) -
          1.0/miu)))/miu));

    const double t12 = ((1.0 + tanh(((-C4) + Rgb)/C5)));

    const double t13 = (b*C*((4.0/15.0 - pow((double)M_E,(-(Rg02/b2)))*((11.0/15.0 +
          (7.0*b2)/(15.0*q02* Rg2))) +
          (7.0*b2)/(15.0*Rg02))));

    const double t14 = (2.0*b4*(((-1.0) + pow((double)M_E,(-(Rg02/b2))) +
          Rg02/b2))*((1.0 + onehalf*(((-1.0) - tanh(((-C4) +
          Rgb)/C5))))));

    const double t15 = ((C3*pow(((Rgb)),((-3.0)/miu)) +
        C2*pow(((Rgb)),((-2.0)/miu)) +
        C1*pow(((Rgb)),((-1.0)/miu))));


    double yy = (pow(q0,p1)*(((-((b*M_PI)/(L*q0))) +t1/L +t2/(q04*Rg22) +
        onehalf*t3*t4)) + (t5*((pow(q0,(p1 - p2))*
        (((-pow(q0,(-p1)))*(((b2*M_PI)/(L*q02) +t6/L +t7/(2.0*C5) -
        t8/(C5*q04*Rg22) + t9/(q04*Rg22) -t10/(q05*Rg22) + onehalf*t11*t12)) -
        b*p1*pow(q0,((-1.0) - p1))*(((-((b*M_PI)/(L*q0))) + t13/L +
        t14/(q04*Rg22) + onehalf*t15*((1.0 + tanh(((-C4) +
        Rgb)/C5)))))))))));

    return (yy);
}

static double
a2long(double q, double L, double b, double p1, double p2, double q0)
{
    double C;
    const double onehalf = 1.0/2.0;

    if( L/b > 10.0) {
        C = 3.06/pow((L/b),0.44);
    } else {
        C = 1.0;
    }

    const double C1 = 1.22;
    const double C2 = 0.4288;
    const double C3 = -1.651;
    const double C4 = 1.523;
    const double C5 = 0.1477;
    const double miu = 0.585;

    const double Rg2 = Rgsquare(q,L,b);
    const double Rg22 = Rg2*Rg2;
    const double b2 = b*b;
    const double b3 = b*b*b;
    const double b4 = b3*b;
    const double q02 = q0*q0;
    const double q03 = q0*q0*q0;
    const double q04 = q03*q0;
    const double q05 = q04*q0;
    const double Rg = sqrt(Rg2);
    const double Rgb = Rg*q0/b;
    const double Rg02 = Rg2*q02;

    const double t1 = (1.0/(b* p1*pow(q0,((-1.0) - p1 - p2)) -
         b*p2*pow(q0,((-1.0) - p1 - p2)) ));

    const double t2 = (b*C*(((-1.0*((14.0*b3)/(15.0*q03*Rg2))) +
         (14.0*b3*pow((double)M_E,(-(Rg02/b2))))/(15.0*q03*Rg2) +
         (2.0*pow((double)M_E,(-(Rg02/b2)))*q0*((11.0/15.0 +
         (7*b2)/(15.0*Rg02)))*Rg2)/b)))/L;

    const double t3 = (Rg*((C3*pow(((Rgb)),((-3.0)/miu)) +
         C2*pow(((Rgb)),((-2.0)/miu)) +
         C1*pow(((Rgb)),((-1.0)/miu))))*
         pow(sech_WR(((-C4) +Rgb)/C5),2.0))/(2.0*C5);

    const double t4 = (b4*Rg*(((-1.0) + pow((double)M_E,(-(Rg02/b2))) +
         Rg02/b2))*pow(sech_WR(((-C4) +
         Rgb)/C5),2))/(C5*q04*Rg22);

    const double t5 = (2.0*b4*(((2.0*q0*Rg2)/b -
         (2.0*pow((double)M_E,(-(Rg02/b2)))*q0*Rg2)/b))*
         ((1.0 + onehalf*(((-1.0) - tanh(((-C4) +
         Rgb)/C5))))))/(q04*Rg22);

    const double t6 = (8.0*b4*b*(((-1.0) + pow((double)M_E,(-(Rg02/b2))) +
         Rg02/b2))*((1.0 + onehalf*(((-1) - tanh(((-C4) +
         Rgb)/C5))))))/(q05*Rg22);

    const double t7 = (((-((3.0*C3*Rg*pow(((Rgb)),((-1.0) -
         3.0/miu)))/miu)) - (2.0*C2*Rg*pow(((Rgb)),
         ((-1.0) - 2.0/miu)))/miu - (C1*Rg*pow(((Rgb)),
         ((-1.0) - 1.0/miu)))/miu));

    const double t8 = ((1.0 + tanh(((-C4) + Rgb)/C5)));

    const double t9 = (b*C*((4.0/15.0 - pow((double)M_E,(-(Rg02/b2)))*((11.0/15.0 +
         (7.0*b2)/(15*Rg02))) + (7.0*b2)/(15.0*Rg02))))/L;

    const double t10 = (2.0*b4*(((-1) + pow((double)M_E,(-(Rg02/b2))) +
          Rg02/b2))*((1.0 + onehalf*(((-1) - tanh(((-C4) +
          Rgb)/C5))))))/(q04*Rg22);

    const double yy = ((-1.0*(t1* ((-pow(q0,-p1)*(((b2*M_PI)/(L*q02) +
         t2 + t3 - t4 + t5 - t6 + onehalf*t7*t8)) - b*p1*pow(q0,((-1.0) - p1))*
         (((-((b*M_PI)/(L*q0))) + t9 + t10 +
         onehalf*((C3*pow(((Rgb)),((-3.0)/miu)) +
         C2*pow(((Rgb)),((-2.0)/miu)) +
         C1*pow(((Rgb)),((-1.0)/miu))))*
         ((1.0 + tanh(((-C4) + Rgb)/C5))))))))));

    return (yy);
}

//
static double
a_short(double q, double L, double b, double p1short,
        double p2short, double factor, double pdiff, double q0)
{
    const double Rg2_sh = Rgsquareshort(q,L,b);
    const double Rg2_sh2 = Rg2_sh*Rg2_sh;
    const double b3 = b*b*b;
    const double t1 = ((q0*q0*Rg2_sh)/(b*b));
    const double Et1 = exp(t1);
    const double Emt1 = 1.0/Et1;
    const double q02 = q0*q0;
    const double q0p = pow(q0,(-4.0 + p1short) );

    double yy = ((factor/(L*(pdiff)*Rg2_sh2)*
        ((b*Emt1*q0p*((8.0*b3*L - 8.0*b3*Et1*L - 2.0*b3*L*p2short +
        2.0*b3*Et1*L*p2short + 4.0*b*L*q02*Rg2_sh + 4.0*b*Et1*L*q02*Rg2_sh -
        2.0*b*Et1*L*p2short*q02*Rg2_sh - Et1*M_PI*q02*q0*Rg2_sh2 +
        Et1*p2short*M_PI*q02*q0*Rg2_sh2))))));

    return(yy);
}
static double
a1short(double q, double L, double b, double p1short, double p2short, double q0)
{

    double factor = 1.0;
    return a_short(q, L, b, p1short, p2short, factor, p1short - p2short, q0);
}

static double
a2short(double q, double L, double b, double p1short, double p2short, double q0)
{
    double factor = -1.0;
    return a_short(q, L, b, p2short, p1short, factor, p1short-p2short, q0);
}

//WR named this w (too generic)
static double
w_WR(double x)
{
    return 0.5*(1 + tanh((x - 1.523)/0.1477));
}

//
static double
u1(double q, double L, double b)
{
    return Rgsquareshort(q,L,b)*q*q;
}

static double
u_WR(double q, double L, double b)
{
    return Rgsquare(q,L,b)*q*q;
}

static double
Sdebye_kernel(double arg)
{
    // ORIGINAL
    double result = 2.0*(exp(-arg) + arg -1.0)/(arg*arg);

    // CONVERSION 1 from http://herbie.uwplse.org/
    //
    // exhibits discontinuity - needs more investigation
    //double a1 = 1.0/6.0;
    //double a2 = 1.0/72.0;
    //double a3 = 1.0/24.0;
    //double result = pow((1.0 - a1*arg - (a2+a3)*arg*arg), 2);

    return result;
}
static double
Sdebye(double q, double L, double b)
{
    double arg = u_WR(q,L,b);
    return Sdebye_kernel(arg);
}

//
static double
Sdebye1(double q, double L, double b)
{
    double arg = u1(q,L,b);
    return Sdebye_kernel(arg);

}

//
static double
Sexv(double q, double L, double b)
{

    const double C1=1.22;
    const double C2=0.4288;
    const double C3=-1.651;
    const double miu = 0.585;
    const double qRg = q*sqrt(Rgsquare(q,L,b));
    const double x = pow(qRg, -1.0/miu);

    double yy = (1.0 - w_WR(qRg))*Sdebye(q,L,b) +
            w_WR(qRg)*x*(C1 + x*(C2 + x*C3));
    return (yy);
}


static double
Sexvnew(double q, double L, double b)
{
    double yy;

    const double C1 =1.22;
    const double C2 =0.4288;
    const double C3 =-1.651;
    const double miu = 0.585;
    const double del=1.05;
    const double qRg = q*sqrt(Rgsquare(q,L,b));
    const double x = pow(qRg, -1.0/miu);


    //calculating the derivative to decide on the corection (cutoff) term?
    // I have modified this from WRs original code
    const double qdel = (Sexv(q*del,L,b)-Sexv(q,L,b))/(q*del - q);
    const double C_star2 =(qdel >= 0.0) ? 0.0 : 1.0;

    yy = (1.0 - w_WR(qRg))*Sdebye(q,L,b) +
         C_star2*w_WR(qRg)*
         x*(C1 + x*(C2 + x*C3));

    return (yy);
}

double Sk_WR(double q, double L, double b)
{
    //
    const double p1 = 4.12;
    const double p2 = 4.42;
    const double p1short = 5.36;
    const double p2short = 5.62;
    const double q0 = 3.1;

    double q0short = fmax(1.9/sqrt(Rgsquareshort(q,L,b)),3.0);
    double Sexvmodify, ans;

    const double C =(L/b > 10.0) ? 3.06/pow((L/b),0.44) : 1.0;

    if( L > 4*b ) { // Longer Chains
       if (q*b <= 3.1) {   //Modified by Yun on Oct. 15,
         Sexvmodify = Sexvnew(q, L, b);
         ans = Sexvmodify + C * (4.0/15.0 + 7.0/(15.0*u_WR(q,L,b)) -
            (11.0/15.0 + 7.0/(15.0*u_WR(q,L,b)))*exp(-u_WR(q,L,b)))*(b/L);

       } else { //q(i)*b > 3.1
         ans = a1long(q, L, b, p1, p2, q0)/(pow((q*b),p1)) +
               a2long(q, L, b, p1, p2, q0)/(pow((q*b),p2)) + M_PI/(q*L);
       }
    } else { //L <= 4*b Shorter Chains
       if (q*b <= fmax(1.9/sqrt(Rgsquareshort(q,L,b)),3.0) ) {
         if (q*b<=0.01) {
            ans = 1.0 - Rgsquareshort(q,L,b)*(q*q)/3.0;
         } else {
            ans = Sdebye1(q,L,b);
         }
       } else {  //q*b > max(1.9/sqrt(Rgsquareshort(q(i),L,b)),3)
         ans = a1short(q,L,b,p1short,p2short,q0short)/(pow((q*b),p1short)) +
               a2short(q,L,b,p1short,p2short,q0short)/(pow((q*b),p2short)) +
               M_PI/(q*L);
       }
    }

  return(ans);
}
