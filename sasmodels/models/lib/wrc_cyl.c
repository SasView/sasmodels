/*
    Functions for WRC implementation of flexible cylinders
*/
static double
AlphaSquare(double x)
{
    return pow( (1.0 + (x/3.12)*(x/3.12) +
         (x/8.67)*(x/8.67)*(x/8.67)),(0.176/3.0) );
}

//
static double
Rgsquarezero(double q, double L, double b)
{
    return (L*b/6.0) * (1.0 - 1.5*(b/L) + 1.5*pow((b/L),2) -
         0.75*pow((b/L),3)*(1.0 - exp(-2.0*(L/b))));
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

static double
sech_WR(double x)
{
    return(1/cosh(x));
}

static double
a1long(double q, double L, double b, double p1, double p2, double q0)
{
    double yy,C,C1,C2,C3,C4,C5,miu,Rg2;
    double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15;
    double E,pi;
    double b2,b3,b4,q02,q03,q04,q05,Rg22;

    pi = 4.0*atan(1.0);
    E = 2.718281828459045091;

    if( L/b > 10.0) {
        C = 3.06/pow((L/b),0.44);
    } else {
        C = 1.0;
    }

    C1 = 1.22;
    C2 = 0.4288;
    C3 = -1.651;
    C4 = 1.523;
    C5 = 0.1477;
    miu = 0.585;

    Rg2 = Rgsquare(q,L,b);
    Rg22 = Rg2*Rg2;
    b2 = b*b;
    b3 = b*b*b;
    b4 = b3*b;
    q02 = q0*q0;
    q03 = q0*q0*q0;
    q04 = q03*q0;
    q05 = q04*q0;

    t1 = (b*C*((4.0/15.0 - pow(E,(-((q02*Rg2)/b2)))*((11.0/15.0 +
        (7.0*b2)/(15.0*q02*Rg2))) + (7.0*b2)/(15.0*q02*Rg2))));

    t2 = (2.0*b4*(((-1.0) + pow(E,(-((q02*Rg2)/b2))) +
        (q02*Rg2)/b2))*((1.0 + 1.0/2.0*(((-1.0) -
        tanh(((-C4) + (sqrt(Rg2)*q0)/b)/C5))))));

    t3 = ((C3*pow((((sqrt(Rg2)*q0)/b)),((-3.0)/miu)) +
        C2*pow((((sqrt(Rg2)*q0)/b)),((-2.0)/miu)) +
        C1*pow((((sqrt(Rg2)*q0)/b)),((-1.0)/miu))));

    t4 = ((1.0 + tanh(((-C4) + (sqrt(Rg2)*q0)/b)/C5)));

    t5 = (1.0/(b*p1*pow(q0,((-1.0) - p1 - p2)) -
        b*p2*pow(q0,((-1.0) - p1 - p2))));

    t6 = (b*C*(((-((14.0*b3)/(15.0*q03*Rg2))) +
        (14.0*b3*pow(E,(-((q02*Rg2)/b2))))/(15.0*q03*Rg2) +
        (2.0*pow(E,(-((q02*Rg2)/b2)))*q0*((11.0/15.0 +
        (7.0*b2)/(15.0*q02*Rg2)))*Rg2)/b)));

    t7 = (sqrt(Rg2)*((C3*pow((((sqrt(Rg2)*q0)/b)),((-3.0)/miu)) +
        C2*pow((((sqrt(Rg2)*q0)/b)),((-2.0)/miu)) +
        C1*pow((((sqrt(Rg2)*q0)/b)),((-1.0)/miu))))*pow(sech_WR(((-C4) +
        (sqrt(Rg2)*q0)/b)/C5),2));

    t8 = (b4*sqrt(Rg2)*(((-1.0) + pow(E,(-((q02*Rg2)/b2))) +
        (q02*Rg2)/b2))*pow(sech_WR(((-C4) + (sqrt(Rg2)*q0)/b)/C5),2));

    t9 = (2.0*b4*(((2.0*q0*Rg2)/b -
        (2.0*pow(E,(-((q02*Rg2)/b2)))*q0*Rg2)/b))*((1.0 + 1.0/2.0*(((-1.0) -
        tanh(((-C4) + (sqrt(Rg2)*q0)/b)/C5))))));

    t10 = (8.0*b4*b*(((-1.0) + pow(E,(-((q02*Rg2)/b2))) +
        (q02*Rg2)/b2))*((1.0 + 1.0/2.0*(((-1.0) - tanh(((-C4) +
        (sqrt(Rg2)*q0)/b)/C5))))));

    t11 = (((-((3.0*C3*sqrt(Rg2)*pow((((sqrt(Rg2)*q0)/b)),((-1.0) -
        3.0/miu)))/miu)) - (2.0*C2*sqrt(Rg2)*pow((((sqrt(Rg2)*q0)/b)),((-1.0) -
        2.0/miu)))/miu - (C1*sqrt(Rg2)*pow((((sqrt(Rg2)*q0)/b)),((-1.0) -
        1.0/miu)))/miu));

    t12 = ((1.0 + tanh(((-C4) + (sqrt(Rg2)*q0)/b)/C5)));

    t13 = (b*C*((4.0/15.0 - pow(E,(-((q02*Rg2)/b2)))*((11.0/15.0 +
        (7.0*b2)/(15.0*q02* Rg2))) +
        (7.0*b2)/(15.0*q02*Rg2))));

    t14 = (2.0*b4*(((-1.0) + pow(E,(-((q02*Rg2)/b2))) +
        (q02*Rg2)/b2))*((1.0 + 1.0/2.0*(((-1.0) - tanh(((-C4) +
        (sqrt(Rg2)*q0)/b)/C5))))));

    t15 = ((C3*pow((((sqrt(Rg2)*q0)/b)),((-3.0)/miu)) +
        C2*pow((((sqrt(Rg2)*q0)/b)),((-2.0)/miu)) +
        C1*pow((((sqrt(Rg2)*q0)/b)),((-1.0)/miu))));


    yy = (pow(q0,p1)*(((-((b*pi)/(L*q0))) +t1/L +t2/(q04*Rg22) +
        1.0/2.0*t3*t4)) + (t5*((pow(q0,(p1 - p2))*
        (((-pow(q0,(-p1)))*(((b2*pi)/(L*q02) +t6/L +t7/(2.0*C5) -
        t8/(C5*q04*Rg22) + t9/(q04*Rg22) -t10/(q05*Rg22) + 1.0/2.0*t11*t12)) -
        b*p1*pow(q0,((-1.0) - p1))*(((-((b*pi)/(L*q0))) + t13/L +
        t14/(q04*Rg22) + 1.0/2.0*t15*((1.0 + tanh(((-C4) +
        (sqrt(Rg2)*q0)/b)/C5)))))))))));

    return (yy);
}

static double
a2long(double q, double L, double b, double p1, double p2, double q0)
{
    double yy,C1,C2,C3,C4,C5,miu,C,Rg2;
    double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,pi;
    double E,b2,b3,b4,q02,q03,q04,q05,Rg22;

    pi = 4.0*atan(1.0);
    E = 2.718281828459045091;
    if( L/b > 10.0) {
        C = 3.06/pow((L/b),0.44);
    } else {
        C = 1.0;
    }

    C1 = 1.22;
    C2 = 0.4288;
    C3 = -1.651;
    C4 = 1.523;
    C5 = 0.1477;
    miu = 0.585;

    Rg2 = Rgsquare(q,L,b);
    Rg22 = Rg2*Rg2;
    b2 = b*b;
    b3 = b*b*b;
    b4 = b3*b;
    q02 = q0*q0;
    q03 = q0*q0*q0;
    q04 = q03*q0;
    q05 = q04*q0;

    t1 = (1.0/(b* p1*pow(q0,((-1.0) - p1 - p2)) -
        b*p2*pow(q0,((-1.0) - p1 - p2)) ));

    t2 = (b*C*(((-1.0*((14.0*b3)/(15.0*q03*Rg2))) +
        (14.0*b3*pow(E,(-((q02*Rg2)/b2))))/(15.0*q03*Rg2) +
        (2.0*pow(E,(-((q02*Rg2)/b2)))*q0*((11.0/15.0 +
        (7*b2)/(15.0*q02*Rg2)))*Rg2)/b)))/L;

    t3 = (sqrt(Rg2)*((C3*pow((((sqrt(Rg2)*q0)/b)),((-3.0)/miu)) +
        C2*pow((((sqrt(Rg2)*q0)/b)),((-2.0)/miu)) +
        C1*pow((((sqrt(Rg2)*q0)/b)),((-1.0)/miu))))*
        pow(sech_WR(((-C4) +(sqrt(Rg2)*q0)/b)/C5),2.0))/(2.0*C5);

    t4 = (b4*sqrt(Rg2)*(((-1.0) + pow(E,(-((q02*Rg2)/b2))) +
        (q02*Rg2)/b2))*pow(sech_WR(((-C4) +
        (sqrt(Rg2)*q0)/b)/C5),2))/(C5*q04*Rg22);

    t5 = (2.0*b4*(((2.0*q0*Rg2)/b -
        (2.0*pow(E,(-((q02*Rg2)/b2)))*q0*Rg2)/b))*
        ((1.0 + 1.0/2.0*(((-1.0) - tanh(((-C4) +
        (sqrt(Rg2)*q0)/b)/C5))))))/(q04*Rg22);

    t6 = (8.0*b4*b*(((-1.0) + pow(E,(-((q02*Rg2)/b2))) +
        (q02*Rg2)/b2))*((1.0 + 1.0/2.0*(((-1) - tanh(((-C4) +
        (sqrt(Rg2)*q0)/b)/C5))))))/(q05*Rg22);

    t7 = (((-((3.0*C3*sqrt(Rg2)*pow((((sqrt(Rg2)*q0)/b)),((-1.0) -
        3.0/miu)))/miu)) - (2.0*C2*sqrt(Rg2)*pow((((sqrt(Rg2)*q0)/b)),
        ((-1.0) - 2.0/miu)))/miu - (C1*sqrt(Rg2)*pow((((sqrt(Rg2)*q0)/b)),
        ((-1.0) - 1.0/miu)))/miu));

    t8 = ((1.0 + tanh(((-C4) + (sqrt(Rg2)*q0)/b)/C5)));

    t9 = (b*C*((4.0/15.0 - pow(E,(-((q02*Rg2)/b2)))*((11.0/15.0 +
        (7.0    *b2)/(15*q02*Rg2))) + (7.0*b2)/(15.0*q02*Rg2))))/L;

    t10 = (2.0*b4*(((-1) + pow(E,(-((q02*Rg2)/b2))) +
        (q02*Rg2)/b2))*((1.0 + 1.0/2.0*(((-1) - tanh(((-C4) +
        (sqrt(Rg2)*q0)/b)/C5))))))/(q04*Rg22);

    yy = ((-1.0*(t1* ((-pow(q0,-p1)*(((b2*pi)/(L*q02) +
        t2 + t3 - t4 + t5 - t6 + 1.0/2.0*t7*t8)) - b*p1*pow(q0,((-1.0) - p1))*
        (((-((b*pi)/(L*q0))) + t9 + t10 +
        1.0/2.0*((C3*pow((((sqrt(Rg2)*q0)/b)),((-3.0)/miu)) +
        C2*pow((((sqrt(Rg2)*q0)/b)),((-2.0)/miu)) +
        C1*pow((((sqrt(Rg2)*q0)/b)),((-1.0)/miu))))*
        ((1.0 + tanh(((-C4) + (sqrt(Rg2)*q0)/b)/C5))))))))));

    return (yy);
}

//
static double
a1short(double q, double L, double b, double p1short, double p2short, double q0)
{
    double yy,Rg2_sh;
    double t1,E,Rg2_sh2,Et1,Emt1,q02,q0p,b3;
    double pi;

    E = 2.718281828459045091;
    pi = 4.0*atan(1.0);
    Rg2_sh = Rgsquareshort(q,L,b);
    Rg2_sh2 = Rg2_sh*Rg2_sh;
    b3 = b*b*b;
    t1 = ((q0*q0*Rg2_sh)/(b*b));
    Et1 = pow(E,t1);
    Emt1 =pow(E,-t1);
    q02 = q0*q0;
    q0p = pow(q0,(-4.0 + p1short) );

    yy = ((1.0/(L*((p1short - p2short))*Rg2_sh2)*
        ((b*Emt1*q0p*((8.0*b3*L - 8.0*b3*Et1*L - 2.0*b3*L*p2short +
        2.0*b3*Et1*L*p2short + 4.0*b*L*q02*Rg2_sh + 4.0*b*Et1*L*q02*Rg2_sh -
        2.0*b*Et1*L*p2short*q02*Rg2_sh - Et1*pi*q02*q0*Rg2_sh2 +
        Et1*p2short*pi*q02*q0*Rg2_sh2))))));

    return(yy);
}

static double
a2short(double q, double L, double b, double p1short, double p2short, double q0)
{
    double yy,Rg2_sh;
    double t1,E,Rg2_sh2,Et1,Emt1,q02,q0p;
    double pi;

    E = 2.718281828459045091;
    pi = 4.0*atan(1.0);
    Rg2_sh = Rgsquareshort(q,L,b);
    Rg2_sh2 = Rg2_sh*Rg2_sh;
    t1 = ((q0*q0*Rg2_sh)/(b*b));
    Et1 = pow(E,t1);
    Emt1 =pow(E,-t1);
    q02 = q0*q0;
    q0p = pow(q0,(-4.0 + p2short) );

    //E is the number e
    yy = ((-(1.0/(L*((p1short - p2short))*Rg2_sh2)*
        ((b*Emt1*q0p*((8.0*b*b*b*L - 8.0*b*b*b*Et1*L - 2.0*b*b*b*L*p1short +
        2.0*b*b*b*Et1*L*p1short + 4.0*b*L*q02*Rg2_sh + 4.0*b*Et1*L*q02*Rg2_sh -
        2.0*b*Et1*L*p1short*q02*Rg2_sh - Et1*pi*q02*q0*Rg2_sh2 +
        Et1*p1short*pi*q02*q0*Rg2_sh2)))))));

    return (yy);
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
    double result = 2.0*(exp(-arg) + arg -1.0)/(pow((arg),2));

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
    double yy,C1,C2,C3,miu,Rg2;
    C1=1.22;
    C2=0.4288;
    C3=-1.651;
    miu = 0.585;

    Rg2 = Rgsquare(q,L,b);

    yy = (1.0 - w_WR(q*sqrt(Rg2)))*Sdebye(q,L,b) +
        w_WR(q*sqrt(Rg2))*(C1*pow((q*sqrt(Rg2)),(-1.0/miu)) +
        C2*pow((q*sqrt(Rg2)),(-2.0/miu)) +
        C3*pow((q*sqrt(Rg2)),(-3.0/miu)));

    return (yy);
}


static double
Sexvnew(double q, double L, double b)
{
    double yy,C1,C2,C3,miu;
    double del=1.05,C_star2,Rg2;

    C1=1.22;
    C2=0.4288;
    C3=-1.651;
    miu = 0.585;

    //calculating the derivative to decide on the corection (cutoff) term?
    // I have modified this from WRs original code

    if( (Sexv(q*del,L,b)-Sexv(q,L,b))/(q*del - q) >= 0.0 ) {
        C_star2 = 0.0;
    } else {
        C_star2 = 1.0;
    }

    Rg2 = Rgsquare(q,L,b);

    yy = (1.0 - w_WR(q*sqrt(Rg2)))*Sdebye(q,L,b) +
        C_star2*w_WR(q*sqrt(Rg2))*(C1*pow((q*sqrt(Rg2)),(-1.0/miu)) +
        C2*pow((q*sqrt(Rg2)),(-2.0/miu)) + C3*pow((q*sqrt(Rg2)),(-3.0/miu)));

    return (yy);
}

double Sk_WR(double q, double L, double b);
double Sk_WR(double q, double L, double b)
{
  //
  double p1,p2,p1short,p2short,q0;
  double C,ans,q0short,Sexvmodify,pi;

  pi = 4.0*atan(1.0);

  p1 = 4.12;
  p2 = 4.42;
  p1short = 5.36;
  p2short = 5.62;
  q0 = 3.1;

  q0short = fmax(1.9/sqrt(Rgsquareshort(q,L,b)),3.0);


  if(L/b > 10.0) {
    C = 3.06/pow((L/b),0.44);
  } else {
    C = 1.0;
  }
  //

  //
  if( L > 4*b ) { // Longer Chains
    if (q*b <= 3.1) {   //Modified by Yun on Oct. 15,

      Sexvmodify = Sexvnew(q, L, b);

      ans = Sexvmodify + C * (4.0/15.0 + 7.0/(15.0*u_WR(q,L,b)) -
            (11.0/15.0 + 7.0/(15.0*u_WR(q,L,b)))*exp(-u_WR(q,L,b)))*(b/L);

    } else { //q(i)*b > 3.1
      ans = a1long(q, L, b, p1, p2, q0)/(pow((q*b),p1)) +
            a2long(q, L, b, p1, p2, q0)/(pow((q*b),p2)) + pi/(q*L);
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
            pi/(q*L);
    }
  }

  return(ans);
}
