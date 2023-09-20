/*Lamellar_ParaCrystal - July 2019 RKH - hack previous sasview version to become Kotlarchyk & Ritzau 
J.Appl.Cryst. 24(1991)753-758 with some "corrections" as per version from FISH .
SHOULD put paraCryst_kr and rlayff_kr into a library as they are also used by the original lamellar_stack_paracrystal_kr.c
*/
double paraCryst_kr(double qval, double RM, double D, double GD);
double rlayff_kr(double qval, double RLBAR, double GL, double TT, double RSIG, double f1);

static double
Iq(double qval,
   double thickness,
   double sigma_t_by_t,
   double interface_t,
   double rsig_lorentz,
   double scale1,
   double scale2,
   double scale3,
   double scale4,
   double scale5,
   double davg,
   double sigma_d_by_d,
   double sld,
   double solvent_sld)
{
      double f1, Znq2,Znq3,Znq4,Znq5;
      const double fp_2 = 2.0;
      const double fp_3 = 3.0;
      const double fp_4 = 4.0;
      const double fp_5 = 5.0;
      
/*      //get the fractional part of Nlayers, to determine the "mixing" of N's
// have removed this averaging of two values of Nlayers, as the pattern changes 
// dramatically if e.g. go from 3 to 3.5 to 4 layers  by this method, 
// the original varies more smoothly, even if a non-integer number of layers is "unphysical"!
      int n1 = (int)(fp_Nlayers);            //truncate towards zero
      int n2 = n1 + 1;
      const double xn = (double)n2 - fp_Nlayers;      //fractional contribution of n1
      //
      //calculate the n1 contribution

      Znq = xn*paraCryst_kr(qval,(double)n1,davg,sigma_d_by_d);
      //
      //calculate the n2 contribution
      Znq += (1.0-xn)*paraCryst_kr(qval,(double)n2,davg,sigma_d_by_d);
*/      
       
      Znq2 = paraCryst_kr(qval,fp_2,davg,sigma_d_by_d);
      Znq3 = paraCryst_kr(qval,fp_3,davg,sigma_d_by_d);
      Znq4 = paraCryst_kr(qval,fp_4,davg,sigma_d_by_d);
      Znq5 = paraCryst_kr(qval,fp_5,davg,sigma_d_by_d);
      //
      // will one day separate the S(Q) and P(Q) here into individual models
      // meanwhile it should be relatively simple to make other versions with shell/core/shell layers etc.
      // I believe that <F> and F^2> are being properly averaged for the polydisperse layer
      const double f2 = rlayff_kr(qval,thickness,sigma_t_by_t,interface_t,rsig_lorentz,f1);
      //
      const double contr = sld - solvent_sld;
      // BEWARE may need to rescale f1 also for sasview 5
      // do not know yet why need a rescale of 0.25 here to get same scale parameter as FISH
      // need check against lamellar model!
      // 12June2020 add 5 terms 
      const double inten = 0.25*contr*contr*(scale1 + scale2*Znq2 + scale3*Znq3 + scale4*Znq4 + scale5*Znq5)*f2;
      
}

// functions for the lamellar paracrystal model, hacked from FISH fortran code, July 2019, RKH
double 
paraCryst_kr(double QQ, double RM, double D, double GD) {
    // presume RMU is mu=1 = cos(angle between Q and axis of stack) in eqn (13)
    const double RMU = 1.0;
    const double QDM = QQ*D*RMU;
    const double AA = -square(GD*QDM)*0.5;
    //      F=EXPSPEC(-0.5*(GD*QDM)**2)
    //      F2=F*F
    //      FPM=F**RM
    //  take the powers of F inside the exponential
    const double F = exp(AA);
    const double F2 = exp(2.0*AA);
    const double FPM = exp(RM*AA);
    //
    // this hits severe problems when Q goes small and F goes up towards
    //  1.0 from below.  We may have to argue that the cos(M*QDM) oscillates
    // faster than cos(QDM) so averages to zero ?
    //
    //      PARCRYS = (1.0-F2)/(1.0-2.*F*COS(QDM) +F2) +(1.0/RM)*
    //     > ( -2.0*F*( (1.0+F2)*COS(QDM) -2.0*F -FPM*COS(QDM*(RM+1.)) +
    //     > 2.0*FPM*F*COS(RM*QDM) -FPM*F2*COS( (RM-1.0)*QDM ) ) )/
    //     > ( 1.0 - 2.0*F*COS(QDM) + F2)**2
    //
    // rearranged, 13/2/92, gives same results
    // tested 14/2/92 in double precision, single prec. gave same results except
    // small errors at very low Q ( <.002), and through peak (<.00002),
    // with RM=10, D=100, GD=0.15 as typical parameters.
    //
    double Snq = (1.0/RM)*(4.0*F2 + RM*(1. - F2*F2) - 4.0*F2*FPM*cos(QDM*RM) + 
                  cos(QDM)*2.*F*(-RM*(1. - F2) - (1. + F2) + FPM*(1. + F2)*cos(RM*QDM)) +
                  sin(QDM)*2.*F*FPM*(F2 - 1.)*sin(RM*QDM))/square(1. - 2.0*F*cos(QDM) + F2);
    // for very low GD, and always at low Q the results tend to the 
    // equation for a perfectly ordered system:
    //      PARCRYS=(1.0/RM)*(SIN(0.5*RM*QDM)/SIN(0.5*QDM))**2
    return Snq;
}

double
rlayff_kr(double QQ, double RLBAR, double GL, double TT, double RSIG, double f1) {
    //// RKH July 2019, FORM FACTOR FOR POLYDISPERSE SHEET - Kotlarchyk & Ritzau, 
    ////  J.Appl.Cryst. 24(1991)753-758,  eqns (16) & (17) with 1/Q**2
    // included as per eqn (13)          12/2/92 RKH
    // interfacial thickness TT = (2*pi)**0.5 * SIG    CORRECTED 7/4/92
    // GL = SIGMA(L)/LBAR = (Z+1)**-0.5
    // 16/3/93 add extra parameter RSIG to model, as Lorentz factor,
    // see Skipper et.al. J.Chem.Phys 94(1991)5751-5760
        double ZP1 = 1.0/square(GL);
        const double Z = ZP1 - 1.0;
        double ZM1 = Z - 1.0;
        const double AL = ZP1/(QQ*RLBAR);
        const double AL2 = 2.0*AL;
        const double DEN = (1.0 + square(QQ*RSIG)*0.5);
        const double Q2S2 = -square(QQ*TT)/(2.0*M_PI);
        //
        //     F=(AL2**ZP1)*SIN(Z*atan(1.0/AL2))*EXP(0.5*Q2S2)/
        //    >  (Z* ( ( AL2**2 +1)**(0.5*Z) )*QQ )
        //rearrange in attempt to avoid overflows ( note AL can be 100, Z say 25 )
        
        // sasview 4 does not need <F> but it is here ready .
        // I believe that <F> and F^2> are being properly averaged for the polydisperse layer
        double AL2SQ = square(AL2);
        // f1 = (1. + 1.0/AL2SQ)**(-0.5*ZP1)*sqrt(AL2SQ + 1)*sin(Z*atan(1.0/AL2))*exp(0.5*Q2S2)/Z;
        f1 = pow( 1.0 + 1.0/AL2SQ, -0.5*ZP1);
        f1 *= sqrt(AL2SQ + 1)*sin(Z*atan(1.0/AL2))*exp(0.5*Q2S2)/Z;
        f1 /= sqrt(DEN);
        
        //
        //     RLAYFF= (AL2**ZP1)*( (AL2**(1.0-Z)) -((AL2**2 +4.0)**(-0.5*ZM1))*
        //    >   COS( ZM1*atan(1.0/AL) ) )*EXP(Q2S2)/ (2.0*Z*ZM1*QQ*QQ)
        
        //rearrange to avoid overflows ( note AL can be 100, Z say 25 )
        //included an extra factor of 4.0 to get beta(Q) to go to 1.0 at low
        // Q and P(Q) to go to 1.0/Q**2
        // f2 = 4.0*(AL2SQ)*(1. - ((1. + 4.0/(AL2SQ))**(-0.5*ZM1))*cos(ZM1*atan(1.0/AL)))*exp(Q2S2)/(2.0*Z*ZM1);
        double f2 = 4.0*(AL2SQ)*(1. - pow(1. + 4.0/AL2SQ,-0.5*ZM1)*cos(ZM1*atan(1.0/AL))) * exp(Q2S2)/(2.0*Z*ZM1);
        f2 /= DEN;
        return f2;
}

