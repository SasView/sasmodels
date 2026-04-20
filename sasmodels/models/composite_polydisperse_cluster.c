#include <math.h>
#include <stdio.h>

/* GAMMLN ------------------------------------------------------------ */
double GAMMLN(double xx)
{
    static double cof[6] = {
        76.18009173, -86.50532033, 24.01409822,
        -1.231739516, 0.00120858003, -0.0000536382
    };
    double stp = 2.50662827465;
    double half = 0.5, one = 1.0, fpf = 5.5;
    double x = xx - one;
    double tmp = x + fpf;
    tmp = (x + half) * log(fabs(tmp)) - tmp;

    double ser = one;
    for(int j=0;j<6;j++){
        x += one;
        ser += cof[j]/x;
    }
    return tmp + log(stp * ser);
}

/* SCHULZ ------------------------------------------------------------ */
double SCHULZ(double R, double RA, double Z)
{
    double dum =
        (Z+1.0)*log(fabs((Z+1.0)/RA))
      + Z*log(fabs(R))
      - (Z+1.0)*R/RA
      - GAMMLN(Z+1.0);

    return exp(dum);
}

/* FI ---------------------------------------------------------------- */
double FI(double X)
{
    if(X > 0.05)
        return 3.0*(sin(X)-X*cos(X)) / (X*X*X);
    else
        return 1.0 - 0.1*X*X;
}

/* S_HS -------------------------------------------------------------- */
double S_HS(double Q, double RHS, double ETA)
{
    double ALN = pow(1.0 - ETA, 4.0);
    double AL = pow(1.0 + 2.0*ETA, 2.0) / ALN;
    double BE = -6.0 * ETA * pow(1.0 + 0.5*ETA, 2.0) / ALN;
    double GA = 0.5 * ETA * AL;

    double AR = 2.0 * RHS * Q;
    double GG;

    if (AR < 0.4)
    {
        GG =  AL*(1.0/3.0 - AR*AR/30.0)
            + BE*(1.0/4.0 - AR*AR/36.0)
            + GA*(1.0/6.0 - AR*AR/48.0);
    }
    else
    {
        double SA = sin(AR);
        double CA = cos(AR);

        GG =  AL*(SA - AR*CA)/pow(AR,3.0)
            + BE*(2.0*AR*SA + (2.0 - AR*AR)*CA - 2.0)/pow(AR,4.0)
            + GA*(-pow(AR,4.0)*CA
                 + 4.0*((3.0*AR*AR - 6.0)*CA
                      + (pow(AR,3.0) - 6.0*AR)*SA + 6.0)) / pow(AR,6.0);
    }

    return 1.0 / (1.0 + 24.0*ETA*GG);
}
/* PSCHULZ ----------------------------------------------------------- */
double PSCHULZ(double Q, double RAV, double SIGMA)
{
    double Z = 1.0/(SIGMA*SIGMA) - 1.0;
    double pi = acos(-1.0);
    double FIPI3 = 4.0*pi/3.0;

    double ALP  = (Z+1.0)/(Q*RAV);
    double ALP1 = (Q*RAV)/(Z+1.0);
    double ATTA = atan(2.0*ALP1);

    double PSI1 = (Z+1.0)*ATTA;
    double PSI2 = (Z+2.0)*ATTA;
    double PSI3 = (Z+3.0)*ATTA;

    double ALP2   = ALP*ALP;
    double ALP2_4 = 1.0/(4.0 + ALP2);

    double C1 = pow(ALP2*ALP2_4,0.5*Z+0.5)*pow(ALP/(Z+1.0),6.0);
    double T1 = pow(ALP/(Z+1.0),6.0) - C1*cos(PSI1);

    double C2 = pow(ALP2*ALP2_4,0.5*Z+1.5)*pow(ALP/(Z+1.0),4.0);
    double T2 = (Z+2.0)/(Z+1.0)*( pow(ALP/(Z+1.0),4.0) + C2*cos(PSI3) );

    double C3 = pow(ALP2*ALP2_4,0.5*Z+1.0)*pow(ALP/(Z+1.0),5.0);
    double T3 = -2.0*C3*sin(PSI2);

    double V2 = pow(FIPI3,2.0) *(Z+6.0)*(Z+5.0)*(Z+4.0)*
                (Z+3.0)*(Z+2.0)/pow(Z+1.0,5.0);

    return (8.0*pi*pi*(T1+T2+T3))/V2;
}

/* P_POLY ------------------------------------------------------------ */
double P_POLY(double Q, double RL, double SIGL, double d, double N_agg)
{
    double PI43 = 4.1887892;
    double RS = 0.5 * d;
    double RDIF = RL - RS;

    double Z = 1.0/(SIGL*SIGL) - 1.0;

    double SVL  = PI43*pow(RDIF,3.0)*(Z+3.0)*(Z+2.0)/((Z+1.0)*(Z+1.0));
    double SVL2 = PI43*PI43 * pow(RDIF,6.0)
                  * (Z+6.0)*(Z+5.0)*(Z+4.0)*(Z+3.0)*(Z+2.0)
                  / pow(Z+1.0,5.0);

    double VS   = PI43*pow(RS,3.0);
    double ETA = N_agg * (SVL * VS)/SVL2;
    double ALPP = ETA/VS;

    

    int NPOI = 100;
    double RMAX = RS + 6.0*RL*SIGL;
    double DELR = RMAX/NPOI;

    double SSQ = 0.0;
    double SSQN= 0.0;

    for(int i=1;i<=NPOI;i++)
    {
        double R  = RS + (i-0.5)*DELR;
        double RD = R - RS;
        double D  = SCHULZ(RD, RDIF, Z);
        
        double RH = R - 2.0*RS; 
        if(RH < 0) RH = 0.0;
        
        double FR = (pow(RD,3.0) - pow(RH,3.0))/pow(RD,3.0);
    if(FR < 0.0) FR = 0.0;
    if(FR > 1.0) FR = 1.0;

    double ETAEFF = ETA*(FR*0.5 + (1.0-FR));

        SSQ  += D * pow(RD,3.0) * S_HS(Q,RS,ETAEFF);
        SSQN += D * pow(RD,3.0) + 1.0e-10;
    }

    SSQ = SSQ/SSQN;

    double PS   = VS*VS * FI(Q*RS)*FI(Q*RS);
    double PLAV = PSCHULZ(Q, RDIF, SIGL);

    double RN  = SVL  * ALPP;
    double RNX = sqrt(SVL2)*ALPP;

    return (RN*PS*SSQ + PS*RNX*(RNX-1.0)*PLAV)/(RN*PS);
}

/* ----------- SASVIEW ENTRY POINT ----------- */
double SQ_COMP(double Q, double RL, double SIGL, double d, double N_agg)
{
    return P_POLY(Q, RL, SIGL, d, N_agg);
}

/* SasView required function */
double Iq(double Q, double R_clust, double sig_rel_R, double dist_points, double N_agg)
{
    double intensity;
    
    intensity = SQ_COMP(Q, R_clust, sig_rel_R, dist_points, N_agg);
    
    return intensity;
}