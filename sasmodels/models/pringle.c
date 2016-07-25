double form_volume(double radius,
          double thickness);

double Iq(double q,
          double radius,
          double thickness,
          double alpha,
          double beta,
          double sld_pringle,
          double sld_solvent);

static
double pringleC(double radius,
                double alpha,
                double beta,
                double q,
                double phi,
                double n) {

    double va, vb;
    double bessargs, cosarg, bessargcb;
    double r, retval, yyy;


    va = 0;
    vb = radius;

    // evaluate at Gauss points
    // remember to index from 0,size-1

    double summ = 0.0;		// initialize integral
    int ii = 0;
    do {
        // Using 76 Gauss points
        r = (Gauss76Z[ii] * (vb - va) + vb + va) / 2.0;

        bessargs = q*r*sin(phi);
        cosarg = q*r*r*alpha*cos(phi);
        bessargcb = q*r*r*beta*cos(phi);

        yyy = Gauss76Wt[ii]*r*cos(cosarg)
                *sas_JN(n, bessargcb)
                *sas_JN(2*n, bessargs);
        summ += yyy;

        ii += 1;
    } while (ii < N_POINTS_76);			// end of loop over quadrature points
    //
    // calculate value of integral to return

    retval = (vb - va) / 2.0 * summ;
    retval = retval / pow(r, 2.0);

    return retval;
}

static
double pringleS(double radius,
                double alpha,
                double beta,
                double q,
                double phi,
                double n) {

    double va, vb, summ;
    double bessargs, sinarg, bessargcb;
    double r, retval, yyy;
    // set up the integration
    // end points and weights

    va = 0;
    vb = radius;

    // evaluate at Gauss points
    // remember to index from 0,size-1

    summ = 0.0;		// initialize integral
    int ii = 0;
    do {
        // Using 76 Gauss points
        r = (Gauss76Z[ii] * (vb - va) + vb + va) / 2.0;

        bessargs = q*r*sin(phi);
        sinarg = q*r*r*alpha*cos(phi);
        bessargcb = q*r*r*beta*cos(phi);

        yyy = Gauss76Wt[ii]*r*sin(sinarg)
                    *sas_JN(n, bessargcb)
                    *sas_JN(2*n, bessargs);
        summ += yyy;

        ii += 1;
    } while (ii < N_POINTS_76);

    // end of loop over quadrature points
    //
    // calculate value of integral to return

    retval = (vb-va)/2.0*summ;
    retval = retval/pow(r, 2.0);

    return retval;
}

static
double _kernel(double thickness,
               double radius,
               double alpha,
               double beta,
               double q,
               double phi) {

    const double sincarg = q * thickness * cos(phi) / 2.0;
    const double sincterm = pow(sin(sincarg) / sincarg, 2.0);

    //calculate sum term from n = -3 to 3
    double sumterm = 0.0;
    for (int nn = -3; nn <= 3; nn++) {
        double powc = pringleC(radius, alpha, beta, q, phi, nn);
        double pows = pringleS(radius, alpha, beta, q, phi, nn);
        sumterm += pow(powc, 2.0) + pow(pows, 2.0);
    }
    double retval = 4.0 * sin(phi) * sumterm * sincterm;

    return retval;

}

static double pringles_kernel(double q,
          double radius,
          double thickness,
          double alpha,
          double beta,
          double sld_pringle,
          double sld_solvent)
{

    //upper and lower integration limits
    const double lolim = 0.0;
    const double uplim = M_PI / 2.0;

    double summ = 0.0;			//initialize integral

    double delrho = sld_pringle - sld_solvent; //make contrast term

    for (int i = 0; i < N_POINTS_76; i++) {
        double phi = (Gauss76Z[i] * (uplim - lolim) + uplim + lolim) / 2.0;
        summ += Gauss76Wt[i] * _kernel(thickness, radius, alpha, beta, q, phi);
    }

    double answer = (uplim - lolim) / 2.0 * summ;
    answer *= delrho*delrho;

    return answer;
}

double form_volume(double radius,
        double thickness){

        return 1.0;
}

double Iq(double q,
          double radius,
          double thickness,
          double alpha,
          double beta,
          double sld_pringle,
          double sld_solvent)
{
    const double form = pringles_kernel(q,
                  radius,
                  thickness,
                  alpha,
                  beta,
                  sld_pringle,
                  sld_solvent);

    return 1.0e-4*form*M_PI*radius*radius*thickness;
}
