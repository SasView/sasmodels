double form_volume(double thick_core,
                   double thick_layer,
                   double radius,
                   double n_stacking);

double Iq(double q,
          double thick_core,
          double thick_layer,
          double radius,
          double n_stacking,
          double sigma_dnn,
          double core_sld,
          double layer_sld,
          double solvent_sld);

double Iqxy(double qx, double qy,
          double thick_core,
          double thick_layer,
          double radius,
          double n_stacking,
          double sigma_dnn,
          double core_sld,
          double layer_sld,
          double solvent_sld,
          double theta,
          double phi);

static
double _kernel(double q,
               double radius,
               double core_sld,
               double layer_sld,
               double solvent_sld,
               double halfheight,
               double thick_layer,
               double sin_alpha,
               double cos_alpha,
               double sigma_dnn,
               double d,
               double n_stacking)

{
    // q is the q-value for the calculation (1/A)
    // radius is the core radius of the cylinder (A)
    // *_sld are the respective SLD's
    // halfheight is the *Half* CORE-LENGTH of the cylinder = L (A)
    // zi is the dummy variable for the integration (x in Feigin's notation)

    const double besarg1 = q*radius*sin_alpha;
    //const double besarg2 = q*radius*sin_alpha;

    const double sinarg1 = q*halfheight*cos_alpha;
    const double sinarg2 = q*(halfheight+thick_layer)*cos_alpha;

    const double be1 = sas_J1c(besarg1);
    //const double be2 = sas_J1c(besarg2);
    const double be2 = be1;
    const double si1 = sinc(sinarg1);
    const double si2 = sinc(sinarg2);

    const double dr1 = core_sld - solvent_sld;
    const double dr2 = layer_sld - solvent_sld;
    const double area = M_PI*radius*radius;
    const double totald = 2.0*(thick_layer + halfheight);

    const double t1 = area * (2.0*halfheight) * dr1 * si1 * be1;
    const double t2 = area * dr2 * (totald*si2 - 2.0*halfheight*si1) * be2;

    double pq = square(t1 + t2);

    // loop for the structure factor S(q)
    double qd_cos_alpha = q*d*cos_alpha;
    double debye_arg = -0.5*square(qd_cos_alpha*sigma_dnn);
    double sq=0.0;
    for (int kk=1; kk<n_stacking; kk++) {
        sq += (n_stacking-kk) * cos(qd_cos_alpha*kk) * exp(debye_arg*kk);
    }
    // end of loop for S(q)
    sq = 1.0 + 2.0*sq/n_stacking;

    return pq * sq;
}


static
double stacked_disks_kernel(double q,
                            double thick_core,
                            double thick_layer,
                            double radius,
                            double n_stacking,
                            double sigma_dnn,
                            double core_sld,
                            double layer_sld,
                            double solvent_sld)
{
/*    StackedDiscsX  :  calculates the form factor of a stacked "tactoid" of core shell disks
like clay platelets that are not exfoliated
*/
    double summ = 0.0;    //initialize integral

    double d = 2.0*thick_layer+thick_core;
    double halfheight = 0.5*thick_core;

    for(int i=0; i<N_POINTS_76; i++) {
        double zi = (Gauss76Z[i] + 1.0)*M_PI_4;
        double sin_alpha, cos_alpha; // slots to hold sincos function output
        SINCOS(zi, sin_alpha, cos_alpha);
        double yyy = _kernel(q,
                           radius,
                           core_sld,
                           layer_sld,
                           solvent_sld,
                           halfheight,
                           thick_layer,
                           sin_alpha,
                           cos_alpha,
                           sigma_dnn,
                           d,
                           n_stacking);
        summ += Gauss76Wt[i] * yyy * sin_alpha;
    }

    double answer = M_PI_4*summ;

    //Convert to [cm-1]
    return 1.0e-4*answer;
}

double form_volume(double thick_core,
                   double thick_layer,
                   double radius,
                   double n_stacking){
    double d = 2.0 * thick_layer + thick_core;
    return M_PI * radius * radius * d * n_stacking;
}

double Iq(double q,
          double thick_core,
          double thick_layer,
          double radius,
          double n_stacking,
          double sigma_dnn,
          double core_sld,
          double layer_sld,
          double solvent_sld)
{
    return stacked_disks_kernel(q,
                    thick_core,
                    thick_layer,
                    radius,
                    n_stacking,
                    sigma_dnn,
                    core_sld,
                    layer_sld,
                    solvent_sld);
}


double
Iqxy(double qx, double qy,
     double thick_core,
     double thick_layer,
     double radius,
     double n_stacking,
     double sigma_dnn,
     double core_sld,
     double layer_sld,
     double solvent_sld,
     double theta,
     double phi)
{
    double q, sin_alpha, cos_alpha;
    ORIENT_SYMMETRIC(qx, qy, theta, phi, q, sin_alpha, cos_alpha);

    double d = 2.0 * thick_layer + thick_core;
    double halfheight = 0.5*thick_core;
    double answer = _kernel(q,
                     radius,
                     core_sld,
                     layer_sld,
                     solvent_sld,
                     halfheight,
                     thick_layer,
                     sin_alpha,
                     cos_alpha,
                     sigma_dnn,
                     d,
                     n_stacking);

    //convert to [cm-1]
    answer *= 1.0e-4;

    return answer;
}

