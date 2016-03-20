double form_volume(double radius, double thickness);

double Iq(double q,
          double radius,
          double thickness,
          double core_sld,
          double shell_sld,
          double solvent_sld,
          double volfraction,
          double frac_dim,
          double cor_length);

double Iqxy(double qx, double qy,
            double radius,
            double thickness,
            double core_sld,
            double shell_sld,
            double solvent_sld,
            double volfraction,
            double frac_dim,
            double cor_length);

double form_volume(double radius, double thickness)
{
    return 4.0 * M_PI / 3.0 * pow((radius + thickness), 3);
}

double Iq(double q,
          double radius,
          double thickness,
          double core_sld,
          double shell_sld,
          double solvent_sld,
          double volfraction,
          double frac_dim,
          double cor_length) {


   double intensity = core_shell_kernel(q,
                              radius,
                              thickness,
                              core_sld,
                              shell_sld,
                              solvent_sld);
    //calculate S(q)
    double frac_1 = frac_dim-1.0;
    double qr = q*radius;

    double t1 = frac_dim*sas_gamma(frac_1)*sin(frac_1*atan(q*cor_length));
    double t2 = (1.0 + 1.0/(q*cor_length)/(q*cor_length));
    double t3 = pow(qr, frac_dim) * pow(t2, (frac_1/2.0));
    double sq = t1/t3;
    sq += 1.0;

    return sq*intensity*volfraction;
}

double Iqxy(double qx, double qy,
            double radius,
            double thickness,
            double core_sld,
            double shell_sld,
            double solvent_sld,
            double volfraction,
            double frac_dim,
            double cor_length) {

    const double q = sqrt(qx*qx+qy*qy);
    return Iq(q,
              radius,
              thickness,
              core_sld,
              shell_sld,
              solvent_sld,
              volfraction,
              frac_dim,
              cor_length);

}


