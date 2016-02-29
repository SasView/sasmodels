double form_volume(double length, double kuhn_length, double radius);
double Iq(double q, double length, double kuhn_length, double radius,
          double sld, double solvent_sld);
double Iqxy(double qx, double qy, double length, double kuhn_length,
            double radius, double sld, double solvent_sld);
double flexible_cylinder_kernel(double q, double length, double kuhn_length,
                                double radius, double sld, double solvent_sld);


double form_volume(double length, double kuhn_length, double radius)
{

      return 0.0;
}

double flexible_cylinder_kernel(double q,
          double length,
          double kuhn_length,
          double radius,
          double sld,
          double solvent_sld)
{

    const double cont = sld-solvent_sld;
    const double qr = q*radius;
    //const double crossSect = (2.0*J1(qr)/qr)*(2.0*J1(qr)/qr);
    const double crossSect = J1c(qr);
    double flex = Sk_WR(q,length,kuhn_length);
    flex *= crossSect*crossSect;
    flex *= M_PI*radius*radius*length;
    flex *= cont*cont;
    flex *= 1.0e-4;
    return flex;
}

double Iq(double q,
          double length,
          double kuhn_length,
          double radius,
          double sld,
          double solvent_sld)
{

    double result = flexible_cylinder_kernel(q, length, kuhn_length, radius, sld, solvent_sld);
    return result;
}

double Iqxy(double qx, double qy,
            double length,
            double kuhn_length,
            double radius,
            double sld,
            double solvent_sld)
{
    double q;
    q = sqrt(qx*qx+qy*qy);
    double result = flexible_cylinder_kernel(q, length, kuhn_length, radius, sld, solvent_sld);

    return result;
}
