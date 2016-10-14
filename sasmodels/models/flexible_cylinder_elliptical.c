double form_volume(double length, double kuhn_length, double radius);
double Iq(double q, double length, double kuhn_length, double radius,
          double axis_ratio, double sld, double solvent_sld);
double flexible_cylinder_ex_kernel(double q, double length, double kuhn_length,
                                double radius, double axis_ratio, double sld,
                                double solvent_sld);
double elliptical_crosssection(double q, double a, double b);

double form_volume(double length, double kuhn_length, double radius)
{
    return 1.0;
}

double
elliptical_crosssection(double q, double a, double b)
{
    double summ=0.0;

    for(int i=0;i<N_POINTS_76;i++) {
        double zi = ( Gauss76Z[i] + 1.0 )*M_PI/4.0;
        double sn, cn;
        SINCOS(zi, sn, cn);
        double arg = q*sqrt(a*a*sn*sn+b*b*cn*cn);
        double yyy = pow((double)sas_J1c(arg),2);
        yyy *= Gauss76Wt[i];
        summ += yyy;
    }

    summ /= 2.0;
    return(summ);

}

double flexible_cylinder_ex_kernel(double q,
          double length,
          double kuhn_length,
          double radius,
          double axis_ratio,
          double sld,
          double solvent_sld)
{

    double flex,crossSect, cont;

    cont = sld - solvent_sld;
    crossSect = elliptical_crosssection(q,radius,(radius*axis_ratio));

    flex = Sk_WR(q,length,kuhn_length);
    flex *= crossSect;
    flex *= M_PI*radius*radius*axis_ratio*axis_ratio*length;
    flex *= cont*cont;
    flex *= 1.0e-4;

    return flex;
}

double Iq(double q,
          double length,
          double kuhn_length,
          double radius,
          double axis_ratio,
          double sld,
          double solvent_sld)
{

    double result = flexible_cylinder_ex_kernel(q,
                    length,
                    kuhn_length,
                    radius,
                    axis_ratio,
                    sld,
                    solvent_sld);

    return result;
}

