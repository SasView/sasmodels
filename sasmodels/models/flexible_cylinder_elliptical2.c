#include <stdio.h>

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
integrand_elliptical_crosssection(double x, double q, double a,  double b, double* out, int n, int i){
    double sn, cn;
    get_sin_x(n, i, &sn);
    get_cos_x(n, i, &cn);

    const double arg = q * sqrt(a*a*sn*sn + b*b*cn*cn);
    const double yyy = sas_2J1x_x(arg);
    *out = yyy * yyy;
}

double
elliptical_crosssection(double q, double a, double b)
{
    double sum;
    integrate(integrand_elliptical_crosssection, 0, M_PI_2, q, a, b, &sum);
    return sum;

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
