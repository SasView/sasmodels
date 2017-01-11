double form_volume(double radius, double thickness, double length);
double Iq(double q, double radius, double thickness, double length, double sld,
	double solvent_sld);
double Iqxy(double qx, double qy, double radius, double thickness, double length, double sld,
	double solvent_sld, double theta, double phi);

//#define INVALID(v) (v.radius_core >= v.radius)

// From Igor library
static double
_hollow_cylinder_scaling(double integrand, double delrho, double volume)
{
    return 1.0e-4 * square(volume * delrho) * integrand;
}


static double
_hollow_cylinder_kernel(double q,
    double radius, double thickness, double length, double sin_val, double cos_val)
{
    const double qs = q*sin_val;
    const double lam1 = sas_J1c((radius+thickness)*qs);
    const double lam2 = sas_J1c(radius*qs);
    const double gamma_sq = square(radius/(radius+thickness));
    //Note: lim_{thickness -> 0} psi = J0(radius*qs)
    //Note: lim_{radius -> 0} psi = sas_J1c(thickness*qs)
    const double psi = (lam1 - gamma_sq*lam2)/(1.0 - gamma_sq);	//SRK 10/19/00
    const double t2 = sas_sinx_x(0.5*q*length*cos_val);
    return psi*t2;
}

double
form_volume(double radius, double thickness, double length)
{
    double v_shell = M_PI*length*(square(radius+thickness) - radius*radius);
    return v_shell;
}


double
Iq(double q, double radius, double thickness, double length,
    double sld, double solvent_sld)
{
    const double lower = 0.0;
    const double upper = 1.0;		//limits of numerical integral

    double summ = 0.0;			//initialize intergral
    for (int i=0;i<76;i++) {
        const double cos_val = 0.5*( Gauss76Z[i] * (upper-lower) + lower + upper );
        const double sin_val = sqrt(1.0 - cos_val*cos_val);
        const double inter = _hollow_cylinder_kernel(q, radius, thickness, length,
                                                     sin_val, cos_val);
        summ += Gauss76Wt[i] * inter * inter;
    }

    const double Aq = 0.5*summ*(upper-lower);
    const double volume = form_volume(radius, thickness, length);
    return _hollow_cylinder_scaling(Aq, solvent_sld - sld, volume);
}

double
Iqxy(double qx, double qy,
    double radius, double thickness, double length,
    double sld, double solvent_sld, double theta, double phi)
{
    double q, sin_alpha, cos_alpha;
    ORIENT_SYMMETRIC(qx, qy, theta, phi, q, sin_alpha, cos_alpha);
    const double Aq = _hollow_cylinder_kernel(q, radius, thickness, length,
        sin_alpha, cos_alpha);

    const double vol = form_volume(radius, thickness, length);
    return _hollow_cylinder_scaling(Aq*Aq, solvent_sld-sld, vol);
}

