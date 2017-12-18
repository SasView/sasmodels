//#define INVALID(v) (v.radius_core >= v.radius)

// From Igor library
static double
_hollow_cylinder_scaling(double integrand, double delrho, double volume)
{
    return 1.0e-4 * square(volume * delrho) * integrand;
}

static double
_fq(double qab, double qc,
    double radius, double thickness, double length)
{
    const double lam1 = sas_2J1x_x((radius+thickness)*qab);
    const double lam2 = sas_2J1x_x(radius*qab);
    const double gamma_sq = square(radius/(radius+thickness));
    //Note: lim_{thickness -> 0} psi = sas_J0(radius*qab)
    //Note: lim_{radius -> 0} psi = sas_2J1x_x(thickness*qab)
    const double psi = (lam1 - gamma_sq*lam2)/(1.0 - gamma_sq);    //SRK 10/19/00
    const double t2 = sas_sinx_x(0.5*length*qc);
    return psi*t2;
}

static double
form_volume(double radius, double thickness, double length)
{
    double v_shell = M_PI*length*(square(radius+thickness) - radius*radius);
    return v_shell;
}


static double
Iq(double q, double radius, double thickness, double length,
    double sld, double solvent_sld)
{
    const double lower = 0.0;
    const double upper = 1.0;        //limits of numerical integral

    double summ = 0.0;            //initialize intergral
    for (int i=0;i<GAUSS_N;i++) {
        const double cos_theta = 0.5*( GAUSS_Z[i] * (upper-lower) + lower + upper );
        const double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
        const double form = _fq(q*sin_theta, q*cos_theta,
                                radius, thickness, length);
        summ += GAUSS_W[i] * form * form;
    }

    const double Aq = 0.5*summ*(upper-lower);
    const double volume = form_volume(radius, thickness, length);
    return _hollow_cylinder_scaling(Aq, solvent_sld - sld, volume);
}

static double
Iqac(double qab, double qc,
    double radius, double thickness, double length,
    double sld, double solvent_sld)
{
    const double form = _fq(qab, qc, radius, thickness, length);

    const double vol = form_volume(radius, thickness, length);
    return _hollow_cylinder_scaling(form*form, solvent_sld-sld, vol);
}
