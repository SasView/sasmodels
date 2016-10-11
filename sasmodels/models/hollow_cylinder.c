double form_volume(double radius, double thickness, double length);

double Iq(double q, double radius, double thickness, double length, double sld,
	double solvent_sld);
double Iqxy(double qx, double qy, double radius, double thickness, double length, double sld,
	double solvent_sld, double theta, double phi);

//#define INVALID(v) (v.radius_core >= v.radius)

// From Igor library
static double hollow_cylinder_scaling(
    double integrand, double delrho, double volume)
{
    double answer;
    // Multiply by contrast^2
    answer = integrand*delrho*delrho;

    //normalize by cylinder volume
    answer *= volume*volume;

    //convert to [cm-1]
    answer *= 1.0e-4;

    return answer;
}


static double _hollow_cylinder_kernel(
    double q, double radius, double thickness, double length, double dum)
{
    const double qs = q*sqrt(1.0-dum*dum);
    const double lam1 = sas_J1c((radius+thickness)*qs);
    const double lam2 = sas_J1c(radius*qs);
    const double gamma_sq = square(radius/(radius+thickness));
    //Note: lim_{r -> r_c} psi = J0(radius_core*qs)
    const double psi = (lam1 - gamma_sq*lam2)/(1.0 - gamma_sq);	//SRK 10/19/00
    const double t2 = sinc(q*length*dum/2.0);
    return square(psi*t2);
}


static double hollow_cylinder_analytical_2D_scaled(
    double q, double q_x, double q_y, double radius, double thickness,
    double length, double sld, double solvent_sld, double theta, double phi)
{
    double cyl_x, cyl_y; //, cyl_z
    //double q_z;
    double vol, cos_val, delrho;
    double answer;
    //convert angle degree to radian
    theta = theta * M_PI_180;
    phi = phi * M_PI_180;
    delrho = solvent_sld - sld;

    // Cylinder orientation
    cyl_x = sin(theta) * cos(phi);
    cyl_y = sin(theta) * sin(phi);
    //cyl_z = -cos(theta) * sin(phi);

    // q vector
    //q_z = 0;

    // Compute the angle btw vector q and the
    // axis of the cylinder
    cos_val = cyl_x*q_x + cyl_y*q_y;// + cyl_z*q_z;

    answer = _hollow_cylinder_kernel(q, radius, thickness, length, cos_val);

    vol = form_volume(radius, thickness, length);
    answer = hollow_cylinder_scaling(answer, delrho, vol);

    return answer;
}


double form_volume(double radius, double thickness, double length)
{
    double v_shell = M_PI*length*((radius+thickness)*(radius+thickness)-radius*radius);
    return(v_shell);
}


double Iq(double q, double radius, double thickness, double length,
    double sld, double solvent_sld)
{
    int i;
    double lower,upper,zi, inter;		//upper and lower integration limits
    double summ,answer,delrho;			//running tally of integration
    double norm,volume;	//final calculation variables

    lower = 0.0;
    upper = 1.0;		//limits of numerical integral

    summ = 0.0;			//initialize intergral
    for (i=0;i<76;i++) {
        zi = ( Gauss76Z[i] * (upper-lower) + lower + upper )/2.0;
        inter = Gauss76Wt[i] * _hollow_cylinder_kernel(q, radius, thickness, length, zi);
        summ += inter;
    }

    norm = summ*(upper-lower)/2.0;
    volume = form_volume(radius, thickness, length);
    delrho = solvent_sld - sld;
    answer = hollow_cylinder_scaling(norm, delrho, volume);

    return(answer);
}


double Iqxy(double qx, double qy, double radius, double thickness,
    double length, double sld, double solvent_sld, double theta, double phi)
{
    const double q = sqrt(qx*qx+qy*qy);
    return hollow_cylinder_analytical_2D_scaled(q, qx/q, qy/q, radius, thickness, length, sld, solvent_sld, theta, phi);
}