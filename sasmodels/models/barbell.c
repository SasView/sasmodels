double form_volume(double bell_radius, double radius, double length);
double Iq(double q, double sld, double solvent_sld, double bell_radius, double radius, double length);
double Iqxy(double qx, double qy, double sld, double solvent_sld,
    double bell_radius, double radius, double length, double theta, double phi);

//barbell kernel - same as dumbell
double _bell_kernel(double q, double h, double bell_radius,
    double length, double sin_alpha, double cos_alpha);
double _bell_kernel(double q, double h, double bell_radius,
    double length, double sin_alpha, double cos_alpha)
{
    const double upper = 1.0;
    const double lower = -1.0*h/bell_radius;

    double total = 0.0;
    for (int i = 0; i < 76; i++){
        const double t = 0.5*(Gauss76Z[i]*(upper-lower)+upper+lower);
	    const double arg1 = q*cos_alpha*(bell_radius*t+h+length*0.5);
	    const double arg2 = q*bell_radius*sin_alpha*sqrt(1.0-t*t);

        const double be = (arg2 == 0.0 ? 0.5 :J1(arg2)/arg2);

	    const double Fq = cos(arg1)*(1.0-t*t)*be;

	    total += Gauss76Wt[i] * Fq;
    }
    const double integral = 0.5*(upper-lower)*total;
    return 4.0*M_PI*bell_radius*bell_radius*bell_radius*integral;
}

double form_volume(double bell_radius,
        double radius,
        double length)
{

    // bell radius should never be less than radius when this is called
    const double hdist = sqrt(bell_radius*bell_radius - radius*radius);
    const double p1 = 2.0*bell_radius*bell_radius*bell_radius/3.0;
    const double p2 = bell_radius*bell_radius*hdist;
    const double p3 = hdist*hdist*hdist/3.0;

    return M_PI*radius*radius*length + 2.0*M_PI*(p1+p2-p3);
}

double Iq(double q, double sld,
    double solvent_sld,
    double bell_radius,
    double radius,
    double length)
{
    double sn, cn; // slots to hold sincos function output

    if (bell_radius < radius) return -1.0;

    const double lower = 0.0;
    const double upper = M_PI_2;
    const double h = sqrt(bell_radius*bell_radius-radius*radius);
    double total = 0.0;
    for (int i = 0; i < 76; i++){
        const double alpha= 0.5*(Gauss76Z[i]*(upper-lower) + upper + lower);
        SINCOS(alpha, sn, cn);

        const double bell_Fq = _bell_kernel(q, h, bell_radius, length, sn, cn);

        const double arg1 = q*length*0.5*cn;
        const double arg2 = q*radius*sn;
        // lim_{x->0} J1(x)/x = 1/2,   lim_{x->0} sin(x)/x = 1
        const double be = (arg2 == 0.0 ? 0.5 :J1(arg2)/arg2);
        const double si = (arg1 == 0.0 ? 1.0 :sin(arg1)/arg1);
        const double cyl_Fq = M_PI*radius*radius*length*si*2.0*be;

        const double Aq = cyl_Fq + bell_Fq;
        total += Gauss76Wt[i] * Aq * Aq * sn;
    }

    const double form = total*(upper-lower)*0.5;

    //Contrast and volume normalization
    const double s = (sld - solvent_sld);
    return form*1.0e-4*s*s; //form_volume(bell_radius,radius,length);
}



double Iqxy(double qx, double qy,
    double sld,
    double solvent_sld,
    double bell_radius,
    double radius,
    double length,
    double theta,
    double phi)
{
     double sn, cn; // slots to hold sincos function output

    // Exclude invalid inputs.
    if (bell_radius < radius) return -1.0;

    // Compute angle alpha between q and the cylinder axis
    SINCOS(theta*M_PI_180, sn, cn);
    // # The following correction factor exists in sasview, but it can't be
    // # right, so we are leaving it out for now.
    const double q = sqrt(qx*qx+qy*qy);
    const double cos_val = cn*cos(phi*M_PI_180)*qx + sn*qy;
    const double alpha = acos(cos_val); // rod angle relative to q
    SINCOS(alpha, sn, cn);

    const double h = sqrt(bell_radius*bell_radius - radius*radius); // negative h
    const double bell_Fq = _bell_kernel(q, h, bell_radius, length, sn, cn)/sn;

    const double besarg = q*radius*sn;
    const double siarg = q*0.5*length*cn;
    // lim_{x->0} J1(x)/x = 1/2,   lim_{x->0} sin(x)/x = 1
    const double bj = (besarg == 0.0 ? 0.5 : J1(besarg)/besarg);
    const double si = (siarg == 0.0 ? 1.0 : sin(siarg)/siarg);
    const double cyl_Fq = M_PI*radius*radius*length*2.0*bj*si;

    // Volume weighted average F(q)
    const double Aq = cyl_Fq + bell_Fq;

    // Multiply by contrast^2, normalize by cylinder volume and convert to cm-1
    const double s = (sld - solvent_sld);
    return 1.0e-4 * Aq * Aq * s * s; // form_volume(radius, cap_radius, length);
}
