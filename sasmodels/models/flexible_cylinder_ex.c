double form_volume(double length, double kuhn_length, double radius);
double Iq(double q, double length, double kuhn_length, double radius,
          double axis_ratio, double sld, double solvent_sld);
double Iqxy(double qx, double qy, double length, double kuhn_length,
            double radius, double axis_ratio, double sld, double solvent_sld);
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
    double uplim,lolim,Pi,summ,arg,zi,yyy,answer;
    int i,nord=76;

    Pi = 4.0*atan(1.0);
    lolim=0.0;
    uplim=Pi/2.0;
    summ=0.0;

    for(i=0;i<nord;i++) {
		zi = ( Gauss76Z[i]*(uplim-lolim) + uplim + lolim )/2.0;
		arg = q*sqrt(a*a*sin(zi)*sin(zi)+b*b*cos(zi)*cos(zi));
		yyy = pow((2.0 * J1(arg) / arg),2);
		yyy *= Gauss76Wt[i];
		summ += yyy;
    }
    answer = (uplim-lolim)/2.0*summ;
    answer *= 2.0/Pi;
    return(answer);

}

double flexible_cylinder_ex_kernel(double q,
          double length,
          double kuhn_length,
          double radius,
          double axis_ratio,
          double sld,
          double solvent_sld)
{

	double Pi,flex,crossSect, cont;

	Pi = 4.0*atan(1.0);
	cont = sld - solvent_sld;
	crossSect = elliptical_crosssection(q,radius,(radius*axis_ratio));

	flex = Sk_WR(q,length,kuhn_length);
	flex *= crossSect;
	flex *= Pi*radius*radius*axis_ratio*axis_ratio*length;
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

double Iqxy(double qx, double qy,
            double length,
            double kuhn_length,
            double radius,
            double axis_ratio,
            double sld,
            double solvent_sld)
{
	double q;
	q = sqrt(qx*qx+qy*qy);
	double result = flexible_cylinder_ex_kernel(q,
	                length,
	                kuhn_length,
	                radius,
	                axis_ratio,
	                sld,
	                solvent_sld);

	return result;
}
