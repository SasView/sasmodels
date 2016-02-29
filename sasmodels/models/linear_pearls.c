double form_volume(double radius, double num_pearls);

double Iq(double q,
            double radius,
            double edge_sep,
            double num_pearls,
            double pearl_sld,
            double solvent_sld);

double Iqxy(double qx, double qy,
            double radius,
            double edge_sep,
            double num_pearls,
            double pearl_sld,
            double solvent_sld);

double linear_pearls_kernel(double q,
            double radius,
            double edge_sep,
            double num_pearls,
            double pearl_sld,
            double solvent_sld);


double form_volume(double radius, double num_pearls)
{
    // Pearl volume
    double pearl_vol = 4.0 /3.0 * M_PI * pow(radius, 3.0);
    // Return total volume
    return num_pearls * pearl_vol;;
}

double linear_pearls_kernel(double q,
            double radius,
            double edge_sep,
            double num_pearls,
            double pearl_sld,
            double solvent_sld)
{
    double n_contrib;
    //relative sld
    double contrast_pearl = pearl_sld - solvent_sld;
    //each volume
    double pearl_vol = 4.0 /3.0 * M_PI * pow(radius, 3.0);
    //total volume
    double tot_vol = num_pearls * pearl_vol;
    //mass
    double m_s = contrast_pearl * pearl_vol;
    //center to center distance between the neighboring pearls
    double separation = edge_sep + 2.0 * radius;

    double x=q*radius;

    // Try Taylor on x*xos(x)
	// double out_cos = x - pow(x,3)/2 + pow(x,5)/24 - pow(x,7)/720 + pow(x,9)/40320;
    // psi -= x*out_cos;

    //sine functions of a pearl
    double psi = sin(q * radius);
    psi -= x * cos(x);
    psi /= pow((q * radius), 3.0);

    // N pearls contribution
    int n_max = num_pearls - 1;
    n_contrib = num_pearls;
    for(int num=1; num<=n_max; num++) {
        n_contrib += (2.0*(num_pearls-num)*sinc(q*separation*num));
    }
    // form factor for num_pearls
    double form_factor = n_contrib;
    form_factor *= pow((m_s*psi*3.0), 2.0);
    form_factor /= (tot_vol * 1.0e4);

    return form_factor;
}

double Iq(double q,
            double radius,
            double edge_sep,
            double num_pearls,
            double pearl_sld,
            double solvent_sld)
{

	double result = linear_pearls_kernel(q,
                    radius,
                    edge_sep,
                    num_pearls,
                    pearl_sld,
                    solvent_sld);

	return result;
}

double Iqxy(double qx, double qy,
            double radius,
            double edge_sep,
            double num_pearls,
            double pearl_sld,
            double solvent_sld)
{
	double q;
	q = sqrt(qx*qx+qy*qy);

	double result = linear_pearls_kernel(q,
                    radius,
                    edge_sep,
                    num_pearls,
                    pearl_sld,
                    solvent_sld);

	return result;
}
