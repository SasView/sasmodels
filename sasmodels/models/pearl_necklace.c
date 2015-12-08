double _pearl_necklace_kernel(double q, double radius, double edge_separation,
	double thick_string, double num_pearls, double sld_pearl,
	double sld_string, double sld_solv);
double form_volume(double radius, double edge_separation,
	double string_thickness, double number_of_pearls);
double sinc(double x);
	
double Iq(double q, double radius, double edge_separation,
	double string_thickness, double number_of_pearls, double sld, 
	double string_sld, double solvent_sld);
double Iqxy(double qx, double qy, double radius, double edge_separation,
	double string_thickness, double number_of_pearls, double sld, 
	double string_sld, double solvent_sld);
	
double ER(double radius, double edge_separation,
	double string_thickness, double number_of_pearls);
double VR(double radius, double edge_separation,
	double string_thickness, double number_of_pearls);

// From Igor library
double _pearl_necklace_kernel(double q, double radius, double edge_separation, double thick_string,
	double num_pearls, double sld_pearl, double sld_string, double sld_solv)
{
	double contrast_pearl = sld_pearl - sld_solv;
	double contrast_string = sld_string - sld_solv;

	//total volume
	double pi = 4.0*atan(1.0);
	double tot_vol = form_volume(radius, edge_separation, thick_string, num_pearls);
	double string_vol = edge_separation * pi * pow((thick_string / 2.0), 2);
	double pearl_vol = 4.0 /3.0 * pi * pow(radius, 3);
	double num_strings = num_pearls - 1;
	//each mass
	double m_r= contrast_string * string_vol;
	double m_s= contrast_pearl * pearl_vol;
	double psi, gamma, beta;
	//form factors
	double sss; //pearls
	double srr; //strings
	double srs; //cross
	double A_s, srr_1, srr_2, srr_3;
	double form_factor;

	//sine functions of a pearl
	psi = sin(q*radius);
	psi -= q * radius * cos(q * radius);
	psi /= pow((q * radius), 3);
	psi *= 3.0;

	// Note take only 20 terms in Si series: 10 terms may be enough though.
	gamma = Si(q* edge_separation);
	gamma /= (q* edge_separation);
	beta = Si(q * (edge_separation + radius));
	beta -= Si(q * radius);
	beta /= (q* edge_separation);

	// center to center distance between the neighboring pearls
	A_s = edge_separation + 2.0 * radius;

	// form factor for num_pearls
	sss = 1.0 - pow(sinc(q*A_s), num_pearls );
	sss /= pow((1.0-sinc(q*A_s)), 2);
	sss *= -sinc(q*A_s);
	sss -= num_pearls/2.0;
	sss += num_pearls/(1.0-sinc(q*A_s));
	sss *= 2.0 * pow((m_s*psi), 2);

	// form factor for num_strings (like thin rods)
	srr_1 = -pow(sinc(q*edge_separation/2.0), 2);

	srr_1 += 2.0 * gamma;
	srr_1 *= num_strings;
	srr_2 = 2.0/(1.0-sinc(q*A_s));
	srr_2 *= num_strings * pow(beta, 2);
	srr_3 = 1.0 - pow(sinc(q*A_s), num_strings);
	srr_3 /= pow((1.0-sinc(q*A_s)), 2);
	srr_3 *= -2.0 * pow(beta, 2);

	// total srr
	srr = srr_1 + srr_2 + srr_3;
	srr *= pow(m_r, 2);

	// form factor for correlations
	srs = 1.0;
	srs -= pow(sinc(q*A_s), num_strings);
	srs /= pow((1.0-sinc(q*A_s)), 2);
	srs *= -sinc(q*A_s);
	srs += (num_strings/(1.0-sinc(q*A_s)));
	srs *= 4.0 * (m_r * m_s * beta * psi);

	form_factor = sss + srr + srs;
	form_factor /= (tot_vol * 1.0e4); // norm by volume and A^-1 to cm^-1

	return (form_factor);
}

double form_volume(double radius, double edge_separation,
	double string_thickness, double number_of_pearls)
{
	double total_vol;

	double pi = 4.0*atan(1.0);
	double number_of_strings = number_of_pearls - 1.0;
	
	double string_vol = edge_separation * pi * pow((string_thickness / 2.0), 2);
	double pearl_vol = 4.0 /3.0 * pi * pow(radius, 3);

	total_vol = number_of_strings * string_vol;
	total_vol += number_of_pearls * pearl_vol;

	return(total_vol);
}

double sinc(double x)
{
	double num = sin(x);
	double denom = x;
	return num/denom;
}


double Iq(double q, double radius, double edge_separation,
	double string_thickness, double number_of_pearls, double sld, 
	double string_sld, double solvent_sld)
{
	double value = 0.0;
	double tot_vol = 0.0;
	
	value = _pearl_necklace_kernel(q, radius, edge_separation, string_thickness,
		number_of_pearls, sld, string_sld, solvent_sld);
	tot_vol = form_volume(radius, edge_separation, string_thickness, number_of_pearls);

	return value*tot_vol;
}

double Iqxy(double qx, double qy, double radius, double edge_separation,
	double string_thickness, double number_of_pearls, double sld, 
	double string_sld, double solvent_sld)
{
    double q = sqrt(qx*qx + qy*qy);
    return(Iq(q, radius, edge_separation, string_thickness, number_of_pearls, 
		sld, string_sld, solvent_sld));
}


double ER(double radius, double edge_separation,
	double string_thickness, double number_of_pearls)
{
	double tot_vol = form_volume(radius, edge_separation, string_thickness, number_of_pearls);
	double pi = 4.0*atan(1.0);

    double rad_out = pow((3.0*tot_vol/4.0/pi), 0.33333);
    
    return rad_out;
}
double VR(double radius, double edge_separation,
	double string_thickness, double number_of_pearls)
{
	return 1.0;
}