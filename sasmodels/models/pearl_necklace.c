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

// From Igor library
double _pearl_necklace_kernel(double q, double radius, double edge_separation, double thick_string,
	double num_pearls, double sld_pearl, double sld_string, double sld_solv)
{
	//relative slds
	double contrast_pearl = sld_pearl - sld_solv;
	double contrast_string = sld_string - sld_solv;
	
	// number of string segments
	num_pearls = floor(num_pearls + 0.5); //Force integer number of pearls
	double num_strings = num_pearls - 1.0;
	
	//Pi
	double pi = 4.0*atan(1.0);
	
	// center to center distance between the neighboring pearls
	double A_s = edge_separation + 2.0 * radius;
	
	// Repeated Calculations
	double sincasq = sinc(q*A_s);
	double oneminussinc = 1 - sincasq;
	double q_r = q * radius;
	double q_edge = q * edge_separation;
	
	// each volume
	double string_vol = edge_separation * pi * thick_string * thick_string / 4.0;
	double pearl_vol = 4.0 / 3.0 * pi * radius * radius * radius;

	//total volume
	double tot_vol;
	//each masses
	double m_r= contrast_string * string_vol;
	double m_s= contrast_pearl * pearl_vol;
	double psi, gamma, beta;
	//form factors
	double sss, srr, srs; //cross
	double srr_1, srr_2, srr_3;
	double form_factor;
	tot_vol = num_strings * string_vol;
	tot_vol += num_pearls * pearl_vol;

	//sine functions of a pearl
	psi = sin(q_r);
	psi -= q_r * cos(q_r);
	psi *= 3.0;
	psi /= q_r * q_r * q_r;

	// Note take only 20 terms in Si series: 10 terms may be enough though.
	gamma = Si(q_edge);
	gamma /= (q_edge);
	beta = Si(q * (A_s - radius));
	beta -= Si(q_r);
	beta /= q_edge;

	// form factor for num_pearls
	sss = 1.0 - pow(sincasq, num_pearls);
	sss /= oneminussinc * oneminussinc;
	sss *= -sincasq;
	sss -= num_pearls / 2.0;
	sss += num_pearls / oneminussinc;
	sss *= 2.0 * m_s * psi * m_s * psi;

	// form factor for num_strings (like thin rods)
	srr_1 = -sinc(q_edge/2.0) * sinc(q_edge/2.0);

	srr_1 += 2.0 * gamma;
	srr_1 *= num_strings;
	srr_2 = 2.0/oneminussinc;
	srr_2 *= num_strings;
	srr_2 *= beta * beta;
	srr_3 = 1.0 - pow(sincasq, num_strings);
	srr_3 /= oneminussinc * oneminussinc;
	srr_3 *= beta * beta;
	srr_3 *= -2.0;

	// total srr
	srr = srr_1 + srr_2 + srr_3;
	srr *= m_r * m_r;

	// form factor for correlations
	srs = 1.0;
	srs -= pow(sincasq, num_strings);
	srs /= oneminussinc * oneminussinc;
	srs *= -sincasq;
	srs += num_strings/oneminussinc;
	srs *= 4.0;
	srs *= (m_r * m_s * beta * psi);

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
	
	double string_vol = edge_separation * pi * string_thickness * string_thickness / 4.0;
	double pearl_vol = 4.0 / 3.0 * pi * radius * radius * radius;

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
	double value, tot_vol;
	
	if (string_thickness >= radius || number_of_pearls <= 0) {
		return NAN;
	}
	
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