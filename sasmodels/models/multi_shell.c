double Iq(double q,
          double core_radius,
          double s_thickness,
          double w_thickness,
          double core_sld,
          double shell_sld,
          double n_pairs);

double Iqxy(double qx, double qy,
          double core_radius,
          double s_thickness,
          double w_thickness,
          double core_sld,
          double shell_sld,
          double n_pairs);

static
double multi_shell_kernel(double q,
          double core_radius,
          double s_thickness,
          double w_thickness,
          double core_sld,
          double shell_sld,
          double n_pairs)
{
	//calculate with a loop, two shells at a time
    int ii=0;
    double fval=0.0;
    double voli = 0.0;
    const double sldi = core_sld-shell_sld;

    do {
        double ri = core_radius + (double)ii*(s_thickness + w_thickness);

        // layer 1
        voli = 4.0*M_PI/3.0*ri*ri*ri;
        fval += voli*sldi*sph_j1c(ri*q);

        ri += s_thickness;

        // layer 2
        voli = 4.0*M_PI/3.0*ri*ri*ri;
        fval -= voli*sldi*sph_j1c(ri*q);

        //do 2 layers at a time
        ii+=1;

	} while(ii<=n_pairs-1);  //change to make 0 < n_pairs < 2 correspond to
	                         //unilamellar vesicles (C. Glinka, 11/24/03)

    fval *= 1.0e-4*fval/voli;

    return(fval);
}

double Iq(double q,
          double core_radius,
          double s_thickness,
          double w_thickness,
          double core_sld,
          double shell_sld,
          double n_pairs)
{
    return multi_shell_kernel(q,
           core_radius,
           s_thickness,
           w_thickness,
           core_sld,
           shell_sld,
           n_pairs);
}

double Iqxy(double qx, double qy,
          double core_radius,
          double s_thickness,
          double w_thickness,
          double core_sld,
          double shell_sld,
          double n_pairs)
{
    double q = sqrt(qx*qx + qy*qy);

    return multi_shell_kernel(q,
           core_radius,
           s_thickness,
           w_thickness,
           core_sld,
           shell_sld,
           n_pairs);
}

