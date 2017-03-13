static double
form_volume(double radius,
          double thick_shell,
          double thick_solvent,
          double fp_n_shells)
{
    int n_shells = (int)(fp_n_shells + 0.5);
    double R_N = radius + n_shells*(thick_shell+thick_solvent) - thick_solvent;
    return M_4PI_3*cube(R_N);
}

static double
multilayer_vesicle_kernel(double q,
          double volfraction,
          double radius,
          double thick_shell,
          double thick_solvent,
          double sld_solvent,
          double sld,
          int n_shells)
{
    //calculate with a loop, two shells at a time
    int ii = 0;
    double fval = 0.0;
    double voli = 0.0;
    const double sldi = sld_solvent-sld;

    do {
        double ri = radius + (double)ii*(thick_shell + thick_solvent);

        // layer 1
        voli = M_4PI_3*ri*ri*ri;
        fval += voli*sldi*sas_3j1x_x(ri*q);

        ri += thick_shell;

        // layer 2
        voli = M_4PI_3*ri*ri*ri;
        fval -= voli*sldi*sas_3j1x_x(ri*q);

        //do 2 layers at a time
        ii++;

    } while(ii <= n_shells-1);  //change to make 0 < n_shells < 2 correspond to
                               //unilamellar vesicles (C. Glinka, 11/24/03)

    return 1.0e-4*volfraction*fval*fval;  // Volume normalization happens in caller
}

static double
Iq(double q,
          double volfraction,
          double radius,
          double thick_shell,
          double thick_solvent,
          double sld_solvent,
          double sld,
          double fp_n_shells)
{
    int n_shells = (int)(fp_n_shells + 0.5);
    return multilayer_vesicle_kernel(q,
           volfraction,
           radius,
           thick_shell,
           thick_solvent,
           sld_solvent,
           sld,
           n_shells);
}

