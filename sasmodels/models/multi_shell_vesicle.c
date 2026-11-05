/*
 * Scattering from a multi-shell vesicle: a solvent-filled core surrounded
 * by up to 10 concentric shells of individually specified thickness and SLD.
 * The SLD of the core is equal to that of the solvent, so the scattering
 * volume is defined by the shells alone.
 *
 * For n=1 with a uniform shell SLD this reduces to the vesicle model.
 */

static double
outer_radius(double radius, double fp_n, double thickness[])
{
    double r = radius;
    const int n = (int)(fp_n+0.5);
    for (int i = 0; i < n; i++) {
        r += thickness[i];
    }
    return r;
}

static double
form_volume(double radius, double fp_n, double thickness[])
{
    // Enclosed volume (shell material + solvent core)
    return M_4PI_3 * cube(outer_radius(radius, fp_n, thickness));
}

static double
shell_volume(double radius, double fp_n, double thickness[])
{
    // Volume of the shell material only (used by sasmodels for the
    // hollow-form normalization so that volfraction refers to the
    // shell material rather than the enclosed sphere)
    return M_4PI_3 * (cube(outer_radius(radius, fp_n, thickness)) - cube(radius));
}

static double
radius_effective(int mode, double radius, double fp_n, double thickness[])
{
    switch (mode) {
    default:
    case 1: // outer radius
        return outer_radius(radius, fp_n, thickness);
    case 2: // core radius
        return radius;
    }
}

static void
Fq(double q, double *F1, double *F2,
   double radius, double sld_solvent, double fp_n,
   double sld[], double thickness[], double volfraction)
{
    const int n = (int)(fp_n+0.5);
    double vol, r, contrast, f, last_sld;

    r = radius;
    last_sld = sld_solvent;  // core is solvent, so no contrast at r < radius
    f = 0.;
    for (int i = 0; i < n; i++) {
        // Amplitude contribution at the inner edge of shell i
        // (contrast between shell i and the previous layer).
        // No volume normalization as this is done by the caller.
        contrast = sld[i] - last_sld;
        vol = M_4PI_3 * cube(r);
        f += vol * sas_3j1x_x(q*r) * contrast;

        last_sld = sld[i];
        r += thickness[i];
    }

    // Contrast between the outermost shell and the solvent
    contrast = sld_solvent - last_sld;
    vol = M_4PI_3 * cube(r);
    f += vol * sas_3j1x_x(q*r) * contrast;

    //rescale to [cm-1].
    // With volume fraction as part of the model in the dilute limit need
    // to return F2 = Vf <fq^2>.  In order for beta approx. to work correctly
    // need F1^2/F2 equal to <fq>^2 / <fq^2>.  By returning F1 = sqrt(Vf) <fq>
    // and F2 = Vf <fq^2> both conditions are satisfied.
    // Since Vf is the volume fraction of vesicles of all radii, it is
    // constant when averaging F1 and F2 over radii and so pops out of the
    // polydispersity loop, so it is safe to apply it inside the model
    // (albeit conceptually ugly).
    *F1 = 1e-2 * sqrt(volfraction) * f;
    *F2 = 1.0e-4 * volfraction * f * f;
}
