
static double
f_constant(double q, double r, double sld)
{
  const double bes = sas_3j1x_x(q * r);
  const double vol = M_4PI_3 * cube(r);
  return sld * vol * bes;
}

static double
form_volume(double core_radius, double fp_n, double thickness[])
{
  double r = core_radius;
  int n = (int)(fp_n+0.5);
  for (int i=0; i < n; i++) {
    r += thickness[i];
  }
  return M_4PI_3 * cube(r);
}

static double
Iq(double q, double core_sld, double core_radius,
   double solvent_sld, double fp_n, double sld[], double thickness[])
{
  const int n = (int)(fp_n+0.5);
  double f, r, last_sld;
  r = core_radius;
  last_sld = core_sld;
  f = 0.;
  for (int i=0; i<n; i++) {
    f += M_4PI_3 * cube(r) * (sld[i] - last_sld) * sas_3j1x_x(q*r);
    last_sld = sld[i];
    r += thickness[i];
  }
  f += M_4PI_3 * cube(r) * (solvent_sld - last_sld) * sas_3j1x_x(q*r);
  return f * f * 1.0e-4;
}
